#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include<math.h>
#include<random>

//#include "config_io.h"

using namespace std;

/*
void read_cnf_atoms()
{
    ifstream myfile("/Users/zlw/ClionProjects/MCCPS Luban/input.txt");
    string temp;
    string para[5]={"ensemble","temperature", "nmolty", "nummols", "numboxes"};
    string ensemble;
    double temperature;
    int nmolty, nummols, numboxes;

    if(!myfile.is_open())
    {
        cout<<"file open error"<<endl;
    }

    int i = 0;//全局变量要小心


    for(string temp; getline(myfile,temp);)
    {
        //outfile<<temp<<endl;
        //由于整行整行读入，读入道str中时，文件中的每个换行符都被丢掉了，为了照原样输出，在out流上输出时需要再补上一个回车<<endl
        int res = temp.compare(para[i++]);
        //下面三个输出用来测试

        /*
        cout<<temp<<endl;//temperature
        cout<<para[i]<<endl;//ensemble
        cout<<"res = "<<res<<endl;
        */
/*
        if(res == 0)
        {
            cout<<temp<<":"<<endl;

            getline(myfile,temp);

            //将字符串参数传回main函数


            cout<<temp<<endl;
        }
        else
        {
            cout<<"parameter error"<<endl;
        }

        //如果是可选参数，不能傻乎乎的比对字符串读取下一行，需要判断有这个可选参数对号入座读下一行，如果没有这个参数就跳过


    }

    myfile.close();


    return;
}

*/
//为了测试给出二维数组形式的原子坐标，全局变量
double r1[][3] = {{1.0,1.3,4.3},{3.7,2.4,8.6},{2.9,1.0,0.4},{2.0,3.0,4.0},{8.0,9.0,5.0}};
double r2[][3] = {{3.0,1.2,5.6}, {5.2,1.5,3.0}, {4.3,6.3,8.9}, {1.2, 4.0, 3.0}, {1.0,3.7,2.9}};

class VariableType{
    public:
    string name;
    string val;
};

class PotentialType{
    public:
    double pot;
    double vir;

    PotentialType potential_l(int i, double ri[i][3], int box_length, double r_cut, double r[][3]);
    PotentialType potential(double box_length, double r_cut, double r[][3]);

};

PotentialType PotentialType::potential_l(int i, double ri[i][3], int box_length, double r_cut, double r[][3])
{
    //ri是一个原子的三维坐标（x，y，z）
    double sr2_ovr = 1.77;
    double r_cut_box = r_cut / box_length;
    double r_cut_box_sq = pow(r_cut_box, 2);
    double box_sq = pow(box_length, 2);

    double rj;

    PotentialType partial;
    partial.pot = 0.0;
    partial.vir = 0.0;
/*
    int count = sizeof(r) / sizeof(r[0]);
    for(int j = 0; j < count; j++)
    {
        rj = ri[j] - r[j];

    }

 */
    random_device seed_device;
    default_random_engine  engine;
    engine.seed(seed_device());
    uniform_int_distribution<int> distr(-10, 10);
    partial.pot = distr(engine);
    partial.vir = distr(engine);

    return partial;
}

PotentialType PotentialType::potential(double box_length, double r_cut, double r[][3])
{
    PotentialType total;
    PotentialType partial;
    total.pot = 0.0;
    total.vir = 0.0;
    int count = sizeof(r) / sizeof(r[0]);
    for(int i = 0; i < count; i++)
    {
      partial = potential_l(r[i][3], box_length, r_cut, r);
    }
    total.pot = total.pot + partial.pot;
    total.vir = total.vir + partial.vir;

    return total;
}

void introductioin()
{
    printf("Lennard-Jones potential\n");
    printf("Cut (but not shifted)\n");
    printf("Diameter, Sigma = 1\n");
    printf("Well depth, epsilon = 1\n");


}
int main() {


   // read_cnf_atoms();
   /*for(int i = 0; i < 5; i++)
   {
       cout<<para_inp[i]<<endl;
   }
    */

    string ensemble = "gemc_nvt";
    double temperature = 1.0;
    int n1 = 256; //原子个数
    int n2 = 256;
    double box_length1 = 5.0;
    double box_length2 = 5.0;
    double r_cut = 2.0;

    //从输入数据（原子个数、box长度）计算体积和密度
    double vol1 = pow(box_length1,3);
    double vol2 = pow(box_length2,3);
    double density1 = n1 / vol1;
    double density2 = n2 / vol2;

    //从输入脚本中读入的模拟相关值
    int nblock = 10;
    int nstep = 100;
    int nswap = 20;
    double dr_max = 0.15;
    double dv_max = 10.0;


    //移动、交换、体积变化接受概率
    VariableType move1_r;
    VariableType move2_r;
    VariableType swap12_r;
    VariableType swap21_r;
    VariableType n_1;
    VariableType n_2;
    VariableType density_1;
    VariableType density_2;
    VariableType energy_1;
    VariableType energy_2;
    VariableType pressure_1;
    VariableType pressure_2;

    PotentialType total1;
    PotentialType total2;




    move1_r.name = "Move ratio(1)";
    move1_r.val = 0.0;
    move2_r.name = "Move ratio(2)";
    move2_r.val = 0.0;

    swap12_r.name = "Swap ratio(1->2)";
    swap12_r.val = 0.0;
    swap21_r.name = "Swap ratio(2->1)";
    swap21_r.val = 0.0;

    n_1.name = "Number (1)";
    n_1.val = (float)n1;
    n_2.name = "Number (2)";
    n_2.val = (float)n2;

    density_1.name = "Density (1)";
    density_1.val = density1;
    density_2.name = "Density (2)";
    density_2.val = density2;

    PotentialType total1_energy = total1.potential(box_length1, r_cut, r1);
    PotentialType total2_energy = total2.potential(box_length2,r_cut,r2);
    energy_1.name = "E/N cut (1)";
    energy_1.val = 1.5 * temperature + total1_energy.pot/n1;
    energy_2.name = "E/N cut (2)";
    energy_2.val = 1.5 * temperature + total2_energy.pot/n2;

    pressure_1.name = "P cut (1)";
    pressure_1.val = 1.5 * temperature + total1_energy.pot/n1;
    pressure_2.name = "P cut (2)";
    pressure_2.val = 1.5 * temperature + total2_energy.pot/n2;

    printf("Gibbs Ensemble Monte Carlo program\n");
    printf("Simulation uses cut(but not shifted) potential\n");

    introductioin();

    printf("Number of blocks:   %d\n", &nblock);





















    return 0;
}
