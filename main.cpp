#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include<math.h>
#include<random>
#include<algorithm>

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
//static double r1[][3] = {{-3.0575316266, -3.0575316266, -3.0575316266},{-0.4367902324, -3.0575316266, 1.3103706971},{ -1.3103706971, 2.1839511618, -1.3103706971},{0.4367902324, -1.3103706971, -1.3103706971},{ 2.1839511618, -1.3103706971, -3.0575316266}};
//static double r2[][3] = {{3.0575316266, 0.4367902324, 1.3103706971}, {3.0575316266, -0.4367902324, -3.0575316266}, {1.3103706971, 3.0575316266, 2.1839511618}, {-2.1839511618, -0.4367902324, -1.3103706971}, {-2.1839511618, -2.1839511618, -1.3103706971}};

//blk_begin中声明的global变量
double blk_nrm, blk_avg, blk_msd;

//run_begin中声明的global变量
int n_avg, line_width;
double run_nrm, run_avg, run_err;



class VariableType{
    public:
    string name;
    string val;
};

class PotentialType{
    public:
    double pot;
    double vir;

    PotentialType potential_l(int i, double r[i][3], int box_length, double r_cut);
    PotentialType potential(double box_length, double r_cut, double r[][3]);
    //friend PotentialType operator+(PotentialType& pot, PotentialType& vir);

};

PotentialType PotentialType::potential_l(int i, double r[i][3], int box_length, double r_cut)
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
      partial = potential_l(i,r, box_length, r_cut);
    }
    total.pot = total.pot + partial.pot;
    total.vir = total.vir + partial.vir;

    return total;
}

PotentialType operator+(const PotentialType& a, const PotentialType& b)
{
    //加friend关键字的原因：operator+函数都要用Potential对象进行私有数据访问，作为普通函数这是不允许的，即编译器会报错。
    //加入friend关键字声明后，编译就对或者两个函数防伪Potential的私有数据网开一面了。
    PotentialType s;
    s.pot = a.pot + b.pot;
    s.vir = a.vir + b.vir;

    return s;
}

PotentialType operator-(const PotentialType& a, const PotentialType& b)
{
    //加friend关键字的原因：operator+函数都要用Potential对象进行私有数据访问，作为普通函数这是不允许的，即编译器会报错。
    //加入friend关键字声明后，编译就对或者两个函数防伪Potential的私有数据网开一面了。
    PotentialType s;
    s.pot = a.pot - b.pot;
    s.vir = a.vir - b.vir;

    return s;
}

void introductioin()
{
    printf("Lennard-Jones potential\n");
    printf("Cut (but not shifted)\n");
    printf("Diameter, Sigma = 1\n");
    printf("Well depth, epsilon = 1\n");


}

void run_begin()
{
    //基于提供的列表名字建立平均值以及其他属性
    //该函数关键点在于画出了表格，表头、数据等，且是动态的



}

void blk_begin()
{
    blk_nrm = 0.0;
    blk_avg = 0.0;
}

double* random_translate_vector(double dr_max, double r[][3])
{
    //生成3个0-1的随机浮点数
    double num;
    double *zeta = (double *)calloc(3, sizeof(double *));
    for(int i = 0; i < 3; i++)
    {
        srand(time(NULL));
        num = (rand() % 1000) * 0.001;
        zeta[i] = num;
        zeta[i] = 2.0 * zeta[i] - 1.0;
    }


    //返回一个二维数组，选中的哪一行每一个分量加上zeta的值
    return zeta;
}

bool metropolis(double delta)
{
    double exponent_guard = 75.0;
    double zeta;
    if(delta > exponent_guard)
        return false;
    else
    {
        srand(time(NULL));
        zeta = (rand() % 1000) * 0.001;
        return exp(-delta) > zeta;
    }
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
    int n1 = 5; //原子个数
    int n2 = 5;
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

    //写出相应参数
    printf("Number of blocks:   %d\n", &nblock);
    printf("Number of steps per block:  %d\n", &nstep);
    printf("Swap attempts per step:     %d\n", &nswap);
    printf("Specified temperature:      %d\n", &temperature);
    printf("Potential cutoff distance:      %lf\n", &r_cut);
    printf("Maximum displacement:   %lf\n", &dr_max);
    printf("Maximum volume change:      %d\n", &dv_max);

    printf("Number of particles:    %d      %d\n", &n1, &n2);
    printf("Simulation box length:      %d      %d\n", &box_length1, &box_length2);
    printf("Density:    %lf      %d", &density1, &density2);


    double m1_ratio = 0.0;
    run_begin();


    //从构型文件中读入原子坐标

    //将坐标转换为模拟盒单位
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            r1[i][j] = r1[i][j] / box_length1;
            r2[i][j] = r2[i][j] / box_length2;
        }
    }

    for(int i = 0; i < nblock; i++)
    {
        //blk_begin


        for(int j = 0; j < nstep; j++)
        {
            double m_acc = 0.0;
            double *zeta = (double *) calloc(3,sizeof(double));

            for(int k = 0; k < n1; k++)
            {
                PotentialType partial_old;
                PotentialType partial_new;
                double delta;
                partial_old.potential_l(k,r1,box_length1,r_cut);
                zeta = random_translate_vector(dr_max / box_length1, r1);
                for (int i = 0; i < 3; ++i) {
                    r1[k][i] = r1[k][i] + zeta[i];
                }

                partial_new.potential_l(k,r1,box_length1,r_cut);
                delta = partial_new.pot - partial_old.pot;
                delta = delta / temperature;
                if(metropolis(delta))
                {
                    total1_energy = operator+(total1_energy,(operator-(partial_new,partial_old))) ;
                    m_acc++;
                }

            m1_ratio / m_acc / n1;



            }
        }
    }


























    return 0;
}
