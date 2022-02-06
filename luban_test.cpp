//
// Created by 张力文 on 2022/1/13.
//

/*
 * 输入1：运行脚本run.inp （先默认手动输入，因为涉及的参数并不多）
 * 输入2：构型文件cnf.inp，因为含有多个原子（256个）的构型坐标，需要重点处理
 * 原始的输入格式需要改进，如果要求输入.xyz格式文件，缺少第一列原子类型
 * 或者，根据输入的原子个数来创建对应的二维数组（n * 3）维，然后分别赋值（最简单）
 *
 */

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include<fstream>
#include<math.h>
#include<random>
#include<algorithm>


//输入input.txt文件，第一行为原子个数n，后面n行为每个原子的坐标
using namespace std;

vector<double> read_cnf_atoms()
{
    ifstream myfile("/Users/zlw/ClionProjects/MCCPS Luban/input.txt");
    //string temp;
    //string para[5]={"ensemble","temperature", "nmolty", "nummols", "numboxes"};
    //string ensemble;
    //double temperature;
    int nmolty;
    //nummols, numboxes;


    if(!myfile.is_open())
    {
        cout<<"file open error"<<endl;
    }

    //myfile>>nmolty;


    //cout<<nmolty<<endl; for test
/*
    double **array;

    for(int i = 0; i < nmolty; i++)
    {
        array[i] = new double[4];
    }

    for(int i = 0; i < nmolty; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            myfile>>array[i][j];
            cout<<array[i][j];
        }
        cout<<endl;
    }
    */
    //int i = 0;//全局变量要小心
    vector<double> V;
    //vector<double>::iterator it;

    double d;
    while(myfile >> d)
    {
        V.push_back(d);
    }

    myfile.close();

    //double aa[100][4];
    /*
    it = V.begin();
    nmolty = (int)*it;
    it++;

    double** aa = new double*[nmolty];
    for(int i = 0; i < nmolty; i++)
    {
        aa[i] = new double[4];
    }

    for(it; it != V.end();)
    {
        for(int i = 0; i < nmolty; i++)
        {
            for (int j = 0; j < 4; j++)
            {

                aa[i][j] = *it;
                it++;
                //cout<<' '<<aa[i][j];
            }

            //cout<<endl;
        }

    }
*/





    return V;
}

double potential_1(double *ri, double box_length, double r_cut, double **atoms, int n)
{
    double sr2_ovr = 1.77;
    double r_cut_box = r_cut / box_length;
    double r_cut_box_sq = pow(r_cut_box,2);
    double box_sq = pow(box_length,2);

    double potential = 0.0;
    double rij_sq;
    double sr2;
    double sr6;
    double sr12;
    bool ovr;
    //int n = (sizeof(atoms)/sizeof(atoms[0][0]));
    //int n = 10;
    //cout<<"potential_1中的n："<<n<<endl;
    double partial = 0.0;
    double* rij = new double[3];
    //针对盒子中的某个原子，当前的atoms_sub是除了该原子外的所有原子坐标二维数组
    for(int j = 0; j < n; j++)
    {
        if((ri[0] == atoms[j][0])&&(ri[1] == atoms[j][1])&&(ri[2] == atoms[j][2]))
        {
            continue;
        }

        //cout<<"当前的j值："<<j<<endl;
        rij[0] = ri[0] - atoms[j][0];
        rij[1] = ri[1] - atoms[j][1];
        rij[2] = ri[2] - atoms[j][2];
        //cout<<rij[0]<<' '<<rij[1]<<' '<<rij[2]<<endl;
        //在模拟盒尺寸为1时，应用周期性边界条件
        rij[0] = rij[0] - round(rij[0]);
        rij[1] = rij[1] - round(rij[1]);
        rij[2] = rij[2] - round(rij[2]);
        //cout<<rij[0]<<' '<<rij[1]<<' '<<rij[2]<<endl;
        rij[0] = pow(rij[0],2);
        rij[1] = pow(rij[1],2);
        rij[2] = pow(rij[2],2);
        //cout<<rij[0]<<' '<<rij[1]<<' '<<rij[2]<<endl;
        rij_sq = rij[0] + rij[1] + rij[2];

        if(rij_sq < r_cut_box_sq)
        {
            rij_sq = rij_sq * box_sq;
            //cout<<"rij_sp: "<<rij_sq<<endl;
            sr2 = 1.0 / rij_sq;
            //cout<<"sr2: "<<sr2<<endl;
            /*ovr = sr2 > sr2_ovr;
            if(ovr)
            {

            }
             */
            sr6 = pow(sr2, 3);
            //cout<<"sr6: "<<sr6<<endl;
            sr12 = pow(sr6, 2);
            //cout<<"sr12: "<<sr12<<endl;
            potential = sr12 - sr6;
            partial = partial + potential;
            //cout<<"potential_1中的partial： "<<partial<<endl;
        }


    }
    partial = 4.0 * partial;

    return partial;

}


double potential(double box_length, double r_cut, double **atoms,int n)
{
    double total = 0.0;
    double partial;
    //int n = (sizeof(atoms)/sizeof(atoms[0][0]));
    //int n = 10;
    //double **temp = new double*[n];
    //double **atoms_1 = new double*[n];
    /*for(int i = 0; i < n; i++)
    {
        temp[i] = new double[3];
    }
    for(int i = 0; i < n - 1; i++)
    {
        atoms_1[i] = new double[3];
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            temp[i][j] = atoms[i][j];
            //cout<<temp[i][j]<<endl;
        }
        //cout<<endl;
    }
    */
    //cout<<"potential中的n："<<n<<endl;

    for(int i = 0; i < n; i++)
    {
        //这里当前的atoms[i]要从总atoms里删除，具体操作：先整体拷贝到另外一个数组，再在另外一个数组中删除
       /*
        if(i == 0)
        {
            for(int m = i + 1; m < n; m++)
            {
                for(int n = 0; n < 3; n++)
                {
                    atoms_1[m][n] = temp[m - 1][n];
                    cout<<atoms_1[m][n];
                }
                cout<<endl;
            }


        }
        else
        {
            //int m = i - 1; i=3
            for(int m = 0; m < i; m++)
            {
                for(int n = 0; n < 3; n++)
                {
                    atoms_1[m][n] = temp[m][n];
                }
            }


            for(int m = i; m < n - 1; m++)
            {
                for(int n = 0; n < 3; n++)
                {
                    atoms_1[m][n] = temp[m + 1][n];
                }
            }
        }
*/
        partial = potential_1(atoms[i],box_length,r_cut,atoms,n);
        //cout<<"partial:"<<partial<<endl;
        total = total + partial;
        //cout<<"total:"<<total<<endl;
    }

   /* for(int i = 0; i < n; i++)
    {
        delete[] temp[i];
    }
    delete[] temp;
    for(int i = 0; i < n; i++)
    {
        delete[] atoms_1[i];
    }
    delete[] atoms_1;
    */
    return total;
}

double* random_translate_vector(double dr_max, double *old)
{
    random_device seed_device;
    default_random_engine  engine;
    engine.seed(seed_device());
    uniform_real_distribution<double> distr(0, 1);
    double* zeta = new double[3];
    for(int i = 0; i < 3; i++)
    {
        zeta[i] = distr(engine);
        zeta[i] = 2.0 * zeta[i] - 1.0;
        //cout<<"zeta[i]:"<<zeta[i]<<endl;
        old[i] = old[i] + zeta[i]*dr_max;
        //cout<<"old[i]:"<<old[i]<<endl;
    }

    return old;

}

bool metropolis(double delta)
{
    double exponent_guard = 0.75;
    if(delta > exponent_guard)
        return false;
    else if(delta <0.0)
        return true;
    else
    {
        double zeta = 0.0;
        random_device seed_device;
        default_random_engine  engine;
        engine.seed(seed_device());
        uniform_real_distribution<double> distr(0, 1);
        zeta = distr(engine);
        return (exp(-delta)> zeta);
    }
}

int main()
{

    vector<double> v1 = read_cnf_atoms();
    vector<double>::iterator it;

    int nmolty;
    int nstep = 1;
    double box_length = 6.98864372;
    double r_cut = 2.0;
    double dr_max = 0.15;

    it = v1.begin();
    //cout<<"nmolty: "<<*it<<endl;

    nmolty = *it;
    //cout<<"模拟原子个数： "<<nmolty<<endl;
    it++;

    double** atoms = new double*[nmolty];
    double** cors = new double* [nmolty];
    for(int i = 0; i < nmolty; i++)
    {
        atoms[i] = new double[4];
        cors[i] = new double[3];
    }

    for(it; it != v1.end();)
    {
        for(int i = 0; i < nmolty; i++)
        {
            for (int j = 0; j < 4; j++)
            {

                atoms[i][j] = *it;
                it++;
                //cout<<' '<<atoms[i][j];
            }

            //cout<<endl;
        }

    }

    for(int i = 0; i < nmolty; i++)
    {
        for(int j = 0; j < 3; j++)
        {

            cors[i][j] = atoms[i][j+1];
            cout<<cors[i][j]<<' ';
        }
        cout<<endl;
    }


    double total1 = potential(box_length,r_cut,cors,nmolty);

    double partial_old = 0.0;
    double partial_new = 0.0;
    double delta = 0.0;
    double *ri = new double[3];

    cout<<"total1:"<<total1<<endl;

    int m_acc = 0;
    double m1_ratio = 0.0;
    for(int i =0; i < nstep; i++)
    {
        for(int j = 0; j < nmolty; j++)
        {
            partial_old = potential_1(cors[i],box_length,r_cut,cors,nmolty);
            ri = random_translate_vector(dr_max/box_length,cors[i]);
            ri[0] = ri[0] - round(ri[0]);//周期性边界条件
            ri[1] = ri[1] - round(ri[1]);
            ri[2] = ri[2] - round(ri[2]);
            cout<<ri[0]<<" "<<ri[1]<<" "<<ri[2]<<endl;


            partial_new = potential_1(ri,box_length,r_cut,cors,nmolty);
            cout<<"partial_old:"<<partial_old<<endl;
            cout<<"partial_new:"<<partial_new<<endl;
            delta = partial_old - partial_new;
            //cout<<"delta:"<<delta<<endl;
            if(metropolis(delta))
            {
                //更新能量
                total1 = total1 +partial_new - partial_old;
                //更新坐标
                ri[0] = cors[i][0];
                ri[1] = cors[i][1];
                ri[2] = cors[i][2];

                m_acc = m_acc + 1;
                cout<<"m_acc:"<<m_acc<<endl;
            }


        }

        m1_ratio = m_acc / nmolty;
        cout<<"m1_ratio:"<<m1_ratio<<endl;

    }







    for(int i = 0; i < nmolty; i++)
    {
        delete[] atoms[i];
    }
    delete[] atoms;

    return 0;
}