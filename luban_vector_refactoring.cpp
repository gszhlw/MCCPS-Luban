//
// Created by 张力文 on 2022/2/19.
//

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include<fstream>
#include<math.h>
#include<random>
#include<algorithm>
#include<iomanip>
#include<ctime>


//输入input.txt文件，第一行为原子个数n，后面n行为每个原子的坐标
using namespace std;

class PotentialType{
public:
    double pot;
    bool ovr = false;

    int nmolty;

    //PotentialType(double pot, bool ovr);
    //static PotentialType potential_1(vector<double> *ri, double box_length, double r_cut, vector<vector<double>> *atoms, int n);
    static PotentialType potential_1(double *ri, double box_length, double r_cut, double **atoms, int n);
    //PotentialType operator+(PotentialType &p) const;
    //PotentialType operator*(PotentialType &p1,PotentialType &p2) const;
    PotentialType potential(double box_length, double r_cut, double r[][3]);
    //friend PotentialType operator+(PotentialType& pot, PotentialType& vir);

};

PotentialType PotentialType :: potential_1(double *ri, double box_length, double r_cut, double **atoms, int n)
{
   // partial = PotentialType::potential_1(atoms[i],box_length,r_cut,atoms,n);
    double sr2_ovr = 1.77;
    double r_cut_box = r_cut / box_length;
    double r_cut_box_sq = pow(r_cut_box,2);
    double box_sq = pow(box_length,2);

    double rij_sq;
    double sr2;
    double sr6;
    double sr12;
    bool ovr;

    PotentialType partial;
    PotentialType potential;

    double* rij = new double[3];


    for(int j = 0; j < n; j++)
    {
        if((ri[0] == atoms[j][0])&&(ri[1] == atoms[j][1])&&(ri[2] == atoms[j][2]))
        {
            continue;
        }

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
        //cout<<setprecision(8)<<rij[0]<<' '<<rij[1]<<' '<<rij[2]<<endl;

        rij_sq = rij[0] + rij[1] + rij[2];

        if(rij_sq < r_cut_box_sq)
        {
            rij_sq = rij_sq * box_sq;
            //cout<<"rij_sp: "<<rij_sq<<endl;
            sr2 = 1.0 / rij_sq;
            //cout<<"sr2: "<<sr2<<endl;
            ovr = (sr2 > sr2_ovr);
            if(ovr)
            {
                potential.ovr = true;
                return partial;
            }

            sr6 = pow(sr2, 3);
            //cout<<"sr6: "<<sr6<<endl;
            sr12 = pow(sr6, 2);
            //cout<<"sr12: "<<sr12<<endl;
            //potential = sr12 - sr6;
            potential.pot= sr12 - sr6;
            partial.pot = partial.pot + potential.pot;
            //cout<<"potential_1中的partial： "<<partial.pot<<endl;
        }
    }

    partial.pot = 4.0 * partial.pot;

    return partial;

}

PotentialType potential(double box_length, double r_cut, double **atoms,int n)
{

    PotentialType total;
    total.pot = 0.0;
    PotentialType partial;
    partial.pot = 0.0;

    //for debug
    //int count = 1;
    for(int i = 0; i < n; i++)
    {
        partial = PotentialType::potential_1(atoms[i],box_length,r_cut,atoms,n);
        //重叠判断
        if(partial.ovr)
        {
            total.ovr = true;
            break;
        }
        //cout<<"partial:"<<partial.pot<<endl;
        total.pot = total.pot + partial.pot;
        //cout<<count++<<" "<<total.pot<<endl;
    }

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


vector<double> read_cnf_atoms(string filename)
{
    ifstream myfile("/Users/zlw/ClionProjects/MCCPS Luban/"+filename);


    if(!myfile.is_open())
    {
        cout<<"file open error"<<endl;
    }


    vector<double> V;

    double d;
    while(myfile >> d)
    {
        V.push_back(d);
    }

    myfile.close();

    return V;
}


double **convertVector2array(vector<vector<double>> v)
{
    int rows = v.size();
    int cols = v[0].size();
    double **outdata_array;
    outdata_array = new double *[rows];
    for(int i = 0; i < rows; i++)
    {
        outdata_array[i] = new double[cols];
    }

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            outdata_array[i][j] = v[i][j];
        }
    }

    return outdata_array;
}

//释放数组内存函数
void free_2D_array(double **array, int rows)
{
    for(int i = 0; i < rows; i++)
    {
        delete[] array[i];
    }

    delete[] array;
}



int main() {


    //模拟参数
    double box_length = 6.98864372;
    double r_cut = 2.0;
    double dr_max = 0.15;
    double temperature = 1.0;
    int nstep = 1;
    int nswap = 5;


    vector<double> v1 = read_cnf_atoms("input_box1.txt");
    vector<double> v2 = read_cnf_atoms("input_box2.txt");

    vector<double>::iterator it;

    int nmolty1;
    int nmolty2;

    it = v1.begin();
    //cout<<"nmolty: "<<*it<<endl;

    nmolty1 = *it;
    cout<<"box1模拟原子个数： "<<nmolty1<<endl;

    it++;

    vector<vector<double> > atoms1(nmolty1, vector<double>(4));
    vector<vector<double> > cors1(nmolty1, vector<double>(3));
    //vector<vector<double>>::iterator iter;
    //vector<double>::iterator it_col;
    //vector<double>vec_tmp;
/*
    for (int iter = atoms1.begin(); iter != atoms1.end() ; ++iter) {
        for(it_col = atoms_col.begin(), it_col != atoms1_col.end(),it_col++)
        {
             = *it;
            it++;
        }

    }
    */
    for(int i =0; i < nmolty1;i++)
    {
        for(int j = 0; j < atoms1[0].size(); j++)
        {
            atoms1[i][j]= *it++;
            //cors1.push_back(atoms1[i][j+1]) ;
            //cout<<atoms1[i][j]<<' ';

        }
        //cout<<endl;
    }

    for(int i =0; i < nmolty1;i++)
    {
        for(int j = 0; j < cors1[0].size(); j++)
        {
            //atoms1[i][j]= *it++;
            //cors1.push_back(atoms1[i][j+1]) ;
            cors1[i][j] = atoms1[i][j+1];
            //cout<<setprecision(9)<<cors1[i][j]<<' ';

        }
        //cout<<endl;
    }

    double **atoms1_array = convertVector2array(atoms1);
    double **cors1_array = convertVector2array(cors1);


/*测试用
    for(int i = 0; i < nmolty1; i++)
    {
        for(int j = 0; j < 3; j++)
        {

            //cors1[i][j] = atoms1[i][j+1];
            cout<<cors1_array[i][j]<<' ';
        }
        cout<<endl;
    }

*/

    //box2
    it = v2.begin();
    //cout<<"nmolty: "<<*it<<endl;

    nmolty2 = *it;
    cout<<"box2模拟原子个数： "<<nmolty2<<endl;

    it++;

    vector<vector<double> > atoms2(nmolty1, vector<double>(4));
    vector<vector<double> > cors2(nmolty1, vector<double>(3));

    for(int i =0; i < nmolty2;i++)
    {
        for(int j = 0; j < atoms2[0].size(); j++)
        {
            atoms2[i][j]= *it++;
            //cors1.push_back(atoms1[i][j+1]) ;
            //cout<<atoms1[i][j]<<' ';

        }
        //cout<<endl;
    }

    for(int i =0; i < nmolty2;i++)
    {
        for(int j = 0; j < cors2[0].size(); j++)
        {
            //atoms1[i][j]= *it++;
            //cors1.push_back(atoms1[i][j+1]) ;
            cors2[i][j] = atoms2[i][j+1];
            //cout<<cors1[i][j]<<' ';

        }
        //cout<<endl;
    }

    double **atoms2_array = convertVector2array(atoms1);
    double **cors2_array = convertVector2array(cors1);



    PotentialType total1= potential(box_length,r_cut,cors1_array,nmolty1);
    //重叠判断
    if(total1.ovr) throw "Overlap in initial configuration 1";

    PotentialType total2 = potential(box_length,r_cut,cors2_array,nmolty2);
    if(total2.ovr) throw "Overlap in initial configuration 2";

    PotentialType partial_old;
    PotentialType partial_new;

    double delta = 0.0;
    double *ri = new double[3];

    cout<<"total1:"<<total1.pot<<endl;
    cout<<"total2:"<<total2.pot<<endl;

    int m_acc = 0;
    double m1_ratio = 0.0;
    double m2_ratio = 0.0;
    double x12_ratio = 0.0;
    double x21_ratio = 0.0;

    for(int i =0; i < nstep; i++) {
        //box1
        for (int j = 0; j < nmolty1; j++) {
            partial_old = PotentialType::potential_1(cors1_array[j], box_length, r_cut, cors1_array, nmolty1);
            ri = random_translate_vector(dr_max / box_length, cors1_array[j]);
            ri[0] = ri[0] - round(ri[0]);//周期性边界条件
            ri[1] = ri[1] - round(ri[1]);
            ri[2] = ri[2] - round(ri[2]);

            partial_new = PotentialType::potential_1(ri, box_length, r_cut, cors1_array, nmolty1);
            //cout<<"box1：partial_old:"<<partial_old.pot<<endl;
            //cout<<"box1：partial_new:"<<partial_new.pot<<endl;
            //cout << setprecision(8) << ri[0] << " " << ri[1] << " " << ri[2] << endl;
            if (!partial_new.ovr) {
                delta = partial_old.pot - partial_new.pot;
                delta = delta / temperature;
            }

            if (metropolis(delta)) {
                //更新能量
                total1.pot = total1.pot + partial_new.pot - partial_old.pot;

                //ri[0] = cors1[i][0];
                //ri[1] = cors1[i][1];
                //ri[2] = cors1[i][2];

                //更新坐标
                cors1_array[i][0] = ri[0];
                cors1_array[i][1] = ri[1];
                cors1_array[i][2] = ri[2];

                m_acc = m_acc + 1;
                cout << "box1：m_acc:" << m_acc << endl;
            }


        }

        m1_ratio = (double) m_acc / nmolty1;
        cout << "m1_ratio:" << m1_ratio << endl;

        m_acc = 0;

        //box2
        for (int j = 0; j < nmolty2; j++) {
            partial_old = PotentialType::potential_1(cors2_array[j], box_length, r_cut, cors2_array, nmolty2);
            ri = random_translate_vector(dr_max / box_length, cors2_array[j]);
            ri[0] = ri[0] - round(ri[0]);//周期性边界条件
            ri[1] = ri[1] - round(ri[1]);
            ri[2] = ri[2] - round(ri[2]);

            partial_new = PotentialType::potential_1(ri, box_length, r_cut, cors2_array, nmolty2);
            //cout<<"box2：partial_old:"<<partial_old.pot<<endl;
            //cout<<"box2：partial_new:"<<partial_new.pot<<endl;

            if (!partial_new.ovr) {
                delta = partial_old.pot - partial_new.pot;
                delta = delta / temperature;
            }

            if (metropolis(delta)) {
                //更新能量
                total2.pot = total2.pot + partial_new.pot - partial_old.pot;
                //更新坐标
                //ri[0] = cors2_array[i][0];
                //ri[1] = cors2_array[i][1];
                //ri[2] = cors2_array[i][2];
                cors2_array[i][0] = ri[0];
                cors2_array[i][1] = ri[1];
                cors2_array[i][2] = ri[2];


                m_acc = m_acc + 1;
                cout << "box2：m_acc:" << m_acc << endl;
            }

        }

        m2_ratio = (double) m_acc / nmolty2;
        cout <<setprecision(4)<< "m2_ratio:" << m2_ratio << endl;

        //swap move

        int x12_try = 0, x21_try = 0;
        int x12_acc = 0, x21_acc = 0;

        for(int j = 0; j < nswap; j++)
        {
            //3个（0，1）均匀分布的随机数
            PotentialType  partial_old;
            PotentialType partial_new;


            //这里的逻辑有点问题，如果是随机生成的一组坐标，不能确定他就是盒子1中的某个粒子

            double* ri = new double[3];

            //随机选择box1中的任意一个粒子

            random_device seed_device;
            default_random_engine  engine;
            uniform_real_distribution<double> distr(0, 1);
           /* engine.seed(seed_device());
            uniform_int_distribution<int> distr(0, nmolty1);
            int i = distr(engine);
             */
            for(int i = 0; i < 3; i++)
            {
                //ri为-0.5～0.5之间的数，3个为一组坐标,为欲交换到的目标坐标

                ri[i] = distr(engine);
                ri[i] = ri[i] - 0.5;


            }
            /*
            for(int j = 0; j < 3; j++)
            {
                ri[j] = cors1_array[i][j];
            }
            */
            //random_device seed_device;
            //default_random_engine  engine;
            engine.seed(seed_device());
            uniform_real_distribution<double> distr_rand(0, 1);
            double rand_num = distr_rand(engine);
           // double rand = 0.4;
            cout<<"rand:"<<rand_num<<endl;


            //swapping1->2
            double *curr = NULL;
            if(rand_num < 0.5)
            {
                x12_try = x12_try + 1;
                if(nmolty1 > 1)
                {

                    random_device seed_device;
                    default_random_engine  engine;
                    engine.seed(seed_device());
                    uniform_int_distribution<int> distr(0, nmolty1);


                    /*
                    int i;
                    srand((int)time(0));
                    i = rand()%nmolty1;
                     */
                    //随机选择box1中的任意一个粒子
                    int i = distr(engine);//i就是box1中要删去的第i行，坐标

                    //box1中需要删去的坐标必须是已经在box1中的坐标，但是新插入box2的坐标是任意的（因为插入任意位置）

                    //debug
                    cout<<"cors1_array:"<<cors1_array[i][0]<<" "<<cors1_array[i][1]<<" "<<cors1_array[i][2]<<endl;
                    cout<<"ri:"<<ri[0]<<" "<<ri[1]<<" "<<ri[2]<<endl;
                    cout<<"i:"<<i<<endl;


                    partial_old = PotentialType::potential_1(cors1_array[i],box_length,r_cut,cors1_array,nmolty1);
                    //if(!partial_old.ovr) throw "Overlap on particle removal";

                    partial_new = PotentialType::potential_1(ri,box_length,r_cut,cors2_array,nmolty2);
                    //cout<<"work"<<endl;
                    cout<<"box1：partial_old:"<<partial_old.pot<<endl;
                    cout<<"box1：partial_new:"<<partial_new.pot<<endl;

                    if(!partial_new.ovr)
                    {
                        delta = (partial_new.pot - partial_old.pot) / temperature;
                        delta = delta - log(pow(box_length,3)/(nmolty2 + 1));
                        delta = delta + log(pow(box_length,3)/nmolty1);
                        cout<<"delta:"<<delta<<endl;
                        if(metropolis(delta))
                        {
                            //ri的三个坐标用push_back推入二维向量vector cors2中
                            cors2[nmolty2].push_back(ri[0]);
                            cors2[nmolty2].push_back(ri[1]);
                            cors2[nmolty2].push_back(ri[2]);



                            //先释放原来的cors2_array的内存，在把新的cors2 vector转换成增加了一行数据后的数组
                            //先记录原来的二维数组位置，防止野指针情况，试试
                            //curr = *cors2_array;

                            cors2_array= convertVector2array(cors2);

                            nmolty2++;
                            cout<<"new:cors2_array:"<<cors2_array[nmolty2 - 1][0]<<" "<<cors2_array[nmolty2 - 1][1]<<" "<<cors2_array[nmolty2 - 1][2]<<endl;

                            //free_2D_array(cors2_array,nmolty2);
                            //double **cors2_array = convertVector2array(cors2);


                            //从cors1_array中删除当前i粒子
                            vector<vector<double>>::iterator it;
                            it = cors1.begin();
                            int count = 0;
                            while(count != i)
                            {
                                it++;
                                count++;
                            }
                            it++;
                            cors1.erase(it);
                            cors1_array = convertVector2array(cors1);
                            //free_2D_array(cors1_array,nmolty1);
                            //double **cors1_array = convertVector2array(cors1);
                            nmolty1--;

                            total1.pot = total1.pot - partial_old.pot;
                            total2.pot = total2.pot + partial_new.pot;
                            x12_acc = x12_acc + 1;
                            cout<<"x12_acc："<<x12_acc<<endl;


                        }
                    }


                }

            }
            else
            {
                //try swapping 2->1
                x21_try = x21_try + 1;
                if(nmolty2 > 1)
                {
                    random_device seed_device;
                    default_random_engine  engine;
                    engine.seed(seed_device());
                    uniform_int_distribution<double> distr(0, nmolty2);


                    /*
                    int i;
                    srand((int)time(0));
                    i = rand()%nmolty1;
                     */
                    //随机选择box1中的任意一个粒子
                    int i = distr(engine);//i就是box1中要删去的第i行，坐标

                    //box1中需要删去的坐标必须是已经在box1中的坐标，但是新插入box2的坐标是任意的（因为插入任意位置）

                    //debug
                    cout<<"cors2_array:"<<cors2_array[i][0]<<" "<<cors2_array[i][1]<<" "<<cors2_array[i][2]<<endl;
                    cout<<"ri:"<<ri[0]<<" "<<ri[1]<<" "<<ri[2]<<endl;
                    cout<<"i:"<<i<<endl;


                    partial_old = PotentialType::potential_1(cors2_array[i],box_length,r_cut,cors2_array,nmolty2);
                    //if(!partial_old.ovr) throw "Overlap on particle removal";

                    partial_new = PotentialType::potential_1(ri,box_length,r_cut,cors1_array,nmolty1);
                    //cout<<"work"<<endl;
                    cout<<"box2：partial_old:"<<partial_old.pot<<endl;
                    cout<<"box2：partial_new:"<<partial_new.pot<<endl;

                    if(!partial_new.ovr)
                    {
                        delta = (partial_new.pot - partial_old.pot) / temperature;
                        delta = delta - log(pow(box_length,3)/(nmolty1 + 1));
                        delta = delta + log(pow(box_length,3)/nmolty2);
                        cout<<"delta:"<<delta<<endl;
                        if(metropolis(delta))
                        {
                            //ri的三个坐标用push_back推入二维向量vector cors1中
                            cors1[nmolty1].push_back(ri[0]);
                            cors1[nmolty1].push_back(ri[1]);
                            cors1[nmolty1].push_back(ri[2]);




                            //先记录原来的二维数组位置，防止野指针情况，试试
                            //curr = *cors2_array;

                            cors1_array= convertVector2array(cors1);

                            //free_2D_array(cors2_array,nmolty2);
                            //double **cors2_array = convertVector2array(cors2);
                            nmolty1++;

                            //从cors1_array中删除当前i粒子
                            vector<vector<double>>::iterator it;
                            it = cors2.begin();
                            int count = 0;
                            while(count != i)
                            {
                                it++;
                                count++;
                            }
                            it++;
                            cors2.erase(it);
                            cors2_array = convertVector2array(cors2);
                            //free_2D_array(cors1_array,nmolty1);
                            //double **cors1_array = convertVector2array(cors1);
                            nmolty2--;

                            total1.pot = total1.pot - partial_old.pot;
                            total2.pot = total2.pot + partial_new.pot;
                            x21_acc = x21_acc + 1;
                            cout<<"x21_acc："<<x21_acc<<endl;


                        }
                    }
                    /*
                    random_device seed_device;
                    default_random_engine  engine;
                    engine.seed(seed_device());
                    uniform_int_distribution<int> distr(0, nmolty2);

                    //随机选择box1中的任意一个粒子
                    int i = distr(engine);
                    cout<<"cors1_array:"<<cors1_array[i][0]<<" "<<cors1_array[i][1]<<" "<<cors1_array[i][2]<<endl;
                    cout<<"ri:"<<ri[0]<<" "<<ri[1]<<" "<<ri[2]<<endl;

                    cout<<"i:"<<i<<endl;
                    partial_old = PotentialType::potential_1(cors2_array[i],box_length,r_cut,cors2_array,nmolty2);
                    //if(!partial_old.ovr) throw "Overlap on particle removal";

                    partial_new = PotentialType::potential_1(ri,box_length,r_cut,cors2_array,nmolty2);
                    cout<<"box2：partial_old:"<<partial_old.pot<<endl;
                    cout<<"box2：partial_new:"<<partial_new.pot<<endl;

                    if(!partial_new.ovr)
                    {
                        delta = (partial_new.pot - partial_old.pot) / temperature;
                        delta = delta - log(pow(box_length,3)/(nmolty1 + 1));
                        delta = delta + log(pow(box_length,3) / nmolty2);

                        if(metropolis(delta))
                        {
                            //ri的三个坐标用push_back推入二维向量vector cors1中
                            cors1[nmolty2].push_back(ri[0]);
                            cors1[nmolty2].push_back(ri[1]);
                            cors1[nmolty2].push_back(ri[2]);

                            //先释放原来的cors1_array的内存，在把新的cors1 vector转换成增加了一行数据后的数组
                            free_2D_array(cors1_array,nmolty1);
                            double **cors1_array = convertVector2array(cors1);
                            nmolty1++;

                            //从cors2_array中删除当前i粒子
                            vector<vector<double>>::iterator it;
                            it = cors2.begin();
                            int count = 0;
                            while(count != i)
                            {
                                it++;
                                count++;
                            }
                            cors2.erase(it);

                            free_2D_array(cors2_array,nmolty2);
                            double **cors2_array = convertVector2array(cors2);
                            nmolty2--;

                            total1.pot = total1.pot - partial_old.pot;
                            total2.pot = total2.pot + partial_new.pot;
                            x21_acc = x21_acc + 1;
                            cout<<"x21_acc"<<x21_acc<<endl;


                        }
                    }
                */
                }
            }

        }

        if(x12_try > 0)
            x12_ratio = x12_acc/x12_try;
        else x12_ratio = 0.0;

        cout<<"x12_ratio:"<<x12_ratio<<endl;

        if(x21_try > 0)
            x21_ratio = x21_acc/x21_try;
        else x21_ratio = 0.0;
        cout<<"x21_ratio:"<<x21_ratio<<endl;

    }

        //volume move


        free_2D_array(atoms1_array,nmolty1);
         free_2D_array(cors1_array,nmolty1);
         vector<vector<double>>().swap(atoms1);
        vector<vector<double>>().swap(cors1);

    return 0;

}