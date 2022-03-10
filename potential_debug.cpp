//
// Created by 张力文 on 2022/3/1.
//
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include<math.h>
#include<random>
#include<algorithm>
#include<iomanip>

using namespace std;

#define Kb 1.3806488e-23

class PotentialType{
public:
    double pot;
    bool ovr = false;
    static PotentialType potential_1(double *ri, double box_length, double r_cut, double **atoms, int n);
    PotentialType potential(double box_length, double r_cut, double **atoms,int n);

    static PotentialType potential_vec(double box_length, double r_cut, vector<vector<double>> &atoms,int n);
    static PotentialType potential_1_vec(vector<double> &ri, double box_length, double r_cut, vector<vector<double>> &atoms, int n);

};

void display_cors(vector<vector<double>>& cors)
{
    int rows = cors.size();
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cors[i].size(); j++)
        {
            cout<<cors[i][j]<<" ";
        }
        cout<<endl;
    }

}

PotentialType PotentialType:: potential_1_vec(vector<double> &ri, double box_length, double r_cut, vector<vector<double>> &atoms, int n)
{
    /*
     * 函数功能：输入一个原子的坐标，计算其原子间相互作用。
     * 输入：欲计算的某原子坐标、模拟盒的尺寸、截断范围、构型中原子坐标数组、原子个数
     * 输出：Potenetial类变量：partial
     */
    double sr2_ovr = 1.77;
    //double r_cut_box = r_cut / box_length;
    double r_cut_box = r_cut;
    double r_cut_box_sq = pow(r_cut_box,2);
    double box_sq = pow(box_length,2);





    double rij_sq;
    double sr2;
    double sr6;
    double sr12;
    bool ovr;

    PotentialType partial;
    PotentialType potential;

    partial.pot = 0.0;
    partial.ovr = false;

    vector<double> rij(3,0);

    for(int j = 0; j < n; j++)
    {

        if((ri[0] == atoms[j][0])&&(ri[1] == atoms[j][1])&&(ri[2] == atoms[j][2]))
        {
            //cout<<"same case->j:"<<j<<endl;
            continue;

        }
        else
        {
            //cout<<"different case->j:"<<j<<endl;
            //display_cors(atoms);
            //cout<<"atoms[j]:"<<setprecision(9)<<atoms[j][0]<<' '<<atoms[j][1]<<' '<<atoms[j][2]<<endl;
            rij[0] = ri[0] - atoms[j][0];
            rij[1] = ri[1] - atoms[j][1];
            rij[2] = ri[2] - atoms[j][2];

            //cout<<setprecision(12)<<"rij before PBC:"<<rij[0]<<' '<<rij[1]<<' '<<rij[2]<<endl;
            //在模拟盒尺寸为1时，应用周期性边界条件
           // rij[0] = rij[0] - round(rij[0]);
            //rij[1] = rij[1] - round(rij[1]);
            //rij[2] = rij[2] - round(rij[2]);

            //模拟盒尺寸不为1时，应用PBC
            rij[0] = rij[0] - box_length * round(rij[0]/box_length);
            rij[1] = rij[1] - box_length * round(rij[1]/box_length);
            rij[2] = rij[2] - box_length * round(rij[2]/box_length);
           // cout<<"rij after PBC:"<<setprecision(9)<<rij[0]<<' '<<rij[1]<<' '<<rij[2]<<endl;
            rij[0] = pow(rij[0],2);
            rij[1] = pow(rij[1],2);
            rij[2] = pow(rij[2],2);


           // rij_sq = sqrt(rij[0] + rij[1] + rij[2]);
            rij_sq = rij[0] + rij[1] + rij[2];
            //cout<<setprecision(9)<<"rij_sq: "<<setprecision(12)<<rij_sq<<endl;
            //cout<<setprecision(9)<<"r_cut_box_sq:"<<r_cut_box_sq<<endl;
            if(rij_sq < r_cut_box_sq)//check within cutoff
            {
                rij_sq = rij_sq * box_sq;
                //cout<<"rij_sq * box_sq: "<<setprecision(12)<<rij_sq<<endl;
                sr2 = 1.0 / rij_sq;//sigma = 1;
                //cout<<"sr2: "<<sr2<<endl;
                ovr = (sr2 > sr2_ovr);

                //check overlap
                if(ovr)
                {
                    partial.ovr = true;
                    return partial;
                }

                sr6 = pow(sr2, 3);
                //cout<<"sr6: "<<sr6<<endl;
                sr12 = pow(sr6, 2);
                //cout<<"sr12: "<<sr12<<endl;
                //potential = sr12 - sr6;
                potential.pot= sr12 - sr6;
                partial.pot = partial.pot + potential.pot;
                //cout<<"potential_1中的partial： ";
                //printf("%e\n",partial.pot);
            }
        }


    }

    partial.pot = 4.0 * partial.pot;//4 * epsilon，epsilon=1


    return partial;

}

PotentialType potential_vec(double box_length, double r_cut, vector<vector<double>> &atoms,int n)
{
    PotentialType total;
    total.pot = 0.0;
    total.ovr = false;

    PotentialType partial;
    partial.pot = 0.0;
    int count = 1;
    for(int i = 0; i < n; i++)
    {
        //cout<<"=====atoms["<<i<<"]:"<<atoms[i][0]<<" "<<atoms[i][1]<<" "<<atoms[i][2]<<"====="<<endl;
        partial = PotentialType::potential_1_vec(atoms[i],box_length,r_cut,atoms,n);

        // cout<<"work"<<endl;
        //cout<<"partial:"<<partial.pot<<endl;
        if(partial.ovr)
        {
            total.ovr = true;
            cout<<"*****************************************overlap!*****************************************"<<endl;
            break;
        }

            total.pot = total.pot + partial.pot;
            //cout<<" total："<<total.pot<<endl;

    }

    return total;
}

vector<double> random_translate_vector_vec(double dr_max, vector<double>& old)
{
    random_device seed_device;
    default_random_engine  engine;
    engine.seed(seed_device());
    uniform_real_distribution<double> distr(0, 1);

    vector<double> zeta(3);
    vector<double> new_vec(3);
    for(int i = 0; i < 3; i++)
    {
        //zeta[i] = distr(engine);
        zeta.push_back(distr(engine));
        zeta[i] = 2.0 * zeta[i] - 1.0;
        //cout<<"zeta[i]:"<<zeta[i]<<endl;
        new_vec[i] = old[i] + zeta[i]*dr_max;
        //cout<<"old[i]:"<<old[i]<<endl;
    }
    return new_vec;
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

int main() {

    double box_length1 = 7.0;
    double box_length2 = 7.0;
    //T*=0.85
    //double sig = 5.67e-8;
    double sig = 1.0;
    double epsilon = 1.62429e-23;
    //LRC, rc*=3.0
    double r_cut = 2.5*sig;
    double dr_max = 0.15;
    double temperature = 1.0;
    int nstep = 1;
    int nswap = 10;

    //double box_length1 = 7.0*sig;

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


    for(int i =0; i < nmolty1;i++)
    {
        for(int j = 0; j < atoms1[0].size(); j++)
        {
            atoms1[i][j]= *it++;

           // cout<<atoms1[i][j]<<' ';

        }
        //cout<<endl;
    }

    for(int i =0; i < nmolty1;i++)
    {
        for(int j = 0; j < cors1[0].size(); j++)
        {
            cors1[i][j] = atoms1[i][j+1];
            //cout<<setprecision(9)<<cors1[i][j]<<' ';

        }
       // cout<<endl;
    }


    //=========================================box2=======================================

    it = v2.begin();
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

            //cout<<atoms1[i][j]<<' ';

        }
        //cout<<endl;
    }

    for(int i =0; i < nmolty2;i++)
    {
        for(int j = 0; j < cors2[0].size(); j++)
        {

            cors2[i][j] = atoms2[i][j+1];
            //cout<<cors1[i][j]<<' ';

        }
         //cout<<endl;
    }


    PotentialType total1= potential_vec(box_length1,r_cut,cors1,nmolty1);
    if(total1.ovr) throw "Overlap in initial configuration 1";


    PotentialType total2 = potential_vec(box_length2,r_cut,cors2,nmolty2);
    if(total2.ovr) throw "Overlap in initial configuration 2";

    cout<<"total1:"<<total1.pot<<endl;
    cout<<"total2:"<<total2.pot<<endl;

    PotentialType partial_old;
    PotentialType partial_new;

    double delta = 0.0;
    vector<double> ri(3);

    int m_acc = 0;
    double m1_ratio = 0.0;
    double m2_ratio = 0.0;
    double x12_ratio = 0.0;
    double x21_ratio = 0.0;

    for(int i =0; i < nstep; i++)
    {
        //=========================================box1:translation move=======================================
        for (int j = 0; j < nmolty1; j++)
        {

            partial_old = PotentialType::potential_1_vec(cors1[j],box_length1,r_cut,cors1,nmolty1);

            //cout<<"before translation move:"<<endl;
            //display_cors(cors1);

            ri = random_translate_vector_vec(dr_max/box_length1,cors1[j]);
            //cout<<"work"<<endl;
            ri[0] = ri[0] - round(ri[0]);//周期性边界条件
            ri[1] = ri[1] - round(ri[1]);
            ri[2] = ri[2] - round(ri[2]);

            //cout<<"after translation move:"<<endl;
            //display_cors(cors1);

           // cout << setprecision(8) << "ri:"<<ri[0] << " " << ri[1] << " " << ri[2] << endl;
            partial_new = PotentialType::potential_1_vec(ri, box_length1, r_cut, cors1, nmolty1);
           // cout<<"========================box1：partial_old=====================:"<<partial_old.pot<<endl;
           // cout<<"=======================box1：partial_new======================:"<<partial_new.pot<<endl;

            if (!partial_new.ovr) {
                delta = partial_old.pot - partial_new.pot;
                delta = delta / temperature;
            }

            if (metropolis(delta)) {
                //更新能量
                total1.pot = total1.pot + partial_new.pot - partial_old.pot;

                //更新坐标
                cors1[i][0] = ri[0];
                cors1[i][1] = ri[1];
                cors1[i][2] = ri[2];

                m_acc = m_acc + 1;
               // cout << "box1：m_acc:" << m_acc << endl;
            }

        }
        m1_ratio = (double) m_acc / nmolty1;
        cout << "m1_ratio:" << m1_ratio << endl;

        //=========================================box1:translation move=======================================






        //=========================================box2:translation move=======================================
        for (int j = 0; j < nmolty2; j++) {

            partial_old = PotentialType::potential_1_vec(cors2[j],box_length2,r_cut,cors2,nmolty2);

            ri = random_translate_vector_vec(dr_max/box_length2,cors2[j]);
            ri[0] = ri[0] - round(ri[0]);//周期性边界条件
            ri[1] = ri[1] - round(ri[1]);
            ri[2] = ri[2] - round(ri[2]);
            //cout << setprecision(8) << "ri:"<<ri[0] << " " << ri[1] << " " << ri[2] << endl;
            //partial_new = PotentialType::potential_1(ri, box_length2, r_cut, cors2_array, nmolty2);
            partial_new = PotentialType::potential_1_vec(ri,box_length2,r_cut,cors2,nmolty2);
            //cout<<"========================box2：partial_old================:"<<partial_old.pot<<endl;
            //cout<<"========================box2：partial_new================:"<<partial_new.pot<<endl;

            if (!partial_new.ovr)
            {
                delta = partial_old.pot - partial_new.pot;
                delta = delta / temperature;
            }

            if (metropolis(delta))
            {
                //更新能量
                total2.pot = total2.pot + partial_new.pot - partial_old.pot;
                //更新坐标
                //ri[0] = cors2_array[i][0];
                //ri[1] = cors2_array[i][1];
                //ri[2] = cors2_array[i][2];
                cors2[i][0] = ri[0];
                cors2[i][1] = ri[1];
                cors2[i][2] = ri[2];


                m_acc = m_acc + 1;
                //cout << "box2：m_acc:" << m_acc << endl;
            }

        }

        m2_ratio = (double) m_acc / nmolty2;
        cout <<setprecision(4)<< "m2_ratio:" << m2_ratio << endl;
        //=========================================box2:translation move=======================================

        //=========================================swap move=======================================
        int x12_try = 0, x21_try = 0;
        int x12_acc = 0, x21_acc = 0;

        PotentialType  partial_old;
        PotentialType partial_new;

        for(int j = 0; j < nswap; j++)
        {
            //3个（0，1）均匀分布的随机数

            //double* ri = new double[3];
            vector <double> ri(3);
            //随机选择box1中的任意一个粒子

            random_device seed_device;
            default_random_engine  engine;
            uniform_real_distribution<double> distr(0, 1);


            for(int i = 0; i < 3; i++)
            {
                //ri为-0.5～0.5之间的数，3个为一组坐标,为欲交换到的目标坐标

                ri[i] = distr(engine);
                ri[i] = ri[i] - 0.5;

            }


            //random_device seed_device;
            //default_random_engine  engine;
            engine.seed(seed_device());
            uniform_real_distribution<double> distr_rand(0, 1);
            double rand_num = distr_rand(engine);
            //double rand_num = 0.4;
            cout<<"==============================="<<endl;
            cout<<"rand:"<<rand_num<<endl;

            //vector<int>rand_collect = randperm(10);
            //int num = 0;

            //swapping1->2

            if(rand_num < 0.5)
            {
                x12_try = x12_try + 1;
                if(nmolty1 > 1)
                {

                    random_device seed_device;
                    default_random_engine  engine;
                    engine.seed(seed_device());
                    uniform_int_distribution<int> distr(0, nmolty1 - 1);
                    /*
                    int i;
                    srand((int)time(0));
                    i = rand()%nmolty1;
                     */
                    //随机选择box1中的任意一个粒子

                    int i = distr(engine);//i就是box1中要删去的第i行，坐标

                    // int i = rand_collect[num++];

                    //vector中erase掉之后并不会紧凑！也就是说，如果下一次又有同样随机值i，vector中那个地方啥也没有，得保证随机值不重复
                    //box1中需要删去的坐标必须是已经在box1中的坐标，但是新插入box2的坐标是任意的（因为插入任意位置）

                    //debug
                    //cout<<setprecision(9)<<"cors1_array:"<<cors1_array[i][0]<<" "<<cors1_array[i][1]<<" "<<cors1_array[i][2]<<endl;
                    //cout<<"ri:"<<ri[0]<<" "<<ri[1]<<" "<<ri[2]<<endl;
                    //cout<<"i:"<<i<<endl;


                    //partial_old = PotentialType::potential_1(cors1_array[i],box_length,r_cut,cors1_array,nmolty1);
                    partial_old = PotentialType::potential_1_vec(cors1[i],box_length1,r_cut,cors1,nmolty1);
                    if(partial_old.ovr)
                    {
                        cout<< "Overlap on particle removal"<<endl;
                        continue;
                    }

                    // partial_new = PotentialType::potential_1(ri,box_length,r_cut,cors2_array,nmolty2);
                    partial_new = PotentialType::potential_1_vec(ri,box_length2,r_cut,cors2,nmolty2);
                    //cout<<"work"<<endl;
                    cout<<"box1：partial_old:"<<partial_old.pot<<endl;
                    cout<<"box1：partial_new:"<<partial_new.pot<<endl;

                    if(!partial_new.ovr)
                    {
                        delta = (partial_new.pot - partial_old.pot) / temperature;
                        delta = delta - log(pow(box_length2,3)/(nmolty2 + 1));
                        delta = delta + log(pow(box_length1,3)/nmolty1);
                        //cout<<"delta:"<<delta<<endl;

                        if(metropolis(delta))
                        {
                            //ri的三个坐标用push_back推入二维向量vector cors2中
                            cors2.resize(nmolty2 + 1);
                            //cout<<"nmolty2:"<<nmolty2<<endl;
                            cout<<"cors2 size:"<<cors2.size()<<endl;
                            //cout<<"cors2_capacity:"<<cors2.capacity()<<endl;

                            cors2[nmolty2].push_back(ri[0]);
                            cors2[nmolty2].push_back(ri[1]);
                            cors2[nmolty2].push_back(ri[2]);

                            cors1.erase(cors1.begin()+i);
                            cors1.resize(nmolty1 - 1);

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

                    //cout<<"cors2_array:"<<cors2_array[i][0]<<" "<<cors2_array[i][1]<<" "<<cors2_array[i][2]<<endl;
                    //cout<<"ri:"<<ri[0]<<" "<<ri[1]<<" "<<ri[2]<<endl;
                    //cout<<"i:"<<i<<endl;


                    //partial_old = PotentialType::potential_1(cors2_array[i],box_length,r_cut,cors2_array,nmolty2);
                    partial_old = PotentialType::potential_1_vec(cors2[i], box_length2, r_cut, cors2, nmolty2);
                    if(partial_old.ovr)
                    {
                        cout<<"Overlap on particle removal";
                        continue;
                    }

                    //partial_new = PotentialType::potential_1(ri,box_length,r_cut,cors1_array,nmolty1);
                    partial_new = PotentialType::potential_1_vec(ri,box_length1,r_cut,cors1,nmolty1);
                    //cout<<"work"<<endl;
                    cout<<"box2：partial_old:"<<partial_old.pot<<endl;
                    cout<<"box2：partial_new:"<<partial_new.pot<<endl;

                    if(!partial_new.ovr)
                    {
                        delta = (partial_new.pot - partial_old.pot) / temperature;
                        delta = delta - log(pow(box_length1,3)/(nmolty1 + 1));
                        delta = delta + log(pow(box_length2,3)/nmolty2);
                        //cout<<"delta:"<<delta<<endl;
                        if(metropolis(delta))
                        {
                            //ri的三个坐标用push_back推入二维向量vector cors1中
                            //nmolty1++;
                            cors1.resize(nmolty1 + 1);
                            //cout<<"nmolty1:"<<nmolty1<<endl;
                            //cout<<"cors1 size:"<<cors1.size()<<endl;
                            cors1[nmolty1].push_back(ri[0]);
                            cors1[nmolty1].push_back(ri[1]);
                            cors1[nmolty1].push_back(ri[2]);

                            //cout<<"cors1:"<<cors1[nmolty1][0]<<" "<<cors1[nmolty1][1]<<" "<<cors1[nmolty1][2]<<endl;

                            nmolty1++;


                            cors2.erase(cors2.begin()+i);
                            cors2.resize(nmolty2 - 1);
                            //cors2_array = convertVector2array(cors2);
                            //free_2D_array(cors1_array,nmolty1);;;;
                            //double **cors1_array = convertVector2array(cors1);
                            nmolty2--;

                            // cout<<"cors1 size:"<<cors1.size()<<endl;
                            //cout<<"cors2 size:"<<cors2.size()<<endl;
                            total1.pot = total1.pot - partial_old.pot;
                            total2.pot = total2.pot + partial_new.pot;
                            x21_acc = x21_acc + 1;
                            cout<<"x21_acc："<<x21_acc<<endl;


                        }
                    }

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

        //=========================================swap move=======================================


        //=========================================volume move=======================================
/*
        double v_ratio = 0.0;

        //dv为（-1，1）间的均匀分布随机数
        random_device seed_device;
        default_random_engine  engine;
        uniform_real_distribution<double> distr(0, 1);
        double dv = 2.0 * distr(engine) - 1;

        double vol1_old = pow(box_length1,3);
        double vol2_old = pow(box_length2, 3);

        double vol1_new = vol1_old - dv;
        double vol2_new = vol2_old + dv;

        double box1_new_length = pow(vol1_new, (1/3));
        double box2_new_length = pow(vol2_new, (1/3));

        if(min(box1_new_length,box2_new_length) < 2.0 * r_cut)
            throw "Box length too small";

        PotentialType total1_new;
        PotentialType total2_new;
        total1_new = PotentialType::potential_vec(box1_new_length,r_cut,cors1,nmolty1);
        total2_new = PotentialType::potential_vec(box2_new_length, r_cut,cors2,nmolty2);

        cout<<"volume move:total1_new:"<<total1_new.pot<<endl;
        cout<<"volume move:total2_new:"<<total2_new.pot<<endl;
        double delta = 0.0;
        if(!(total1_new.ovr or total2_new.ovr))
        {
            delta = total1_new.pot + total2_new.pot - total1.pot - total2.pot;
            delta = delta / temperature;
            delta = delta - nmolty1 * log(vol1_new/vol1_old);
            delta = delta - nmolty2 * log(vol2_new/vol2_old);

            if(metropolis(delta))
            {
                total1 = total1_new;
                total2 = total2_new;

                box_length1 = box1_new_length;
                box_length2 = box2_new_length;
                v_ratio = 1.0;
            }
        }

        //=========================================volume move=======================================

*/
    }

    return 0;
}

