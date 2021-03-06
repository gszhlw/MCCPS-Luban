//
// Created by 张力文 on 2022/2/19.
//

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include<math.h>
#include<random>
#include<algorithm>
#include<iomanip>



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
    PotentialType potential(double box_length, double r_cut, double **atoms,int n);
    //friend PotentialType operator+(PotentialType& pot, PotentialType& vir);
    static PotentialType potential_vec(double box_length, double r_cut, vector<vector<double>> &atoms,int n);
    static PotentialType potential_1_vec(vector<double> &ri, double box_length, double r_cut, vector<vector<double>> &atoms, int n);



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
PotentialType PotentialType:: potential_1_vec(vector<double> &ri, double box_length, double r_cut, vector<vector<double>> &atoms, int n)
{
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

    vector<double> rij(3);
    for(int j = 0; j < n; j++)
    {
        if((ri[0] == atoms[j][0])&&(ri[1] == atoms[j][1])&&(ri[2] == atoms[j][2]))
        {
            continue;
        }

        //rij[0] = ri[0] - atoms[j][0];
        //rij[1] = ri[1] - atoms[j][1];
        //rij[2] = ri[2] - atoms[j][2];
        rij.push_back(ri[0] - atoms[j][0]);
        rij.push_back(ri[1] - atoms[j][1]);
        rij.push_back(ri[2] - atoms[j][2]);
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

PotentialType potential_vec(double box_length, double r_cut, vector<vector<double>> &atoms,int n)
{
    PotentialType total;
    total.pot = 0.0;
    PotentialType partial;
    partial.pot = 0.0;

    for(int i = 0; i < n; i++)
    {

        partial = PotentialType::potential_1_vec(atoms[i],box_length,r_cut,atoms,n);
       // cout<<"work"<<endl;
        if(!partial.ovr)
        {
            //cout<<"partial:"<<partial.pot<<endl;
            total.pot = total.pot + partial.pot;
            //cout<<count++<<" "<<total.pot<<endl;
        }
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

vector<double> random_translate_vector_vec(double dr_max, vector<double>& old)
{
    random_device seed_device;
    default_random_engine  engine;
    engine.seed(seed_device());
    uniform_real_distribution<double> distr(0, 1);

    vector<double> zeta;
    for(int i = 0; i < 3; i++)
    {
        //zeta[i] = distr(engine);
        zeta.push_back(distr(engine));
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

/*
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
*/
//释放数组内存函数
/*
void free_2D_array(double **array, int rows)
{
    for(int i = 0; i < rows; i++)
    {
        delete[] array[i];
    }

    delete[] array;
}
*/
/*
vector<int> randperm(int Num)
{
    vector<int> temp;
    for (int i = 0; i < Num; ++i)
    {
        temp.push_back(i + 1);
    }

    random_shuffle(temp.begin(), temp.end());

    for (int i = 0; i < temp.size(); i++)
    {
        cout << temp[i] << " ";
    }
     */
    //return temp;
//}


int main() {


    //模拟参数
    double box_length1 = 6.98864372;
    double box_length2 = 7.09910868;
    double r_cut = 2.0;
    double dr_max = 0.15;
    double temperature = 1.0;
    int nstep = 1;
    int nswap = 10;


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

    atoms1.reserve(500);
    cors1.reserve(500);


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
/*
    double **atoms1_array = convertVector2array(atoms1);
    double **cors1_array = convertVector2array(cors1);
*/

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

    atoms2.reserve(500);
    cors2.reserve(500);

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
       // cout<<endl;
    }
/*
    double **atoms2_array = convertVector2array(atoms1);
    double **cors2_array = convertVector2array(cors1);

*/

    //PotentialType total1= potential(box_length1,r_cut,cors1_array,nmolty1);
    PotentialType total1= potential_vec(box_length1,r_cut,cors1,nmolty1);
    //重叠判断
    if(total1.ovr) throw "Overlap in initial configuration 1";

    //PotentialType total2 = potential(box_length2,r_cut,cors2_array,nmolty2);
    PotentialType total2 = potential_vec(box_length2,r_cut,cors2,nmolty2);
    if(total2.ovr) throw "Overlap in initial configuration 2";

    PotentialType partial_old;
    PotentialType partial_new;

    double delta = 0.0;
    //double *ri = new double[3];
    vector<double> ri;
    cout<<"total1:"<<total1.pot<<endl;
    cout<<"total2:"<<total2.pot<<endl;

    int m_acc = 0;
    double m1_ratio = 0.0;
    double m2_ratio = 0.0;
    double x12_ratio = 0.0;
    double x21_ratio = 0.0;

    for(int i =0; i < nstep; i++)
    {
        //box1
        for (int j = 0; j < nmolty1; j++) {
            //partial_old = PotentialType::potential_1(cors1_array[j], box_length1, r_cut, cors1_array, nmolty1);

            partial_old = PotentialType::potential_1_vec(cors1[j],box_length1,r_cut,cors1,nmolty1);
            cout<<"box1：partial_old:"<<partial_old.pot<<endl;
            //ri = random_translate_vector(dr_max / box_length1, cors1_array[j]);
            ri = random_translate_vector_vec(dr_max/box_length1,cors1[j]);
            //cout<<"work"<<endl;
            ri[0] = ri[0] - round(ri[0]);//周期性边界条件
            ri[1] = ri[1] - round(ri[1]);
            ri[2] = ri[2] - round(ri[2]);

            //partial_new = PotentialType::potential_1(ri, box_length1, r_cut, cors1_array, nmolty1);
            partial_new = PotentialType::potential_1_vec(ri, box_length1, r_cut, cors1, nmolty1);

            cout<<"box1：partial_new:"<<partial_new.pot<<endl;
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
                cors1[i][0] = ri[0];
                cors1[i][1] = ri[1];
                cors1[i][2] = ri[2];

                m_acc = m_acc + 1;
                //cout << "box1：m_acc:" << m_acc << endl;
            }


        }

        m1_ratio = (double) m_acc / nmolty1;
        cout << "m1_ratio:" << m1_ratio << endl;

        m_acc = 0;

        //box2
        for (int j = 0; j < nmolty2; j++) {
            //partial_old = PotentialType::potential_1(cors2_array[j], box_length2, r_cut, cors2_array, nmolty2);
            partial_old = PotentialType::potential_1_vec(cors2[j],box_length2,r_cut,cors2,nmolty2);
            //ri = random_translate_vector(dr_max / box_length2, cors2_array[j]);
            ri = random_translate_vector_vec(dr_max/box_length2,cors2[j]);
            ri[0] = ri[0] - round(ri[0]);//周期性边界条件
            ri[1] = ri[1] - round(ri[1]);
            ri[2] = ri[2] - round(ri[2]);

            //partial_new = PotentialType::potential_1(ri, box_length2, r_cut, cors2_array, nmolty2);
            partial_new = PotentialType::potential_1_vec(ri,box_length2,r_cut,cors2,nmolty2);
            //cout<<"box2：partial_old:"<<partial_old.pot<<endl;
            //cout<<"box2：partial_new:"<<partial_new.pot<<endl;

            if (!partial_new.ovr)
            {
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
                cors2[i][0] = ri[0];
                cors2[i][1] = ri[1];
                cors2[i][2] = ri[2];


                m_acc = m_acc + 1;
                //cout << "box2：m_acc:" << m_acc << endl;
            }

        }

        m2_ratio = (double) m_acc / nmolty2;
        cout <<setprecision(4)<< "m2_ratio:" << m2_ratio << endl;

        //swap move

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

                            //cout<<"cors2:"<<cors2[nmolty2][0]<<" "<<cors2[nmolty2][1]<<" "<<cors2[nmolty2][2]<<endl;

                            //先释放原来的cors2_array的内存，在把新的cors2 vector转换成增加了一行数据后的数组
                            //先记录原来的二维数组位置，防止野指针情况，试试
                            //curr = *cors2_array;

                            //这里有问题，原因八成应该就是因为直接覆盖的方式来创建新数组行不通
                            //cors2_array= convertVector2array(cors2);

                            //nmolty2++;
                            //cout<<"new:cors2_array:"<<cors2_array[nmolty2 - 1][0]<<" "<<cors2_array[nmolty2 - 1][1]<<" "<<cors2_array[nmolty2 - 1][2]<<endl;

                            //free_2D_array(cors2_array,nmolty2);
                            //double **cors2_array = convertVector2array(cors2);
                            nmolty2++;

                            //从cors1_array中删除当前i粒子
                            /*
                            vector<vector<double>>::iterator it;
                            it = cors1.begin();
                            int count = 0;
                            while(count != i)
                            {
                                it++;
                                count++;
                            }

                            cors1.erase(it);
                             */
                           cors1.erase(cors1.begin()+i);
                            cors1.resize(nmolty1 - 1);

                            //free_2D_array(cors1_array,nmolty1);

                            //double **cors1_array = convertVector2array(cors1);
                            //free_2D_array(cors1_array,nmolty1);
                            //double **cors1_array = convertVector2array(cors1);
                            nmolty1--;

                            //cout<<"cors1 size:"<<cors1.size()<<endl;
                            //cout<<"cors2 size:"<<cors2.size()<<endl;


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


                            //先记录原来的二维数组位置，防止野指针情况，试试
                            //curr = *cors2_array;

                            //cors1_array= convertVector2array(cors1);

                            //free_2D_array(cors2_array,nmolty2);
                            //double **cors2_array = convertVector2array(cors2);
                            nmolty1++;

                            //从cors1_array中删除当前i粒子
                            /*
                            vector<vector<double>>::iterator it;
                            it = cors2.begin();
                            int count = 0;
                            while(count != i)
                            {
                                it++;
                                count++;
                            }
                            it++;
                             */
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

        //free_2D_array(cors1_array,nmolty1);
        //double** cors1_array = convertVector2array(cors1);
        //free_2D_array(cors2_array,nmolty2);
        //double** cors2_array = convertVector2array(cors2);

        //cout<<"cors1_array:"<<cors1_array[0][0]<<" "<<cors1_array[0][1]<<" "<<cors1_array[0][2]<<endl;
        //cout<<"cors2_array:"<<cors2_array[0][0]<<" "<<cors2_array[0][1]<<" "<<cors2_array[0][2]<<endl;

        //display
        /*
        for(int i = 0; i < nmolty1; i++)
        {
            for(int j = 0; j < 3;j++)
            {
                cout<<cors1_array[i][j]<<" ";
            }
            cout<<endl;
        }

        for(int i = 0; i < nmolty2; i++)
        {
            for(int j = 0; j < 3;j++)
            {
                cout<<cors2_array[i][j]<<" ";
            }
            cout<<endl;
        }
    */

        if(x12_try > 0)
            x12_ratio = x12_acc/x12_try;
        else x12_ratio = 0.0;

        cout<<"x12_ratio:"<<x12_ratio<<endl;

        if(x21_try > 0)
            x21_ratio = x21_acc/x21_try;
        else x21_ratio = 0.0;
        cout<<"x21_ratio:"<<x21_ratio<<endl;

        //volume move
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




    }



/*
        free_2D_array(atoms1_array,nmolty1);
         free_2D_array(cors1_array,nmolty1);
         */
         vector<vector<double>>().swap(atoms1);
        vector<vector<double>>().swap(cors1);

    return 0;

}