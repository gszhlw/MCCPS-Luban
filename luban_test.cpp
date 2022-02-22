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

vector<double> read_cnf_atoms(string filename)
{
    ifstream myfile("/Users/zlw/ClionProjects/MCCPS Luban/"+filename);
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

class PotentialType{
public:
    double pot;
    bool ovr = false;
    //PotentialType(double pot, bool ovr);
    static PotentialType potential_1(double *ri, double box_length, double r_cut, double **atoms, int n);
    PotentialType operator+(PotentialType &p) const;
    //PotentialType operator*(PotentialType &p1,PotentialType &p2) const;
    PotentialType potential(double box_length, double r_cut, double r[][3]);
    //friend PotentialType operator+(PotentialType& pot, PotentialType& vir);

};


//PotentialType::PotentialType(double pot, bool ovr) : pot(pot), ovr(ovr){}

PotentialType PotentialType::operator+(PotentialType &p) const
{

    PotentialType temp;
    temp.pot=temp.pot+p.pot;
    temp.ovr=temp.ovr && (p.ovr);
    return temp;
}



PotentialType PotentialType:: potential_1(double *ri, double box_length, double r_cut, double **atoms, int n)
{
    double sr2_ovr = 1.77;
    double r_cut_box = r_cut / box_length;
    double r_cut_box_sq = pow(r_cut_box,2);
    double box_sq = pow(box_length,2);

    //double potential = 0.0;

    double rij_sq;
    double sr2;
    double sr6;
    double sr12;
    bool ovr;
    //int n = (sizeof(atoms)/sizeof(atoms[0][0]));
    //int n = 10;
    //cout<<"potential_1中的n："<<n<<endl;
    //double partial = 0.0;
    PotentialType partial;
    PotentialType potential;

    //partial.pot = 0.0;
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
        rij_sq = pow(rij[0] + rij[1] + rij[2],2);

        if(rij_sq < r_cut_box_sq)
        {
            rij_sq = rij_sq * box_sq;
            //cout<<"rij_sp: "<<rij_sq<<endl;
            sr2 = 1.0 / rij_sq;
            //cout<<"sr2: "<<sr2<<endl;
            ovr = sr2 > sr2_ovr;
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
            //cout<<"potential_1中的partial： "<<partial<<endl;
        }


    }
    partial.pot = 4.0 * partial.pot;

    return partial;


}


PotentialType potential(double box_length, double r_cut, double **atoms,int n)
{
    //double total = 0.0;
    PotentialType total;
    total.pot = 0.0;
    //double partial;
    PotentialType partial;
    partial.pot = 0.0;
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
        partial = PotentialType::potential_1(atoms[i],box_length,r_cut,atoms,n);
        //cout<<"partial:"<<partial<<endl;
        total.pot = total.pot + partial.pot;
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

    vector<double> v1 = read_cnf_atoms("input_box1.txt");
    vector<double> v2 = read_cnf_atoms("input_box2.txt");

    vector<double>::iterator it;

    int nmolty1;
    int nmolty2;

    int nstep = 1;
    int nswap = 10;
    double box_length = 6.98864372;
    double r_cut = 2.0;
    double dr_max = 0.15;

    double temperature = 1.0;
    //box1
    it = v1.begin();
    //cout<<"nmolty: "<<*it<<endl;

    nmolty1 = *it;
    cout<<"box1模拟原子个数： "<<nmolty1<<endl;
    it++;

    double** atoms1 = new double*[nmolty1];
    double** cors1 = new double* [nmolty1];
   //vector< vector<double>> atoms1(nmolty1,4);
    //vector<vector<double> > atoms1(nmolty1, vector<int>(4)); //定义二维动态数组5行6列

    //vector<vector<double>> cors1(nmolty1,3);
    for(int i = 0; i < nmolty1; i++)
    {
        atoms1[i] = new double[4];
        cors1[i] = new double[3];
    }

    for(it; it != v1.end();)
    {
        for(int i = 0; i < nmolty1; i++)
        {
            for (int j = 0; j < 4; j++)
            {

                atoms1[i][j] = *it;
                it++;
                //cout<<' '<<atoms[i][j];
            }

            //cout<<endl;
        }

    }

    for(int i = 0; i < nmolty1; i++)
    {
        for(int j = 0; j < 3; j++)
        {

            cors1[i][j] = atoms1[i][j+1];
            //cout<<cors1[i][j]<<' ';
        }
        //cout<<endl;
    }



    //box2
    it = v2.begin();

    nmolty2 = *it;
    cout<<"box2模拟原子个数： "<<nmolty2<<endl;
    it++;


    double** atoms2 = new double*[nmolty1];
    double** cors2 = new double* [nmolty1];

    for(int i = 0; i < nmolty2; i++)
    {
        atoms2[i] = new double[4];
        cors2[i] = new double[3];
    }

    for(it; it != v2.end();)
    {
        for(int i = 0; i < nmolty2; i++)
        {
            for (int j = 0; j < 4; j++)
            {

                atoms2[i][j] = *it;
                it++;
                //cout<<' '<<atoms[i][j];
            }

            //cout<<endl;
        }

    }

    for(int i = 0; i < nmolty2; i++)
    {
        for(int j = 0; j < 3; j++)
        {

            cors2[i][j] = atoms2[i][j+1];
            //cout<<cors2[i][j]<<' ';
        }
        //cout<<endl;
    }




    //double total1 = potential(box_length,r_cut,cors1,nmolty1);
    //double total2 = potential(box_length,r_cut,cors2,nmolty2);
    PotentialType total1= potential(box_length,r_cut,cors1,nmolty1);
    if(total1.ovr) throw "Overlap in initial configuration 1";
    PotentialType total2 = potential(box_length,r_cut,cors2,nmolty2);
    if(total2.ovr) throw "Overlap in initial configuration 2";
    //double partial_old = 0.0;
    //double partial_new = 0.0;
    PotentialType partial_old;
    PotentialType partial_new;

    double delta = 0.0;
    //double *ri = new double[3];

    cout<<"total1:"<<total1.pot<<endl;
    cout<<"total2:"<<total2.pot<<endl;

    int m_acc = 0;
    double m1_ratio = 0.0;
    double m2_ratio = 0.0;

    for(int i =0; i < nstep; i++)
    {

        //two boxes, tranlation move
        //box1
        for(int j = 0; j < nmolty1; j++)
        {
            partial_old = PotentialType::potential_1(cors1[j],box_length,r_cut,cors1,nmolty1);
            ri = random_translate_vector(dr_max/box_length,cors1[j]);
            ri[0] = ri[0] - round(ri[0]);//周期性边界条件
            ri[1] = ri[1] - round(ri[1]);
            ri[2] = ri[2] - round(ri[2]);
            //cout<<ri[0]<<" "<<ri[1]<<" "<<ri[2]<<endl;


            partial_new =PotentialType:: potential_1(ri,box_length,r_cut,cors1,nmolty1);
            cout<<"box1：partial_old:"<<partial_old.pot<<endl;
            cout<<"box1：partial_new:"<<partial_new.pot<<endl;

            if(!partial_new.ovr)
            {
                delta = partial_old.pot - partial_new.pot;
                delta = delta / temperature;
            }

            //cout<<"delta:"<<delta<<endl;
            if(metropolis(delta))
            {
                //更新能量
                total1.pot = total1.pot + partial_new.pot - partial_old.pot;
                //更新坐标
                ri[0] = cors1[i][0];
                ri[1] = cors1[i][1];
                ri[2] = cors1[i][2];

                m_acc = m_acc + 1;
                cout<<"box1：m_acc:"<<m_acc<<endl;
            }


        }

        m1_ratio = (double)m_acc / nmolty1;
        cout<<"m1_ratio:"<<m1_ratio<<endl;

        //box2
        for(int j = 0; j < nmolty2; j++)
        {
            partial_old = PotentialType::potential_1(cors2[j],box_length,r_cut,cors2,nmolty2);
            ri = random_translate_vector(dr_max/box_length,cors2[j]);
            ri[0] = ri[0] - round(ri[0]);//周期性边界条件
            ri[1] = ri[1] - round(ri[1]);
            ri[2] = ri[2] - round(ri[2]);
            //cout<<ri[0]<<" "<<ri[1]<<" "<<ri[2]<<endl;


            partial_new = PotentialType::potential_1(ri,box_length,r_cut,cors2,nmolty2);

            cout<<"box2：partial_old:"<<partial_old.pot<<endl;
            cout<<"box2：partial_new:"<<partial_new.pot<<endl;

            if(!partial_new.ovr)
            {
                delta = partial_old.pot - partial_new.pot;
                delta = delta / temperature;
            }


            //cout<<"delta:"<<delta<<endl;
            if(metropolis(delta))
            {
                //更新能量
                total2.pot = total2.pot +partial_new.pot - partial_old.pot;
                //更新坐标
                ri[0] = cors2[i][0];
                ri[1] = cors2[i][1];
                ri[2] = cors2[i][2];

                m_acc = m_acc + 1;
                cout<<"box2：m_acc:"<<m_acc<<endl;
            }


        }
        m2_ratio =(double) m_acc / nmolty2;
        cout<<"m2_ratio:"<<m2_ratio<<endl;

        int x12_try = 0, x21_try = 0;
        int x12_acc = 0, x21_acc = 0;

        //swap move

              for(int j = 0; j < nstep; j++)
              {
                  //3个（0，1）均匀分布的随机数
                  random_device seed_device;
                  default_random_engine  engine;
                  engine.seed(seed_device());
                  uniform_real_distribution<double> distr(0, 1);
                  //uniform_real_distribution<double> distr_2(0, 1);

                  PotentialType  partial_old;
                  PotentialType partial_new;

                  double* ri = new double[3];
                  for(int i = 0; i < 3; i++) {
                      //ri为-0.5～0.5之间的数，3个为一组坐标,为欲交换到的目标坐标
                      ri[i] = distr(engine);
                      ri[i] = ri[i] - 0.5;
                  }
                  double rand = distr(engine);


                  //swapping1->2
                  if(rand < 0.5)
                  {
                      x12_try = x12_try + 1;
                      if(nmolty1 > 1)
                      {
                          random_device seed_device;
                          default_random_engine  engine;
                          engine.seed(seed_device());
                          uniform_int_distribution<int> distr(0, nmolty1);
                          int i = distr(engine);

                          partial_old = PotentialType::potential_1(cors1[i],box_length,r_cut,cors1,nmolty1);
                          if(!partial_old.ovr) throw "Overlap on particle removal";

                          partial_new = PotentialType::potential_1(ri,box_length,r_cut,cors1,nmolty1);

                          if(!partial_new.ovr)
                          {
                              delta = (partial_new.pot - partial_old.pot) / temperature;
                              delta = delta - log(pow(box_length,3)/(nmolty2 + 1));
                              delta = delta + log(pow(box_length,3) / nmolty1);

                              if(metropolis(delta))
                              {
                                  //转换成vector扩容后再转换回数组
                                  vector<double>v;
                                  for(int i = 0; i < nmolty1; i++)
                                  {
                                      for(int j = 0; j < 3; j++)
                                      {
                                          v.push_back(cors1[i][j]);
                                      }

                                  }
                                  /*
                                  r2      = np.append ( r2, ri[np.newaxis,:], 0 ) # Add new particle to r2 array
                                  n2      = r2.shape[0]                           # New value of N2
                                  r1      = np.copy(rj)                           # Delete particle from r1 array
                                  n1      = r1.shape[0]                           # New value of N1
                                  total1  = total1 - partial_old                  # Update total values
                                  total2  = total2 + partial_new                  # Update total values
                                  x12_acc = x12_acc + 1                           # Increment 1->2 move counter
                                   */

                                  //cors2[nmolty2][0] =
                              }
                          }

                      }
                  }

                  //swapping 2->1
                  else{

                  }

              }


    }


    for(int i = 0; i < nmolty1; i++)
    {
        delete[] cors1[i];
    }
    delete[] cors1;



    for(int i = 0; i < nmolty2; i++)
    {
        delete[] cors2[i];
    }
    delete[] cors2;


    for(int i = 0; i < nmolty1; i++)
    {
        delete[] atoms1[i];
    }
    delete[] atoms1;


    for(int i = 0; i < nmolty2; i++)
    {
        delete[] atoms2[i];
    }
    delete[] atoms2;

    return 0;
}