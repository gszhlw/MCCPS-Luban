//
// Created by 张力文 on 2022/2/28.
//
#include<vector>
#include<iostream>
#include<random>

using namespace  std;
vector<double> random_translate_vector_vec(double dr_max, vector<double>& old)
{
    random_device seed_device;
    default_random_engine  engine;
    engine.seed(seed_device());
    uniform_real_distribution<double> distr(0, 1);

    vector<double> zeta;

    for(int i = 0; i < 3; i++)
    {

        zeta.push_back(distr(engine));
        cout<<"work"<<endl;
        zeta[i] = 2.0 * zeta[i] - 1.0;
        cout<<"zeta[i]:"<<zeta[i]<<endl;
        old[i] = old[i] + zeta[i]*dr_max;
        cout<<"old[i]:"<<old[i]<<endl;
    }
    return old;
}

int main()
{
   // vector<double> ri;
    double dr_max = 0.15;
    double box_length1 = 6.98864372;

    vector<vector<double>> cors1;
    cors1.resize(10);
    for (int k = 0; k < 10; ++k){
        cors1[k].resize(3);//每行为c列
    }
    double nums[10][3] = {
            {3.277722192,	-2.708689183,	-3.516056409},
            {-0.339364742,	-1.985222044,	2.888865293},
            {-1.972300937,	-0.596871757,	-2.967362437},
            {0.416892535,	0.743878206,	-1.395525409},
            {-2.316870615,	-2.978967151,	-2.705687392},
            {0.570337313,	2.185794358, 0.502898799},
            {1.273890981,	2.166174729,	-0.419976492},
            {-1.489210568,	2.270554197,	1.239211993},
            {-3.140619296,	3.172023134,	-1.933723419},
            {-1.763412159,	1.998193326,	-2.841737868}
    };
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 3; j++){
            cors1[i][j] = nums[i][j];
            cout<<cors1[i][j]<<" ";
        }
        cout<<endl;
    }
    for (int j = 0; j < 10; j++)
    {
        //double partial_old = 8928.51;
        //double partial_new = 9492.71;
        cout<<"cors[j]:"<<cors1[j][0]<<" "<<cors1[j][1]<<" "<<cors1[j][2]<<endl;
        vector<double> ri = random_translate_vector_vec(dr_max/box_length1,cors1[j]);
        //提示可能是数组越界或内存溢出

        cout  <<"ri:"<<ri[0] << " " << ri[1] << " " << ri[2] << endl;

    }


    return 0;
}