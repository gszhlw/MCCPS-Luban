//
// Created by 张力文 on 2021/11/26.
//Metropolis采样算法解决的问题是：从一个复杂的目标分布获取近似的样本
//

#include<iostream>
#include<random>
#define k 1.3806526E-23
#difine T 500.0

using namespace std;
double get_randnum();
double get_energy();
int main()
{
    vector<double> rand_energy(10);
    vector<double> rand_num(10);
    double energy_diff = 0.0;
    double energy_final = 0.0;
    double boltzmann_weight = 0.0;
    double energy_temp = 0.0;

    for(int n = 0; n < 10; ++n)
    {
        rand_energy[n] = get_energy();
        rand_num[n] = get_randnum();
        //cout<<get_energy()<<' ';

    }

    energy_diff = rand_energy[1] - rand_energy[0];
    boltzmann_weight = exp(-energy_diff/k * T);

    if(boltzmann_weight > rand_num[0])
    {
        energy_temp = rand_energy[1];
    }
    else
    {

    }
/*
    for(vector<double>::iterator it = rand_energy.begin(); it != rand_energy.end(); ++it )
    {
        cout<<*it<<' ';
    }
    cout << '\n';

    for(int i = 0; i < rand_num.size(); ++i)
    {
        cout <<rand_num[i]<<" ";
    }
*/



}

double get_randnum()
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<>dis(0,1);
    return dis(gen);
}

double get_energy()
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<>dis(-100.0,100.0);
    return dis(gen);
}