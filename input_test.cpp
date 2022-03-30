//
// Created by 张力文 on 2022/5/2.
//
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cstdio>
#include <math.h>
#include <random>
#include <algorithm>
#include <iomanip>

using namespace std;

#define k_boltz 1.987206504191549E-003

// 设置盒子内的原子个数
int n_atoms = 500;

// 设置执行MC移动的次数
const int num_moves = 2.5E8;

// 设置盒子尺寸 (单位：A)
double box_length = 271;
double box_size[3] = { box_length, box_length, box_length};

//原子的最大位移值
double max_translate = 0.5;  // angstroms

//Ar原子的L-J参数
const double sigma = 3.41; // angstroms
const double epsilon = 0.697;// kcal mol-1

//温度设置为室温
const double temperature = 102;

//截断距离设置为3*sigma
const double r_cut = 3.0*sigma;

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

// 执行周期性边界条件
double make_periodic(double x, double box)
{


    while (x < -0.5*box)
    {
        x = x + box;
    }

    while (x > 0.5*box)
    {
        x = x - box;
    }

    //x = x - box* round(x/box);

    return x;
}

// 保证构型坐标都位于盒子内部
double wrap_into_box(double x, double box)
{
    while (x > box)
    {
        x = x - box;
    }

    while (x < 0)
    {
        x = x + box;
    }

    return x;
}

// 将坐标打印至pdb文件中
void print_pdb(vector<vector<double>>coords, int n_atoms, int move)
{
    char filename[128];

    snprintf(filename, 128, "output%000006d.pdb", move);

    FILE *f = fopen(filename, "w");

    fprintf(f, "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00\n",
            box_size[0], box_size[1], box_size[2]);

    for (int i = 0; i < n_atoms; i = i + 1)
    {
        fprintf(f, "ATOM  %5d  Kr   Kr     1    %8.3f%8.3f%8.3f  1.00  0.00          Kr\n",
                i+1, coords[i][0], coords[i][1], coords[i][2]);
        fprintf(f, "TER\n");
    }

    fclose(f);
}

// 计算原子能量
double calculate_energy(vector<vector<double>>coords, int n_atoms, double *box_size,
                        double sigma, double epsilon, double r_cut)
{
    // 循环所有原子对并计算LJ势能
    double total_energy = 0;
    double vir = 0.0;
    //double r_cut_box = r_cut/box_size[0];
    //double r_cut_box_sq = pow(r_cut_box,2);
    double r_cut_sq = pow(r_cut,2);
    double box_sq = pow(box_size[0],2);

    for (int i = 0; i < n_atoms-1; i = i + 1)
    {
        for (int j = i+1; j < n_atoms; j = j + 1)
        {
            double delta_x = coords[j][0] - coords[i][0];
            double delta_y = coords[j][1] - coords[i][1];
            double delta_z = coords[j][2] - coords[i][2];

            // PBC
            delta_x = make_periodic(delta_x, box_size[0]);
            delta_y = make_periodic(delta_y, box_size[1]);
            delta_z = make_periodic(delta_z, box_size[2]);

            //r2 = delx^2+dely^2+delz^2
            double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                        (delta_z*delta_z);

            // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
            // r_cut = 3*sigma
            // r_cut_sq = r_cut^2
            if(r2 < r_cut_sq)
            {
                //r2 = r2 * box_sq;
                const double sig2_over_r2 = (sigma*sigma) / r2;//reduced units
                const double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
                const double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;
                // const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );
                const double e_lj = 4.0 *  ( sig12_over_r12 - sig6_over_r6 );//U* = U/epsilon

                total_energy = total_energy + e_lj;

            }
        }
    }

    // 返回原子的总能量
    return total_energy;
}

void copy_coordinates(vector<vector<double>>from, vector<vector<double>>to)
{
    for (int i=0; i<n_atoms; ++i)
    {
        to[i][0] = from[i][0];
        to[i][1] = from[i][1];
        to[i][2] = from[i][2];
    }
}
double rand(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}
int main(int argc, const char **argv)
{
    vector<double> v = read_cnf_atoms("input_box.txt");



    vector<double>::iterator it = v.begin();



    //it = v1.begin();
    //cout<<"nmolty: "<<*it<<endl;

    //nmolty1 = *it;
    //cout<<"box1模拟原子个数： "<<nmolty1<<endl;

    //it++;


    vector<vector<double>> cors1(n_atoms, vector<double>(3));
    vector<vector<double>> old_cors1(n_atoms, vector<double>(3));
    while(it != v.end())
    {
        for(int i =0; i < n_atoms;i++)
        {
            for(int j = 0; j < cors1[0].size(); j++)
            {
                cors1[i][j]= *it++;

                //cout<<setprecision(12)<<cors1[i][j]<<' ';

            }
                //cout<<endl;
        }
    }


    // 计算 kT
    double kT = k_boltz * temperature;

    // 接受移动的总个数
    int naccept = 0;

    // 拒绝移动的总个数
    int nreject = 0;

    // 打印初始PDB文件
    print_pdb(cors1, n_atoms, 0);

    double total_energy_sum = 0.0;
    double total_energy_avg = 0.0;
    double move_ratio = 0.0;
    double naccept_ratio = 0.0;
    double move_cal = 0.0;

    //执行随机移动（第一种移动）
    for (int move=1; move<=num_moves; move = move + 1)
    {
         // 计算旧能量
        const double old_energy = calculate_energy(cors1, n_atoms, box_size, sigma, epsilon, r_cut);

         // 随机选择原子
         int atom = int( rand(0, n_atoms) );

        // 保存旧坐标
        copy_coordinates(cors1, old_cors1);

        // 执行随机位移 - 通过delta值在三个维度执行位移
        const double delta_x = rand(-max_translate, max_translate);
        const double delta_y = rand(-max_translate, max_translate);
        const double delta_z = rand(-max_translate, max_translate);

        cors1[atom][0] = cors1[atom][0] + delta_x;
        cors1[atom][1] = cors1[atom][1] + delta_y;
        cors1[atom][2] = cors1[atom][2] + delta_z;

        // 保证坐标在盒子内
        cors1[atom][0] = wrap_into_box(cors1[atom][0], box_size[0]);
        cors1[atom][1] = wrap_into_box(cors1[atom][1], box_size[1]);
        cors1[atom][2] = wrap_into_box(cors1[atom][2], box_size[2]);

        // 计算新能量
        const double new_energy = calculate_energy(cors1, n_atoms, box_size, sigma, epsilon, r_cut);

        bool accept = false;

        // 如果能量降低，则自动接受
        if (new_energy <= old_energy)
        {
            accept = true;
        }
        else
        {
            // 执行Metropolis能量判断 - 比较
            // exp( E_new - E_old / kT ) >= rand(0,1)
            const double x = exp( -(new_energy - old_energy) / kT );

            if (x >= rand(0.0, 1.0))
            {
                accept = true;
            }
            else
            {
                accept = false;
            }

        }

        double total_energy = 0.0;


        if (accept)
        {
            // 接受移动
            naccept = naccept + 1;
            total_energy = new_energy;
        }
        else
        {
            // 拒绝移动 - 恢复旧构型
            nreject = nreject + 1;


            copy_coordinates( old_cors1, cors1 );

            total_energy = old_energy;
        }

        //平均值
        /*
        if(move <= 5E2)
            continue;
        else
        {
            total_energy_sum += total_energy;

            total_energy_avg = total_energy_sum/(move-5e2);
        }
        */

        total_energy_sum += total_energy;
        total_energy_avg = total_energy_sum/(naccept+nreject);

        naccept_ratio = (double)naccept;

        move_cal = (double)move;

        move_ratio = naccept_ratio/move_cal;

        /*
        if(move_ratio > 0.55)
                max_translate = max_translate * 1.05;

        else if(move_ratio <0.45)
                max_translate = max_translate * 0.95;
        */
        // print the energy every 1000 moves每间隔打印一次能量、接受拒绝情况
        if (move % 1000 == 0)
        {
            printf("%d %e  %e   %d  %d %.2lf\n", move, total_energy, total_energy_avg, naccept, nreject,move_ratio);
           // printf("%d %e  %e   %d  %d %.2lf\n", move, total_energy, total_energy_avg, naccept, nreject,move_ratio);
        }
        // print the coordinates every 10000 moves

       /* if (move % 10000 == 0)
        {
             print_pdb(cors1, n_atoms, move);
        }
        */
    }


    return 0;
}

