//
// Created by 张力文 on 2022/4/21.
//
/*
 * 用于在程序中添加约化单位以测试Lennard-Jones基准值的程序
 * reference：chryswoods metropolis.cpp
 * 构型不是提前输入，而是每次程序运行时随机生成的
 * 原始程序中的原子数和盒子尺寸是给定的，需要做一定的修改
 * 输出：每1000次移动输出move, total_energy, naccept, nreject
 *      每10000次移动输出一次pdf文件
 *
 * 2022.5.10新修改：
 * 将输入坐标从二维vector改回静态数组，并且将其大小初始化为1000
 * 连同其他函数的形参类型也改回二维数组
 *
 * 2022.5.11新修改：
 * 约化单位更改：
 * （1）使用的是Ar的epsilon和sigma值
 * （2）对于给定初始构型的case（如N=30，L=8，根据LJ MC网站的基准值表格，调整L=11）
 * （3）根据T* = 0.85， rho* = 0.009，得出T=120
 * （4）加入截断判断r_cut = 3*sigma （os：前后的能量可视化后并看不到很大的差别……）
 * （5）能量的平均值计算
 *
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

using namespace std;

// Set the number of atoms in the box
const int n_atoms = 30;

// Set the number of Monte Carlo moves to perform
const int num_moves = 500000;

// Set the size of the box (in Angstroms)
//const double box_size[3] = { 15.0, 15.0, 15.0 };
const double box_size[3] = { 11.0, 11.0, 11.0 };

// The maximum amount that the atom can be translated by
const double max_translate = 0.5;  // angstroms

// Simulation temperature
//const double temperature = 298.15;   // kelvin
const double temperature = 102;

// Give the Lennard Jones parameters for the atoms
// (these are the OPLS parameters for Krypton)
//const double sigma = 3.624;     // angstroms
//const double epsilon = 0.317;   // kcal mol-1

const double k_boltz = 1.987206504191549E-003;

//Ar    T* = 0.85   rho*=0.009

const double sigma = 3.41;
const double epsilon_boltz = 120;
const double epsilon = epsilon_boltz*k_boltz;
const double r_cut = 3 * sigma;

// function to return a random number between 'start' to 'end'
double rand(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}

// Subroutine to apply periodic boundaries
double make_periodic(double x, const double box)
{
    while (x < -0.5*box)
    {
        x = x + box;
    }

    while (x > 0.5*box)
    {
        x = x - box;
    }

    return x;
}

// Subroutine to wrap the coordinates into a box
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

// Subroutine to print a PDB of the coordinates
void print_pdb(double coords[][3] , const int n_atoms, const int move)
//void print_pdb(vector<vector<double>> &coords, const int n_atoms, const int move)
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

// Subroutine that calculates the energies of the atoms
double calculate_energy(double coords[][3], const int n_atoms, const double *box_size,
                        const double sigma, const double epsilon,const double r_cut)
//double calculate_energy(vector<vector<double>>& coords, const int n_atoms, const double *box_size,
//                        const double sigma, const double epsilon)
{
    // Loop over all pairs of atoms and calculate
    // the LJ energy
    double total_energy = 0;
    double r_cut_box = r_cut / box_size[0];
    double r_cut__box_sq = r_cut_box * r_cut_box;
    double r_cut_sq = r_cut * r_cut;

    for (int i = 0; i < n_atoms-1; i = i + 1)
    {
        for (int j = i+1; j < n_atoms; j = j + 1)
        {
            double delta_x = coords[j][0] - coords[i][0];
            double delta_y = coords[j][1] - coords[i][1];
            double delta_z = coords[j][2] - coords[i][2];

            // Apply periodic boundaries
            delta_x = make_periodic(delta_x, box_size[0]);
            delta_y = make_periodic(delta_y, box_size[1]);
            delta_z = make_periodic(delta_z, box_size[2]);

            //const double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
            //                  (delta_z*delta_z);
            double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                              (delta_z*delta_z);
            if(r2 < r_cut_sq)
            {
                // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
                const double sig2_over_r2 = (sigma*sigma) / r2;
                const double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
                const double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

                //const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );

                //约化单位之后的能量计算公式
                const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );
                total_energy = total_energy + e_lj;
            }

//            // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
//            const double sig2_over_r2 = (sigma*sigma) / r2;
//            const double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
//            const double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;
//
//            const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );
//
//            //约化单位之后的能量计算公式
//            //const double e_lj = 4.0 * ( sig12_over_r12 - sig6_over_r6 );
//            total_energy = total_energy + e_lj;
        }
    }

    // return the total energy of the atoms
    return total_energy * epsilon;
}

void copy_coordinates(double from[][3], double to[][3])
//void copy_coordinates(vector<vector<double>>& from, vector<vector<double>>& to)
{
    for (int i=0; i<n_atoms; ++i)
    {
        to[i][0] = from[i][0];
        to[i][1] = from[i][1];
        to[i][2] = from[i][2];
    }
}

int CountLines(char *filename)
{
    ifstream ReadFile;
    int n=0;
    string tmp;
    ReadFile.open(filename,ios::in);//ios::in 表示以只读的方式读取文件
    if(ReadFile.fail())//文件打开失败:返回0
    {
        return 0;
    }
    else//文件存在
    {
        while(getline(ReadFile,tmp,'\n'))
        {
            n++;
        }
        ReadFile.close();
        return n;
    }
}


int main(int argc, const char **argv)
{


    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

    double coords[1000][3] = {0.0};
    //double **coords1 = new double*[1000];

    double old_coords[1000][3] = {0.0};
    //double **old_coords1 = new double*[1000];

    //vector<vector<double>> coords(n_atoms, vector<double>(3));
    //vector<vector<double>> old_coords(n_atoms, vector<double>(3));

    /*
    for (int i = 0; i < n_atoms; i = i + 1)
    {
        coords1[i] = new double[3];
        old_coords1[i] = new double[3];
    }
     */

    ifstream file;
    int LINES;
    char filename[512]="input_box.txt";
    file.open(filename,ios::in);

    if(file.fail())
    {
        cout<<"文件不存在."<<endl;
        file.close();
    }
    else//文件存在
    {
        LINES = CountLines(filename);
        cout << "The number of lines is:" << LINES<< endl;

        int cols = 3;

        while (!file.eof()) //读取数据到数组
        {
            for (int i = 0; i < LINES; ++i)
            {
                for (int j = 0; j < cols; j++)
                {
                    file >> coords[i][j];
                    //cout << coords[i][j] << ' ';
                }
                //cout << endl;
            }
        }
        file.close(); //关闭文件
    }


    // Randomly generate the coordinates of the atoms in the box
    /*
    for (int i = 0; i < n_atoms; i = i + 1)
    {
       // coords[i] = new double[3];
        //old_coords[i] = new double[3];

        // Note "rand(0,x)" would generate a random number
        // between 0 and $x
        coords[i][0] = rand(0, box_size[0]);
        coords[i][1] = rand(0, box_size[1]);
        coords[i][2] = rand(0, box_size[2]);
    }
     */

    // calculate kT
    //const double k_boltz = 1.987206504191549E-003;  // kcal mol-1 K-1

    const double kT = k_boltz * temperature;

    // The total number of accepted moves
    int naccept = 0;

    // The total number of rejected moves
    int nreject = 0;

    // Print the initial PDB file
    print_pdb(coords, n_atoms, 0);

    for (int move=1; move<=num_moves; move = move + 1)
    {
        // calculate the old energy
        const double old_energy = calculate_energy(coords, n_atoms, box_size, sigma, epsilon,r_cut);

        // Pick a random atom
        int atom = int( rand(0, n_atoms) );

        // save the old coordinates
        copy_coordinates(coords, old_coords);

        // Make the move - translate by a delta in each dimension
        const double delta_x = rand(-max_translate, max_translate);
        const double delta_y = rand(-max_translate, max_translate);
        const double delta_z = rand(-max_translate, max_translate);

        coords[atom][0] = coords[atom][0] + delta_x;
        coords[atom][1] = coords[atom][1] + delta_y;
        coords[atom][2] = coords[atom][2] + delta_z;

        // wrap the coordinates back into the box
        coords[atom][0] = wrap_into_box(coords[atom][0], box_size[0]);
        coords[atom][1] = wrap_into_box(coords[atom][1], box_size[1]);
        coords[atom][2] = wrap_into_box(coords[atom][2], box_size[2]);

        // calculate the new energy
        const double new_energy = calculate_energy(coords, n_atoms, box_size, sigma, epsilon,r_cut);

        bool accept = false;

        // Automatically accept if the energy goes down
        if (new_energy <= old_energy)
        {
            accept = true;
        }
        else
        {
            // Now apply the Monte Carlo test - compare
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

        double total_energy = 0;
        double total_energy_avg = 0;
        double total_energy_sum = 0;

        if (accept)
        {
            // accept the move
            naccept = naccept + 1;
            total_energy = new_energy;
        }
        else
        {
            // reject the move - restore the old coordinates
            nreject = nreject + 1;

            // restore the old coordinates
            copy_coordinates( old_coords, coords );

            total_energy = old_energy;
        }

        //能量平均值计算
        if(move > 5000)
        {
            total_energy_sum += total_energy;
            total_energy_avg = total_energy_sum / (move - 5000);
        }

        // print the energy every 1000 moves
        if (move % 1000 == 0)
        {
            printf("%d %f  %f %d  %d\n", move, total_energy, total_energy_avg, naccept, nreject);
        }

        // print the coordinates every 10000 moves
        /*
        if (move % 10000 == 0)
        {
            print_pdb(coords, n_atoms, move);
        }
         */
    }


    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();

    chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t2 - t1);

    cout << "run time = " << time_used.count() << " seconds. " << endl;


    return 0;
}
