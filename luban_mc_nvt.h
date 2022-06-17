//
// Created by Administrator on 2022/6/9.
//

#ifndef MCCPS_LUBAN_LUBAN_MC_NVT_H
#define MCCPS_LUBAN_LUBAN_MC_NVT_H
//#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
using namespace  std;
//#include </usr/local/Cellar/libomp/14.0.0/include/omp.h>

// Set the number of atoms in the box
//const int n_atoms = 30;
const int n_atoms = 32;
// Set the number of Monte Carlo moves to perform
const int num_moves = 500000;

// Set the size of the box (in Angstroms)
//const double box_size[3] = { 15.0, 15.0, 15.0 };
const double box_size[3] = { 15.0, 15.0, 15.0 };


// The maximum amount that the atom can be translated by
const double max_translate = 0.4;  // angstroms

// Simulation temperature
//const double temperature = 298.15;   // kelvin
const double temperature = 120;

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

    //用盒子尺寸L去判断 x+L x-L x， 返回abs的min

    while (x < -0.5*box)
    {
        x = x + box;
    }

    while (x > 0.5*box)
    {
        x = x - box;
    }
    /*
   double x1 = x + box;
   double x2 = x - box;

   double x1_abs = fabs(x1);
   double x2_abs = fabs(x2);
   double x_abs = fabs(x);

    double min_pre = (x1_abs < x2_abs) ? x1_abs : x2_abs;
    double min = (min_pre < x_abs) ? min_pre : x_abs;

    return min;
     */
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
    double r_cut_box_sq = r_cut_box * r_cut_box;
    double r_cut_sq = r_cut * r_cut;
    double box_sq = box_size[0] * box_size[0];
    double delta_x;
    double delta_y;
    double delta_z;
    double r2;
    double sig2_over_r2;
    double sig6_over_r6;
    double sig12_over_r12;
    double e_lj;

     #pragma omp parallel
    {
         #pragma omp single
        //#pragma omp parallel for collapse(2) shared(n_atoms) private(j)

        {
            for (int i = 0; i < n_atoms-1; i = i + 1)
            {
                //#pragma omp task
               // #pragma omp parallel for collapse(2)
                for (int j = i+1; j < n_atoms; j = j + 1)
                {
                    #pragma omp task
                    {
                        delta_x = coords[j][0] - coords[i][0];
                        delta_y = coords[j][1] - coords[i][1];
                        delta_z = coords[j][2] - coords[i][2];

                        // Apply periodic boundaries
                        //PBC判断 9个距离值的最小值
                        delta_x = make_periodic(delta_x, box_size[0]);
                        delta_y = make_periodic(delta_y, box_size[1]);
                        delta_z = make_periodic(delta_z, box_size[2]);

                        //const double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                        //                  (delta_z*delta_z);
                        r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                                    (delta_z*delta_z);
                        if((r2/box_sq)< r_cut_box_sq)
                        {
                            // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
                            sig2_over_r2 = (sigma*sigma) / r2;
                            sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
                            sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

                            //const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );

                            //约化单位之后的能量计算公式
                             e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );
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
            }
        }

    }
    // return the total energy of the atoms
    return total_energy / epsilon;
}
double calculate_pressure(double coords[][3], const int n_atoms, const double *box_size,
                          const double sigma, const double epsilon,const double r_cut)
{
    // Loop over all pairs of atoms and calculate
    // the LJ pressure
    double vir = 0;

    double r_cut_box = r_cut / box_size[0];
    double r_cut_box_sq = r_cut_box * r_cut_box;
    double r_cut_sq = r_cut * r_cut;
    double box_sq = box_size[0] * box_size[0];
    double delta_x;
    double delta_y;
    double delta_z;
    double r2;
    double sig2_over_r2;
    double sig6_over_r6;
    double sig12_over_r12;


     #pragma omp parallel
    {
        #pragma omp single
        //#pragma omp parallel for collapse(2) shared(n_atoms) private(j)

        {
            for (int i = 0; i < n_atoms-1; i = i + 1)
            {
               // #pragma omp task
               // #pragma omp parallel for collapse(2)
                for (int j = i+1; j < n_atoms; j = j + 1)
                {
                    #pragma omp task
                    {
                        delta_x = coords[j][0] - coords[i][0];
                        delta_y = coords[j][1] - coords[i][1];
                        delta_z = coords[j][2] - coords[i][2];

                        // Apply periodic boundaries

                        delta_x = make_periodic(delta_x, box_size[0]);
                        delta_y = make_periodic(delta_y, box_size[1]);
                        delta_z = make_periodic(delta_z, box_size[2]);

                        //const double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                        //                  (delta_z*delta_z);
                        r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                                    (delta_z*delta_z);
                        if((r2/box_sq) < r_cut_box_sq)
                        {
                            // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
                            sig2_over_r2 = (sigma*sigma) / r2;
                            sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
                            sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

                            //const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );

                            //约化单位之后的能量计算公式
                            //const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );
                            vir = -(2 * sig12_over_r12 - sig6_over_r6) * 24.0 * epsilon / 3.0;
                            //total_energy = total_energy + e_lj;
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
            }
        }



    }


    // return the total energy of the atoms
    // return total_energy / epsilon;
    return vir * sigma *sigma *sigma / epsilon;

}
double calculate_lrc_energy(double coords[][3], const int n_atoms, const double *box_size,
                            const double sigma, const double epsilon,const double r_cut)
{
    /*
     * sr3 = sigma / r_cut ** 3
    return math.pi * ( (8.0/9.0)  * sr3**3  - (8.0/3.0)  * sr3 ) * density
     */
    double potential_lrc = 0;
    double sr3 = (sigma / r_cut)*(sigma / r_cut)*(sigma / r_cut);
    potential_lrc = M_PI * ((8.0/9.0) * sr3 * sr3 * sr3 - (8.0/3.0)  * sr3) *n_atoms / (box_size[0] *box_size[1] * box_size[2]);

    return potential_lrc / epsilon;
}


double calculate_lrc_pressure(double coords[][3], const int n_atoms, const double *box_size,
                              const double sigma, const double epsilon,const double r_cut)
{
    double pressure_lrc = 0;
    double sr3 = (sigma / r_cut)*(sigma / r_cut)*(sigma / r_cut);
    pressure_lrc = M_PI * ((32.0/9.0) * sr3 * sr3 * sr3 - (16.0/3.0)  * sr3) *((n_atoms /(box_size[0] *box_size[1] * box_size[2])) * (n_atoms /(box_size[0] *box_size[1] * box_size[2])));
    // math.pi * ( (32.0/9.0) * sr3**3  - (16.0/3.0) * sr3 ) * density**2
    return pressure_lrc * sigma *sigma *sigma / epsilon;
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

#endif //MCCPS_LUBAN_LUBAN_MC_NVT_H
