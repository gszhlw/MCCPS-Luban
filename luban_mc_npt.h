//
// Created by Administrator on 2022/6/9.
//

#ifndef MCCPS_LUBAN_LUBAN_MC_NPT_H
#define MCCPS_LUBAN_LUBAN_MC_NPT_H
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <iostream>
#include <vector>
// Set the number of atoms in the box
const int n_atoms = 500;

// Set the number of Monte Carlo moves to perform
const int num_moves = 20000000;

// Set the size of the box (in Angstroms)
double box_size[3] = { 60.0, 60.0, 60.0 };

// The maximum amount that the atom can be translated by
const double max_translate = 0.2;   // angstroms

// Simulation temperature
const double temperature = 102;   // kelvin

// Simulation pressure (atmospheres converted to internal
// units - kcal mol-1 A-3)
//double pressure = 1 * 1.458397506863647E-005;   // atmospheres
double pressure = 0.01341956808;
// The maximum amount to change the volume - the
// best value is 10% of the number of atoms
const double max_volume_change = 0.1 * n_atoms;   // Angstroms**3

// Give the Lennard Jones parameters for the atoms
// (these are the OPLS parameters for Krypton)

// calculate kT
const double k_boltz = 1.987206504191549E-003;  // kcal mol-1 K-1

const double sigma = 3.41;     // angstroms
const double epsilon = 120 * k_boltz;   // kcal mol-1

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
//void print_pdb(double **coords, const int n_atoms, const int move)
void print_pdb(double coords[][3], const int n_atoms, const int move)
{
    char filename[128];

    snprintf(filename, 128, "output%00000008d.pdb", move);

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
//double calculate_energy(double **coords, const int n_atoms, const double *box_size,
//                        const double sigma, const double epsilon)
/*
double calculate_energy(vector<vector<double>>& coords, const int n_atoms, const double *box_size,
                        const double sigma, const double epsilon)
{
    // Loop over all pairs of atoms and calculate
    // the LJ energy
    double total_energy = 0;

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

            const double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                              (delta_z*delta_z);

            // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
            const double sig2_over_r2 = (sigma*sigma) / r2;
            const double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
            const double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

            const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );

            total_energy = total_energy + e_lj;
        }
    }

    // return the total energy of the atoms
    return total_energy;
}
*/
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

    // #pragma omp parallel
    {
        // #pragma omp single
        {
            for (int i = 0; i < n_atoms-1; i = i + 1)
            {
                for (int j = i+1; j < n_atoms; j = j + 1)
                {
                    //#pragma omp task
                    {


                        double delta_x = coords[j][0] - coords[i][0];
                        double delta_y = coords[j][1] - coords[i][1];
                        double delta_z = coords[j][2] - coords[i][2];

                        // Apply periodic boundaries
                        //PBC判断 9个距离值的最小值
                        delta_x = make_periodic(delta_x, box_size[0]);
                        delta_y = make_periodic(delta_y, box_size[1]);
                        delta_z = make_periodic(delta_z, box_size[2]);

                        //const double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                        //                  (delta_z*delta_z);
                        double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                                    (delta_z*delta_z);
                        if((r2/box_sq)< r_cut_box_sq)
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
            }
        }

    }
    // return the total energy of the atoms
    return total_energy / epsilon;
}
void copy_coordinates(double from[][3], double to[][3])
{
    for (int i=0; i<n_atoms; ++i)
    {
        to[i][0] = from[i][0];
        to[i][1] = from[i][1];
        to[i][2] = from[i][2];
    }
}

#endif //MCCPS_LUBAN_LUBAN_MC_NPT_H
