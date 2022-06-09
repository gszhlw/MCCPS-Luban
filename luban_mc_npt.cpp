//
// Created by 张力文 on 2022/4/21.
//
/*
 * 2022.5.31
 * (1)NPT模拟需要给定压力值：0.01341956808
 *（2）8.50E-01	9.00E-011.5	5.27E-03	2.2314E+00	2.72E-02
 * (3)标准的box：1,101.4394724361
 */
#include "luban_mc_npt.h"

using namespace std;


int main(int argc, const char **argv)
{
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
   // double **coords = new double*[n_atoms];
   // double **old_coords = new double*[n_atoms];
    double coords[1000][3] = {0.0};
    double old_coords[1000][3] = {0.0};

    //vector<vector<double>> coords(n_atoms, vector<double>(3));
    //vector<vector<double>> old_coords(n_atoms, vector<double>(3));


    // Randomly generate the coordinates of the atoms in the box
    for (int i = 0; i < n_atoms; i = i + 1)
    {
        //coords[i] = new double[3];
        //old_coords[i] = new double[3];

        // Note "rand(0,x)" would generate a random number
        // between 0 and $x
        coords[i][0] = rand(0, box_size[0]);
        coords[i][1] = rand(0, box_size[1]);
        coords[i][2] = rand(0, box_size[2]);
    }



    const double kT = k_boltz * temperature;

    // The total number of accepted moves
    int naccept = 0;

    // The total number of rejected moves
    int nreject = 0;

    // The total number of accepted volume moves
    int nvolaccept = 0;

    // The total number of rejected volume moves
    int nvolreject = 0;

    // Print the initial PDB file
    print_pdb(coords, n_atoms, 0);


    double total_energy_avg = 0;
    double total_energy_sum = 0;

    double move_cal = 0;

    for (int move=1; move<=num_moves; move = move + 1)
    {
        // calculate the old energy
        const double old_energy = calculate_energy(coords, n_atoms, box_size, sigma, epsilon,r_cut);

        // calculate the old volume of the box
        const double V_old = box_size[0] * box_size[1] * box_size[2];

        // Pick a random atom
        int atom = int( rand(0, n_atoms) );

        // save the old coordinates
        copy_coordinates(coords, old_coords);

        // save the old box dimensions
        const double old_box_size[3] = { box_size[0], box_size[1], box_size[2] };

        // Decide if we are performing an atom move, or a volume move
        if (rand(0.0, 1.0) <= 1.0 / n_atoms)
        {
            // 1 in $n_atoms chance of being here. Perform a volume move
            // by changing the volume for a random amount
            double delta_vol = rand(-max_volume_change, max_volume_change);

            double V_new = V_old + delta_vol;

            // Volume is the cube of the box length, so add the cube root
            // of this change onto the box size
            double box_side = pow(V_new, 1.0/3.0);

            // work out how much we need to scale the position of the atoms
            double scale_ratio = std::pow( V_new / V_old, 1.0/3.0 );

            // now translate every atom so that it is scaled from the center
            // of the box
            for (int i=0; i<n_atoms; ++i)
            {
                double dx = coords[i][0] - (0.5*box_size[0]);
                double dy = coords[i][1] - (0.5*box_size[1]);
                double dz = coords[i][2] - (0.5*box_size[2]);

                double length = std::sqrt(dx*dx + dy*dy + dz*dz);

                if (length > 0.01)   // don't scale atoms already near the center
                {
                    dx /= length;
                    dy /= length;
                    dz /= length;

                    length *= scale_ratio;

                    coords[i][0] = (0.5*box_size[0]) + dx * length;
                    coords[i][1] = (0.5*box_size[1]) + dy * length;
                    coords[i][2] = (0.5*box_size[2]) + dz * length;
                }
            }

            // now update the new size of the box
            box_size[0] = box_side;
            box_size[1] = box_side;
            box_size[2] = box_side;
        }
        else
        {
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
        }

        // calculate the new energy
        const double new_energy = calculate_energy(coords, n_atoms, box_size, sigma, epsilon,r_cut);

        // calculate the new volume of the box
        const double V_new = box_size[0] * box_size[1] * box_size[2];

        bool accept = false;

        // Automatically accept if the energy goes down
        if (new_energy <= old_energy)
        {
            accept = true;
        }
        else
        {
            // Now apply the Monte Carlo test - compare
                // exp( -(E_new - E_old + P(V_new - V_old)) / kT
                //             +  N ln (V_new - V_old) ) >= rand(0,1)
            double x = exp( -((new_energy - old_energy + pressure * (V_new - V_old)) / kT)
                            + (n_atoms * (log(V_new) - log(V_old) )) );

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

            box_size[0] = old_box_size[0];
            box_size[1] = old_box_size[1];
            box_size[2] = old_box_size[2];

            // restore the old coordinates
            copy_coordinates( old_coords, coords );

            total_energy = old_energy;
        }


        //能量平均值计算，5000为可视化后确定的能量开始收敛的step
        if(move > 100000)
        {
            //  count++;
            move_cal = 1.0 * (move - 100000);
            total_energy_sum += total_energy;

            //  cout<<"count: "<<count;
            //  cout<<"     total_energy_sum"<<total_energy_sum;
            total_energy_avg = total_energy_sum /move_cal;

        }

        //理想气体贡献部分
        const double ideal_gas_contribution_energy = 1.5 * temperature;


        // print the energy every 50000 moves
        if (move % 5000 == 0)
        {
            printf("%d %e  %e %d  %d ", move, total_energy+ideal_gas_contribution_energy, total_energy_avg+ideal_gas_contribution_energy,
                   naccept, nreject);

            const double vol = box_size[0] * box_size[1] * box_size[2];

            printf("    Box size = (%f,%f,%f). Volume = %f A^3\n",
                   box_size[0], box_size[1], box_size[2], vol);
        }

        // print the coordinates every 500000 moves
      /*
        if (move % 500000 == 0)
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