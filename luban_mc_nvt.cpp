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
 * （5）能量的平均值计算(计算出的能量平均值为-4.563448，与NIST上的-6.2391对不上，怀疑是差在截断势能)
 *
 * 2022.5.12新修改：
 * 根据python example代码中有关势能截断的函数实现添加新的计算截断势能函数potential lrc，每次打印的截断能量大小为-0.001667，对总能量的影响微乎其微
 *
 * 2022.5.13新修改：
 * 模拟原子仍为Ar，根据T*=0.85， rho*=0.9,推出L=28
 * (1)输入构型原子个数改为500
 *（2）构型坐标随机生成
 *
 * 2022.5.18新修改
 * 约化单位更改：
 * （1）使用的是Ne的epsilon和sigma值
 * （2）对于给定初始构型的case（如N=30，L=8，根据LJ MC网站的基准值表格，调整L=8.9）
 * （3）根据T* = 0.85， rho* = 0.009，得出T=30.6
 *  (4)得出平均能量值为0.303973，与NIST仍旧对不上
 *
 * 2022.5.19
 *（1）使用NIST的给定构型坐标（N = 400）
 *（2）根据T* = 0.85， rho* = 0.009，T=120，L=26
 *
 * 2022.5.20
 * (1)添加了计算压力的函数
 *（2）模拟参数为 N=400，对应 T* = 0.85， rho* = 0.82，U*=- E+00， p* = 5.5355E-01，
 *
 * 2022.5.21
 * (1) 8.50E-01	7.00E-03	-7.2834E-02	-- -36.417   	5.6565E-03-- 2.82	L = 141
 * (2)8.50E-01	5.00E-03	-5.1901E-02	-- -25.9  7.53E-05	4.1003E-03	5.05E-07
 *
 * 2022.5.23
 * （1）气态密度下的能量压力不收敛问题（也可能液态也没收敛，只是不明显）
 * （2）可能原因：能量计算过程中迭代值是否清零
 *              尾部、长程校正的问题
 * （3）仍然使用最初的测试数据
 *      8.50E-01	9.00E-01--L=28	-6.2391E+00 -- -3,119.55	5.27E-03	2.2314E+00	2.72E-02 -- 1,115.7
 *
 * 2022.5.23
 * (1) 计算平均能量的过程需要添加ideal gas contribution = 1.5 * temperature
 * (2) 计算压力
 *
 * 2022.5.25
 * （1）最后一组测试数据 T* = 8.50E-01	rho* = 7.80E-01	U* = -5.5386E+00	7.26E-04	p* = 4.7924E-02	3.18E-03
 *
 * 2022.5.30
 *  （1）加入体积改变移动，能量对不上
 *  20000000 -9.600476e-01  -1.963523e+00 4417621  15582379
 *  Box size = (10.144485,10.144485,10.144485). Volume = 1043.974918 A^3
 *  （2）不知道是不正确，还是npt模拟结果不能跟nist比对（应该不是。。。）
 *
 */

#include "luban_mc_nvt.h"

using namespace std;




int main(int argc, const char **argv)
{
    //chrono相关用于计算执行模拟前后花费的时间
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

    //如果需要读入外部构型文件，则需要执行下列代码
    /*
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
*/

    // 如果需要随机生成构型坐标，则使用下列代码 Randomly generate the coordinates of the atoms in the box
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


    // calculate kT
    //const double k_boltz = 1.987206504191549E-003;  // kcal mol-1 K-1

    const double kT = k_boltz * temperature;

    // The total number of accepted moves
    int naccept = 0;

    // The total number of rejected moves
    int nreject = 0;

    // Print the initial PDB file
    print_pdb(coords, n_atoms, 0);


    double total_energy_avg = 0;
    double total_energy_sum = 0;
    double total_pressure_avg = 0;
    double total_pressure_sum = 0;

    double move_cal = 0;


    for (int move=1; move<=num_moves; move = move + 1)
    {
        // calculate the old energy
        const double old_energy = calculate_energy(coords, n_atoms, box_size, sigma, epsilon,r_cut);
        const double old_energy_lrc = calculate_lrc_energy(coords, n_atoms, box_size, sigma, epsilon,r_cut);

        const double old_pressure = calculate_pressure(coords, n_atoms, box_size, sigma, epsilon,r_cut);
        const double old_pressure_lrc = calculate_lrc_pressure(coords, n_atoms, box_size, sigma, epsilon,r_cut);
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
        const double new_energy_lrc = calculate_lrc_energy(coords, n_atoms, box_size, sigma, epsilon,r_cut);

        const double new_pressure = calculate_pressure(coords, n_atoms, box_size, sigma, epsilon,r_cut);
        const double new_pressure_lrc = calculate_lrc_pressure(coords, n_atoms, box_size, sigma, epsilon,r_cut);
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
        double total_energy_lrc = 0;
        double total_pressure = 0;
        double total_pressure_lrc = 0;

        if (accept)
        {
            // accept the move
            naccept = naccept + 1;
            total_energy = new_energy;
            total_energy_lrc = new_energy_lrc;
            total_pressure = new_pressure;
            total_pressure_lrc = new_pressure_lrc;
        }
        else
        {
            // reject the move - restore the old coordinates
            nreject = nreject + 1;

            // restore the old coordinates
            copy_coordinates( old_coords, coords );

            total_energy = old_energy;
            total_energy_lrc = old_energy_lrc;
            total_pressure = old_pressure;
            total_pressure_lrc = old_pressure_lrc;
        }

        //能量平均值计算，5000为可视化后确定的能量开始收敛的step
        if(move > 100000)
        {
          //  count++;
            move_cal = 1.0 * (move - 100000);
            total_energy_sum += total_energy;
            total_pressure_sum += total_pressure;
          //  cout<<"count: "<<count;
          //  cout<<"     total_energy_sum"<<total_energy_sum;
            total_energy_avg = total_energy_sum /move_cal;
            total_pressure_avg = total_pressure_sum/move_cal/(box_size[0] * box_size[1] * box_size[2]);

        //    cout<<"     total_energy_avg"<<total_energy_avg;
        //    cout<<"     move - 5000:"<<move_cal<<endl;
        }

        const double ideal_gas_contribution_energy = 1.5 * temperature;
        const double ideal_gas_contribution_pressure = n_atoms / (box_size[0] * box_size[1] * box_size[2]) * temperature;

        // print the energy every 1000 moves

        if (move % 1000 == 0)
        {
            printf("%d %f  %f %f  %f %f %f %d  %d\n", move, total_energy + ideal_gas_contribution_energy, total_energy_avg + ideal_gas_contribution_energy,total_energy_lrc,total_pressure + ideal_gas_contribution_pressure, total_pressure_avg +ideal_gas_contribution_pressure,total_pressure_lrc,naccept, nreject);
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
