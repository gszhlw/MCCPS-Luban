//
// Created by 张力文 on 2021/11/25.
//
#include <iostream>

using namespace  std;
/*
 * 函数名：monte_carlo
 * 功能：执行蒙特卡洛（MC）模拟的主要功能
 */

int nstep;



for(int n = 0; n < nstep; n++)
{
    for(int i = 0; i < inerctcle; i++)
{
    //实行一个蒙特卡洛移动
    luban_mc_move();

    //累计平均值

    luban_accumulate_averages();

}

void luban_mc_move()
{

}


    //计算总能量
    luban_total_energy();


}

//崩溃恢复，对应Cassendra中的.chk文件的作用，跳过初始化/平衡，直接跳转到生产运行中写入崩溃文件之后的位置

//给采样程序分配内存

for(int i = 0; i < number_of_steps; i++)
{
    //设置选择系统

    //随机选择一个组分

    //随机在已经选择的概率中选择任意MC移动，其实就是以随机顺序执行各种移动
    translation_move();
    rotation_move();
    reinsertion_move();
    identity_change_move();
    volume_move();

    //每个循环，优化接受率

    //优化能量，计算所有系统的总能量



}

/*
 * gemc_control
 * 1.模拟几种组分
 * 2.载入盒子形状，盒子个数，盒子类型，计算盒子的不同性质，包括体积等
 * 3.确定VDW类型，电荷作用模型，以及相关参数，VDW混合规则
 * 4.载入分子连接性以及力场参数，注意在此之前必须已经知道组分个数了
 * 5.如果模拟包括一个用于转移的中转盒子，从那个盒子里读取理想气体分子个数
 *
 * 确定unique atom类型的个数和类别，建立一个VDW相互作用
 *
 * 创建分子内非键标度数组
 *
 * 启动类型
 *
 * 种子信息
 *
 * 获得模拟温度
 *
 * 读入所有移动的概率
 *
 * 获得CBMC信息
 *
 *
 *
 *
 *
 */