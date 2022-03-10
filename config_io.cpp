//
// Created by 张力文 on 2021/12/8.
//
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "config_io.h"

using namespace std;
#define Kb 1.3806488e-23
//用来实现原子、分子构型I/O的例程

int main()
{
    double temperature = 1.0;
    double epsilon = Kb*temperature/0.85;
    cout<<"epsilon:"<<epsilon<<endl;
    return 0;
}

