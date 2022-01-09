//
// Created by 张力文 on 2021/12/8.
//
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "config_io.h"

using namespace std;

//用来实现原子、分子构型I/O的例程

void read_cnf_atoms()
{
    ifstream myfile("/Users/zlw/ClionProjects/MCCPS Luban/hello.txt");
    string temp;
    string para[5]={"ensemble","temperature", "nmolty", "nummols", "numboxes"};
    string ensemble;
    double temperature;
    int nmolty, nummols, numboxes;

    if(!myfile.is_open())
    {
        cout<<"file open error"<<endl;
    }

    int i = 0;//全局变量要小心


    for(string temp; getline(myfile,temp);)
    {
        //outfile<<temp<<endl;
        //由于整行整行读入，读入道str中时，文件中的每个换行符都被丢掉了，为了照原样输出，在out流上输出时需要再补上一个回车<<endl
        int res = temp.compare(para[i++]);
        //下面三个输出用来测试
        /*
        cout<<temp<<endl;//temperature
        cout<<para[i]<<endl;//ensemble
        cout<<"res = "<<res<<endl;
        */

        if(res == 0)
        {
            cout<<temp<<":"<<endl;

            getline(myfile,temp);

           //将字符串参数传回main函数


           cout<<temp<<endl;
        }
        else
        {
            cout<<"parameter error"<<endl;
        }

        //如果是可选参数，不能傻乎乎的比对字符串读取下一行，需要判断有这个可选参数对号入座读下一行，如果没有这个参数就跳过


    }

    myfile.close();


    return;
}

