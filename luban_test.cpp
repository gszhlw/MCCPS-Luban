//
// Created by 张力文 on 2022/1/13.
//

/*
 * 输入1：运行脚本run.inp （先默认手动输入，因为涉及的参数并不多）
 * 输入2：构型文件cnf.inp，因为含有多个原子（256个）的构型坐标，需要重点处理
 * 原始的输入格式需要改进，如果要求输入.xyz格式文件，缺少第一列原子类型
 * 或者，根据输入的原子个数来创建对应的二维数组（n * 3）维，然后分别赋值（最简单）
 *
 */

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include<fstream>
#include<math.h>
#include<random>
#include<algorithm>

//输入input.txt文件，第一行为原子个数n，后面n行为每个原子的坐标
using namespace std;
void read_cnf_atoms()
{
    ifstream myfile("/Users/zlw/ClionProjects/MCCPS Luban/input.txt");
    //string temp;
    //string para[5]={"ensemble","temperature", "nmolty", "nummols", "numboxes"};
    //string ensemble;
    //double temperature;
    int nmolty = 10;
    //nummols, numboxes;



    if(!myfile.is_open())
    {
        cout<<"file open error"<<endl;
    }

    //myfile>>nmolty;


    //cout<<nmolty<<endl; for test
/*
    double **array;

    for(int i = 0; i < nmolty; i++)
    {
        array[i] = new double[4];
    }

    for(int i = 0; i < nmolty; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            myfile>>array[i][j];
            cout<<array[i][j];
        }
        cout<<endl;
    }
    */
    //int i = 0;//全局变量要小心
    vector<double> V;
    vector<double>::iterator it;

    double d;
    while(myfile >> d)
    {
        V.push_back(d);
    }

    myfile.close();

    //double aa[100][4];
    double** aa = new double*[nmolty];
    for(int i = 0; i < nmolty; i++)
    {
        aa[i] = new double[4];
    }

    for(it = V.begin(); it != V.end();)
    {
        for(int i = 0; i < nmolty; i++)
        {
            for (int j = 0; j < 4; j++)
            {

                aa[i][j] = *it;
                it++;
                cout<<' '<<aa[i][j];
            }

            cout<<endl;
        }

    }

/*
    for(string temp; getline(myfile,temp);)
    {
        //outfile<<temp<<endl;
        //由于整行整行读入，读入道str中时，文件中的每个换行符都被丢掉了，为了照原样输出，在out流上输出时需要再补上一个回车<<endl
        //int res = temp.compare(para[i++]);
        //下面三个输出用来测试


        cout<<temp<<endl;//temperature
        cout<<para[i]<<endl;//ensemble
        cout<<"res = "<<res<<endl;
        */
/*
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
*/
for(int i = 0; i < nmolty; i++)
{
    delete[] aa[i];
}
delete[] aa;


    return;
}
int main()
{
    read_cnf_atoms();

    return 0;
}