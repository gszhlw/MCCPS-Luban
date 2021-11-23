#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

int main() {

    ifstream myfile("/Users/zlw/ClionProjects/MCCPS Luban/hello.txt");
    //ofstream  outfile("/Users/zlw/ClionProjects/MCCPS Luban/out.txt");
    string temp;
    string para[5]={"ensemble","temperature", "nmolty", "nummols", "numboxes"};

    if(!myfile.is_open())
    {
        cout<<"file open error"<<endl;
    }

//按行读，读进来判断，有点傻，但是可行
//getline函数的作用为从输入文件流中读入一行数据，放入string变量temp中
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
        getline(myfile,temp);
        cout<<"真正输出的参数值："<<temp<<endl;
    }
    else
    {
        cout<<"parameter error"<<endl;
    }

}


    myfile.close();
    //outfile.close();
    return 0;
}
