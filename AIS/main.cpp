#include <QCoreApplication>
#include "ais.h"

#include <fstream>
#include <iostream>

int size_signal = 32770;
int length_frame = 1024;
int length_data = 1024;

int main(int argc, char *argv[])
{
    char file_name[] = "ais263.dat";
    module_ais(size_signal, length_frame, length_data, file_name);

    //****************** ТЕСТИРОВАНИЕ *******************************
    std::ifstream i1, i2;
    // ФАЙЛ С РЕЗУЛЬТАТОМ MATLAB
    i1.open("C:\\diplom\\test");
    i2.open("test.txt");
    std::string s1, s2;
    while (i1)
    {
       getline(i1,s1);
    }
    while (i2)
    {
       getline(i2,s2);
    }
    int k = s1.compare(s2);
    if(k == 0)
    {
        std::cout << "The result coincided with Matlab!";
    }
    else
        std::cout << "The result IS NOT coincided with Matlab!";
    i1.close();
    i2.close();

    QCoreApplication a(argc, argv);
    return a.exec();
}
