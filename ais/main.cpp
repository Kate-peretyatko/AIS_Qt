#include <QCoreApplication>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <complex>
using namespace std;

int main(int argc, char *argv[])
{
    ifstream f;
    f.open("C:\\diplom\\ais138.dat");
    int i = 0;
    complex <float> data;
    float g;
    if(f)
    {
        int i = 0;
        while(!f.eof())
        {
            f >> data;
            cout << data.real() << " " << data.imag() << "\t";
            i++;
        }
        cout << i;
        f.close();
    }
    else
    {
        cout << "Error openning file" << endl;
    }


    QCoreApplication a(argc, argv);

    return a.exec();
}

