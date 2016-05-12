#include <QCoreApplication>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <complex>
#include <math.h>
using namespace std;

int main(int argc, char *argv[])
{
    ifstream f;
    f.open("C:\\diplom\\my.txt");

    float data[32770];
    if(f)
    {
        int i = 0;
        while(!f.eof())
        {
            f >> data[i];
            i++;
        }
        f.close();
    }
    else
    {
        cout << "Error openning file" << endl;
    }
    complex <float> *com = new complex  <float> [32770];
    int l = 0;
    for( int j = 0; j < 32770; j = j + 2)
    {
        com[l] = complex <float> (data[j], data[j+1]);
        l++;
    }

    complex <float> *diff_com = new complex  <float> [32770];
    for(int k = 0; k < l ; k++)
    {
        diff_com[k] = com[k]*conj(com[k+1]);
        //cout << diff_com[k] << endl;
    }

    complex <double> *twiddles = new complex <double> [18];
    complex <float> one(0, 1);
    for(int t = 1; t < 18; t++)
    {
        twiddles[t] = complex <double> (exp(one.imag()*(t-1)*M_PI/8));
        cout << twiddles[t] << endl;
    }
    //double synchro[] = {twiddles[4], twiddles[3], twiddles[2], twiddles[1]};
    //cout << synchro[0] << " " << synchro[1] << " " << synchro[2] << " " << synchro[3] << " " << endl;

    QCoreApplication a(argc, argv);

    return a.exec();
}

