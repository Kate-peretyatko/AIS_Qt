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

    complex <float> *twiddles = new complex <float> [18];
    complex <float> one(0, 1);
    for(int t = 1; t < 18; t++)
    {
         twiddles[t] = exp(one * float((t - 1) *(M_PI/8)));
         //cout << twiddles[t] << endl;
    }

    /*complex <float> *down = new complex <float> [4];
    int y = 4;
    for(int p = 0; p < 4; p++)
    {
        down[p] = twiddles[y];
        y--;
        cout << down[p] << " ";
    }

    complex <float> *up = new complex <float> [4];
    y = 2;
    for(int p = 0; p < 4; p++)
    {
        up[p] = twiddles[y];
        y++;
        cout << up[p] << " ";
    }*/

    complex <float> synchro[] = {twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[6], twiddles[7], twiddles[8], twiddles[9], twiddles[10], twiddles[11], twiddles[12], twiddles[13],
                                 twiddles[14], twiddles[15], twiddles[16], twiddles[17], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
                                 twiddles[6], twiddles[7], twiddles[8], twiddles[9], twiddles[8], twiddles[7], twiddles[6], twiddles[5]};
    for(int u = 0; u < 128; u++)
    {
        cout << "u = " << u << " " << synchro[u]<< endl;
    }

    QCoreApplication a(argc, argv);

    return a.exec();
}

