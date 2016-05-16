#include <QCoreApplication>
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <complex>
#include <math.h>
#include "fasttransforms.h"
//using namespace std;
//using namespace alglib;

int main(int argc, char *argv[])
{
    std::ifstream f;
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
        std::cout << "Error openning file" << std::endl;
    }

    std::complex <float> *com = new std::complex  <float> [32770];
    int l = 0;
    for( int j = 0; j < 32770; j = j + 2)
    {
        com[l] = std::complex <float> (data[j], data[j+1]);
        l++;
    }

    std::complex <float> *diff_com = new std::complex  <float> [32770];
    for(int k = 0; k < l ; k++)
    {
        diff_com[k] = com[k]*conj(com[k+1]);
        //cout << diff_com[k] << endl;
    }

    std::complex <float> *twiddles = new std::complex <float> [18];
    std::complex <float> one(0, 1);
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

    std::complex <float> synchro[] = {twiddles[4], twiddles[3], twiddles[2], twiddles[1], twiddles[2], twiddles[3], twiddles[4], twiddles[5],
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
        //cout << "u = " << u << " " << synchro[u]<< endl;
    }
    alglib::complex_1d_array z;

        alglib::complex *pContent = new alglib::complex[128];
        for(int u = 0; u < 128; u++){
            pContent[u].x = synchro[u].real();
            pContent[u].y = synchro[u].imag();
        }
        z.setcontent(128, pContent);
        fftc1d(z);
        pContent = z.getcontent();
        for(int o = 0; o < 128; o++)
        {
            pContent[o] = alglib::conj(pContent[o]);
            //printf("%s", pContent[o].tostring(3).c_str());
        }


            alglib::complex_1d_array y;

            alglib::complex *pSignal = new alglib::complex[32770];
            for(int u = 0; u < 32770; u++){
                pSignal[u].x = diff_com[u].real();
                pSignal[u].y = diff_com[u].imag();
            }
            //printf("%s ", pSignal[0].tostring(3).c_str());
            y.setcontent(32770, pSignal);
            fftc1d(y);

            pSignal = y.getcontent();
            for(int o = 0; o < 30; o++)
            {
                //printf("%s ", pSignal[o].tostring(3).c_str());
            }

            alglib::complex *pSynchro = new alglib::complex[32770];
            for(int u = 0; u < 128; u++){
                pSynchro[u].x = pContent[u].x;
                pSynchro[u].y = pContent[u].y;
            }
            for(int u = 128; u < 32770; u++){
                pSynchro[u].x = 0;
                pSynchro[u].y = 0;
            }

            alglib::complex_1d_array e;

            alglib::complex *pMult = new alglib::complex[32770];
            for( int h = 0; h < 32770; h++)
            {
                pMult[h].x = pSynchro[h].x*pSignal[h].x - pSynchro[h].y*pSignal[h].y;
                pMult[h].y = pSynchro[h].x*pSignal[h].y + pSynchro[h].y*pSignal[h].x;
            }
            e.setcontent(32770, pMult);
            fftc1dinv(e);
            printf("%s ", e.tostring(3).c_str());

    std::cout << "USPECH" << std::endl;
    //std::cout << diff_com[0].real() << " " << diff_com[0].imag() << std::endl;
    QCoreApplication a(argc, argv);
    //printf("%s", y.tostring(3).c_str());
    return a.exec();
}

