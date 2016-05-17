#include <QCoreApplication>
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <complex>
#include <math.h>
#include "fasttransforms.h"


int size_signal = 32770;
int length_frame = 1024;

int main(int argc, char *argv[])
{
    //*********** ОТКРЫТИЕ ФАЙЛА ************************************
    std::ifstream f;
    f.open("C:\\diplom\\my.txt");

    float data[size_signal];
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
    //**************** СКЛЕИВАНИЕ В КОМПЛЕКСНЫЕ ЧИСЛА ***************
    std::complex <float> *com = new std::complex  <float> [size_signal];
    int l = 0;
    for( int j = 0; j < 32770; j = j + 2)
    {
        com[l] = std::complex <float> (data[j], data[j+1]);
        l++;
    }
    //***************** ПРИВЕДЕНИЕ К ДИФФ. ВИДУ *********************
    std::complex <float> *diff_com = new std::complex  <float> [size_signal];
    for(int k = 0; k < l ; k++)
    {
        diff_com[k] = com[k]*conj(com[k+1]);
        //std::cout << diff_com[k] << std::endl; // СОВПАДАЕТ С MATLAB
    }
    //****************** ПОДГОТОВКА ПРЕАМБУЛЫ ***********************
    std::complex <float> *twiddles = new std::complex <float> [18];
    std::complex <float> one(0, 1);
    for(int t = 1; t < 18; t++)
    {
         twiddles[t] = exp(one * float((t - 1) *(M_PI/8)));
         //std::cout << twiddles[t] << std::endl; // СОВПАДАЕТ С MATLAB
    }

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
        //std::cout << "u = " << u << " " << synchro[u]<< std::endl; // СОВПАДАЕТ С MATLAB
    }

    //****************** БЫСТРАЯ КОРЕЕЛЯЦИЯ *************************
    //ifft(fft(diff_com).*conj(fft(synchro)));
    //******************* РАЗБИВАЕМ НА ЭТАПЫ ************************

    //******************* FFT(SYNCHRO) ******************************
    alglib::complex_1d_array z;
    alglib::complex *pSynchro = new alglib::complex[1024];
    for(int u = 0; u < 128; u++){
        pSynchro[u].x = synchro[u].real();
        pSynchro[u].y = synchro[u].imag();
    }
    //*******************ЗАПОЛНЕНИЕ ПРЕАМБУЛЫ НУЛЯМИ ****************
    for(int u = 128; u < length_frame; u++){
        pSynchro[u].x = 0;
        pSynchro[u].y = 0;
    }
    z.setcontent(1024, pSynchro);
    fftc1d(z);
    //printf("%s ",z.tostring(3).c_str()); // СОВПАДАЕТ С MATLAB

    //******************* CONJ(FFT(SYNCHRO)) ************************
    pSynchro = z.getcontent();
    for(int o = 0; o < 128; o++)
    {
        pSynchro[o] = alglib::conj(pSynchro[o]);
        //printf("%s ", pContent[o].tostring(3).c_str()); // СОВПАДАЕТ С MATLAB
    }

    //******** IFFT(FFT(DIFF_COM).*CONJ(FFT(SYNCHRO)))***************
    alglib::complex_1d_array y;
    alglib::complex_1d_array e;
    alglib::complex *pSignal = new alglib::complex[size_signal];
    alglib::complex *pMult = new alglib::complex[size_signal];
    alglib::complex *pIFFT = new alglib::complex[size_signal];
    int r = 0;
    //for(int r = 0; r < size_signal; r+length_frame)
    {
        for(int u = r; u < r+length_frame; u++){
        pSignal[u].x = diff_com[u].real();
        pSignal[u].y = diff_com[u].imag();
        }
        //************ FFT(DIFF_COM) ********************************
        y.setcontent(length_frame, pSignal);
        fftc1d(y);
        pSignal = y.getcontent();
        //printf("%s ", y.tostring(3).c_str()); // СОВПАДАЕТ С MATLAB
        //********** FFT(DIFF_COM).*CONJ(FFT(SYNCHRO)) **************
        for( int h = r; h <  r+length_frame; h++)
        {
            //pMult[h] = pSynchro[h]*pSignal[h];
            pMult[h].x = pSynchro[h].x*pSignal[h].x - pSynchro[h].y*pSignal[h].y;
            pMult[h].y = pSynchro[h].y*pSignal[h].x + pSynchro[h].x*pSignal[h].y;
            //printf("%i %s \n",h, pMult[h].tostring(3).c_str()); // СОВПАЛО С MATLAB ДО 128 ОТСЧЁТА
        }
        //*********** IFFT() ****************************************
        e.setcontent(length_frame, pMult);
        fftc1dinv(e);
        pIFFT = e.getcontent();
        for(int i = 0; i < 150; i++)
        {
           // printf("%i %s \n", i, pIFFT[i].tostring(3).c_str());
        }
         alglib::complex_1d_array t = "[0,0,0,0]";
         fftc1d(t);
         printf("%s ", t.tostring(3).c_str());
    }


    QCoreApplication a(argc, argv);
    return a.exec();
}

