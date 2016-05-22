#include <QCoreApplication>
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <complex>
#include <math.h>
#include <algorithm>
#include "fasttransforms.h"


int size_signal = 32770;
int length_frame = 1024;
int length_data = 1024;

void open_read_file(float data[])
{
    std::ifstream f;
    f.open("C:\\diplom\\my.txt");

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
}

void complex_sig(float data[], std::complex <float> *com)
{
    int l = 0;
    for( int j = 0; j < size_signal; j = j + 2)
    {
        com[l] = std::complex <float> (data[j], data[j+1]);
        l++;
    }
}

void diff_complex_sig(std::complex <float> *com, std::complex <float> *diff_com)
{

    for(int k = 0; k < size_signal/2 ; k++)
    {
        diff_com[k] = com[k]*conj(com[k+1]);
        //std::cout << "k = " << k << " : " << diff_com[k] << std::endl; // СОВПАДАЕТ С MATLAB
    }
}

void twd(std::complex <float> *twiddles)
{
    std::complex <float> one(0, 1);
    for(int t = 1; t < 18; t++)
    {
         twiddles[t] = exp(one * float((t - 1) *(M_PI/8)));
         //std::cout << twiddles[t] << std::endl; // СОВПАДАЕТ С MATLAB
    }


}

void correlation(std::complex <float> synchro[], std::complex <float> *com,
                 alglib::complex *pIFFT)
{
    //ifft(fft(diff_com).*conj(fft(synchro)));
    //******************* РАЗБИВАЕМ НА ЭТАПЫ ************************

    //******************* FFT(SYNCHRO) ******************************
    alglib::complex_1d_array z;
    alglib::complex *pSynchro = new alglib::complex[length_frame];
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

    pSynchro = z.getcontent();

    /*std::ofstream fout;
    fout.open("fft_synchro.txt");
    for(int i = 0; i < length_frame; i++)
    {
        fout << "i = " << i+1 << " : " << pSynchro[i].tostring(3).c_str() << std::endl;
    }
    fout.close();*/

    //******************* CONJ(FFT(SYNCHRO)) ************************
    for(int o = 0; o < length_frame; o++)
    {
        pSynchro[o] = alglib::conj(pSynchro[o]);
       // printf("%s ", pSynchro[o].tostring(3).c_str()); // СОВПАДАЕТ С MATLAB
    }

    /*std::ofstream fout;
    fout.open("conj_fft_synchro.txt");
    for(int i = 0; i < length_frame; i++)
    {
        fout << "i = " << i+1 << " : " << pSynchro[i].tostring(3).c_str() << std::endl;
    }
    fout.close();*/

    //******** IFFT(FFT(DIFF_COM).*CONJ(FFT(SYNCHRO)))***************
    alglib::complex_1d_array y;
    alglib::complex_1d_array e;
    alglib::complex *pSignal = new alglib::complex[length_frame];
    alglib::complex *pMult = new alglib::complex[length_frame];
    alglib::complex *p = new alglib::complex[length_frame];

    std::ofstream fout;
    fout.open("ifft.txt");
    int r = 0;
    for(r = 0; r < size_signal/2-1; r = r+length_frame)
    {
        for(int u = 0; u < length_frame; u++){
        pSignal[u].x = com[u+r].real();
        pSignal[u].y = com[u+r].imag();
        }

        //************ FFT(DIFF_COM) ********************************
        y.setcontent(length_frame, pSignal);
        fftc1d(y);
        pSignal = y.getcontent();

        /*std::ofstream fout;
        fout.open("fft_diff_signal.txt");
        for(int i = 0; i < length_frame; i++)
        {
            fout << "i = " << i+1 << " : " << pSignal[i].tostring(3).c_str() << std::endl;
        }
        fout.close();*/
        //printf("%s ", y.tostring(3).c_str()); // СОВПАДАЕТ С MATLAB
        //********** FFT(DIFF_COM).*CONJ(FFT(SYNCHRO)) **************
        for( int h = 0; h <  length_frame; h++)
        {
            //pMult[h] = pSynchro[h]*pSignal[h];
            pMult[h].x = pSynchro[h].x*pSignal[h].x - pSynchro[h].y*pSignal[h].y;
            pMult[h].y = pSynchro[h].y*pSignal[h].x + pSynchro[h].x*pSignal[h].y;
            //printf("%i %s \n",h, pMult[h].tostring(3).c_str()); // СОВПАЛО С MATLAB ДО 128 ОТСЧЁТА
        }
        /*std::ofstream fout;
        fout.open("Mult.txt");
        for(int i = 0; i < length_frame; i++)
        {
             fout << "i = " << i+1 << " : " << pMult[i].tostring(3).c_str() << std::endl;
        }
        fout.close();*/
        //*********** IFFT() ****************************************
        e.setcontent(length_frame, pMult);
        fftc1dinv(e);
        p = e.getcontent();

        //*********** ЗАПИСЫВАЕМ В БОЛЬШОЙ ВЕКТОР *******************
        int k = 0;
        for(int i = r; i < r+length_frame; i++)
        {
            pIFFT[i] = p[k];
            fout << "i = " << i+1 << " : " << pIFFT[i].tostring(3).c_str() << std::endl;
            k++;
        }
    }
    fout.close();
}

void find_max(alglib::complex *pIFFT, alglib::complex *max, int index[])
{
    int u = 0;
    for(int r = 0; r < size_signal/2-1; r = r+length_frame)
    {
        index[u] = 0;
        for(int i = r; i < r+length_frame; i++)
        {
            if(alglib::abscomplex(pIFFT[i]) > alglib::abscomplex(max[u]))
            {
                max[u] = pIFFT[i];
                index[u] = i;
            }
        }
        u++;
    }
}

void norm_max(alglib::complex *max, int index[], alglib::complex *norm_max_corr)
{
    alglib::complex Max_corr = (0, 0);
    int id = 0;
    for(int i = 0; i < 16; i++)
    {
        std::cout << "max : " << max[i].tostring(3).c_str() << "   index = " << index[i]+1 << std::endl;
        if(alglib::abscomplex(max[i]) > alglib::abscomplex(Max_corr))
        {
           Max_corr = max[i];
        }
    }
    std::cout << "Max_corr : " << Max_corr.tostring(3).c_str() << std::endl;

    for(int i = 0; i < 16; i++)
    {
        norm_max_corr[i] = max[i]/Max_corr;
        std::cout << "norm_max_corr : " << norm_max_corr[i].tostring(3).c_str() << std::endl;
    }
}

void pars_package(alglib::complex *norm_max_corr, int index[], std::complex <float> *com, std::complex <float> *package)
{
    int number_index = 0;
    alglib::complex m(0,0);
    for(int i = 0; i < 16; i++)
    {
        if(alglib::abscomplex( norm_max_corr[i]) > alglib::abscomplex(m))
        {
           m =  norm_max_corr[i];
           number_index = index[i];
        }
    }
    std::cout << "m : " << m.tostring(3).c_str() << " " << number_index+1 << std::endl;

    for(int i = 0; i < length_data; i++)
    {
        package[i] = com[number_index+i];
        //std::cout << package[i] << std::endl;
    }
}

void demodulation(std::complex <float> *package)
{
    std::complex <float> one(0, 1);
    std::complex <float> *sig1 = new std::complex <float> [5];
    std::complex <float> *sig2 = new std::complex <float> [5];
    sig1[0] = 1;
    sig1[1] = exp(one*float(M_PI/8));
    sig1[2] = exp(one*float(M_PI/4));
    sig1[3] = exp(one*float(3*M_PI/8));
    sig1[4] = one;
    for(int i = 0; i < 5; i++)
    {
        sig2[i] = std::conj(sig1[i]);
        //std::cout << sig2[i] << std::endl;
    }
    std::complex <float> *x = new std::complex <float> [5];
    int bits[length_data/4];
    int ii = 0;
    for(int i = 2; i < length_data; i = i + 4)
    {
        std::complex <float> pr1(0,0);
        std::complex <float> pr2(0,0);

        x[0] = package[i];
        x[1] = package[i+1];
        x[2] = package[i+2];
        x[3] = package[i+3];
        x[4] = package[i+4];
        for(int step = 0; step < 5; step++)
        {
            pr1 = pr1 + x[step]*sig1[step];
            pr2 = pr2 + x[step]*sig2[step];
        }

        //std::cout << "pr1 : " << pr1 << " pr2 : " << pr2 << std::endl;
        if(pr1.real() > pr2.real())
        {
            bits[ii] = 0;
        }
        else
        {
            bits[ii] = 1;
        }
        ii++;
    }
    for( int i = 0; i < 32; i++)
    {
        if(bits[i] == bits[i+1])
            bits[i] = 1;
        else
            bits[i] = 0;
        std::cout << bits[i] << " ";
    }
}

int main(int argc, char *argv[])
{
    //*********** ОТКРЫТИЕ ФАЙЛА ************************************    
    float data[size_signal];
    open_read_file(data);

    //**************** СКЛЕИВАНИЕ В КОМПЛЕКСНЫЕ ЧИСЛА ***************
    std::complex <float> *com = new std::complex  <float> [size_signal/2];
    complex_sig(data, com);

    /*std::ofstream fout;
    fout.open("signal.txt");
    for(int i = 0; i < size_signal/2; i++)
    {
        fout << "i = " << i+1 << " : " << com[i] << std::endl;
    }
    fout.close();*/
    //***************** ПРИВЕДЕНИЕ К ДИФФ. ВИДУ *********************
    std::complex <float> *diff_com = new std::complex  <float> [size_signal/2];
    diff_complex_sig(com, diff_com);

    /*std::ofstream fout;
    fout.open("diff_signal.txt");
    for(int i = 0; i < size_signal/2; i++)
    {
        fout << "i = " << i+1 << " : " << diff_com[i] << std::endl;
    }
    fout.close();*/
    //****************** ПОДГОТОВКА ПРЕАМБУЛЫ ***********************
    std::complex <float> *twiddles = new std::complex <float> [18];
    twd(twiddles);
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


    //****************** БЫСТРАЯ КОРЕЕЛЯЦИЯ *************************
    alglib::complex *pIFFT = new alglib::complex[size_signal/2];
    correlation(synchro, com, pIFFT);

    //*********************** ПОИСК МАКСИМУМА В КАЖДОМ ФРЕЙМЕ *******
    alglib::complex *max = new alglib::complex[16];
    int index[32];
    find_max(pIFFT, max, index);

    //********************** НОРМИРОВАНИЕ ***************************
    //******************** ПОИСК САМОГО МАКСИМАЛЬНОГО ЭЛЕМЕНТА ******
    alglib::complex *norm_max_corr = new alglib::complex[16];
    norm_max(max, index, norm_max_corr);

    //************ ВЫДЕЛЕНИЕ ПАКЕТА *********************************    
    std::complex <float> *package = new std::complex  <float> [length_data];
    pars_package(norm_max_corr, index, com, package);

    //********************** ДЕМОДУЛЯТОР ****************************
    demodulation(package);

    QCoreApplication a(argc, argv);
    return a.exec();
}

