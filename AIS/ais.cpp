#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <complex>
#include <math.h>
#include <cmath>
#include <algorithm>
#include "fasttransforms.h"

#include <string>
#include "ais.h"


// ФУНКЦИЯ ОСУЩЕСТВЛЯЕТ ОТКРЫТИЕ ФАЙЛА И ЧТЕНИЕ ИЗ НЕГО
void open_read_file(float data[], char file_name[])
{
    std::ifstream f;
    f.open(file_name);

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
//*******************************************************************
// ПОЛУЧЕНИЕ КОМПЛЕСНОГО СИГНАЛА
void complex_sig(float data[], std::complex <float> *com, int size_signal)
{
    int l = 0;
    for( int j = 0; j < size_signal; j = j + 2)
    {
        com[l] = std::complex <float> (data[j], data[j+1]);
        l++;
    }
}
//*******************************************************************
// ПРИВЕДЕНИЕ СИГНАЛА К ДИФФЕРЕНЦИАЛЬНОМУ ВИДУ
void diff_complex_sig(std::complex <float> *com, std::complex <float> *diff_com, int size_signal)
{

    for(int k = 0; k < size_signal/2 ; k++)
    {
        diff_com[k] = com[k]*conj(com[k+1]);
    }
}
//*******************************************************************
// ДЕЛЕНИЕ КОМПЛЕКСНОЙ ПЛОСКОСТИ
void twd(std::complex <float> *twiddles)
{
    std::complex <float> one(0, 1);
    for(int t = 1; t < 18; t++)
    {
         twiddles[t] = exp(one * float((t - 1) *(M_PI/8)));
    }
}
//*******************************************************************
// КОРРЕЛЯЦИОННОЕ ОБНАРУЖЕНИЕ ПАКЕТА
void correlation(std::complex <float> synchro[], std::complex <float> *com,
                 alglib::complex *pIFFT, int length_frame, int size_signal)
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

    //******************* CONJ(FFT(SYNCHRO)) ************************
    for(int o = 0; o < length_frame; o++)
    {
        pSynchro[o] = alglib::conj(pSynchro[o]);
    }

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

        //********** FFT(DIFF_COM).*CONJ(FFT(SYNCHRO)) **************
        for( int h = 0; h <  length_frame; h++)
        {
            pMult[h].x = pSynchro[h].x*pSignal[h].x - pSynchro[h].y*pSignal[h].y;
            pMult[h].y = pSynchro[h].y*pSignal[h].x + pSynchro[h].x*pSignal[h].y;
        }

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
//*******************************************************************
// ПОИСК МАКСИМАЛЬНОГО ЗНАЧЕНИЯ КОРРЕЛЯЦИИ В КАЖДОМ ФРЕЙМЕ
void find_max(alglib::complex *pIFFT, alglib::complex *max, int index[], int size_signal, int length_frame)
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
//*******************************************************************
// ПОИСК МАКСИМАЛЬНОГО ЗНАЧЕНИЯ КОРРЕЛЯЦИИ СРЕДИ ФРЕЙМОВ, НОРМИРОВАНИЕ ЗНАЧЕНИЙ КОРРЕЛЯЦИИ
void norm_max(alglib::complex *max, int index[], alglib::complex *norm_max_corr)
{
    alglib::complex Max_corr = (0, 0);
    int id = 0;
    for(int i = 0; i < 16; i++)
    {
        if(alglib::abscomplex(max[i]) > alglib::abscomplex(Max_corr))
        {
           Max_corr = max[i];
        }
    }
    for(int i = 0; i < 16; i++)
    {
        norm_max_corr[i] = max[i]/Max_corr;
    }
}
//*******************************************************************
// ИЗВЛЕЧЕНИЕ ПАКЕТА ИЗ СИГНАЛА
void pars_package(alglib::complex *norm_max_corr, int index[], std::complex <float> *com,
                  std::complex <float> *package, int length_data)
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
    for(int i = 0; i < length_data; i++)
    {
        package[i] = com[number_index+i];
    }
}
//*******************************************************************
// ДЕМОДУЛЯЦИЯ ПАКЕТА
void demodulation(std::complex <float> *package, int bits[], int length_data)
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
    }
    std::complex <float> *x = new std::complex <float> [5];
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
    std::cout << "preamble ";
    for( int i = 0; i < 32; i++)
    {
        if(bits[i] == bits[i+1])
            bits[i] = 1;
        else
            bits[i] = 0;
        std::cout << bits[i] << " ";
    }
    std::cout << std::endl;
}
//*******************************************************************
// ПЕРЕВОД СЕГМЕНТА ДАННЫХ ПАКЕТА АИС ИЗ 2ой СИСТЕМЫ СЧИСЛЕНИЯ В 10ую
void translation(int coordinate_data[])
{
    int id_data[32];
    int latitude[32];
    int longitude[32];
    for(int i = 0; i < 32; i++)
    {
        id_data[i] = coordinate_data[i];
        latitude[i] = coordinate_data[32+i];
        longitude[i] = coordinate_data[64+i];
    }

    int ID = 0;
    int lat = 0;
    int lon = 0;
    for(int i = 31; i >= 0; i--)
    {
        ID = ID + pow(2,id_data[i]);
        lat = lat + pow(2,latitude[i]);
        lon = lon + pow(2,longitude[i]);
    }
    std::cout << "ID = " << ID << ", latitude = " << lat << ", longitude = " << lon;
    std::ofstream of("test.txt");
    if (!of) {
        std::cout << "Cannot open file.\n";
    }
    of << ID << " " << lat << " " << lon << " ";

    int speed_data[16];
    int course_data[16];
    for(int i = 0; i < 16; i++)
    {
        speed_data[i] = coordinate_data[96+i];
        course_data[i] = coordinate_data[112+i];
    }

    int speed = 0;
    int course = 0;
    for(int i = 15; i >= 0; i--)
    {
        speed = speed + pow(2,speed_data[i]);
        course = course + pow(2, course_data[i]);
    }
    std::cout << " speed = " << speed << ", course = " << course << std::endl;
    of << speed << " " << course;
    of.close();
}
//*******************************************************************
//****************** МОДУЛЬ АИС *************************************
// ВЫПОЛНЯЕТ ПОСЛЕДОВАТЕЛЬНУЮ ПОЛНУЮ ОБРАБОТКУ ПАКЕТА
void module_ais(int size_signal, int length_frame, int length_data, char file_name[])
{
    //*********** ОТКРЫТИЕ ФАЙЛА ************************************
    float data[size_signal];
    open_read_file(data, file_name);

    //**************** СКЛЕИВАНИЕ В КОМПЛЕКСНЫЕ ЧИСЛА ***************
    std::complex <float> *com = new std::complex  <float> [size_signal/2];
    complex_sig(data, com, size_signal);

    //***************** ПРИВЕДЕНИЕ К ДИФФ. ВИДУ *********************
    std::complex <float> *diff_com = new std::complex  <float> [size_signal/2];
    diff_complex_sig(com, diff_com, size_signal);

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
    correlation(synchro, com, pIFFT, length_frame, size_signal);

    //*********************** ПОИСК МАКСИМУМА В КАЖДОМ ФРЕЙМЕ *******
    alglib::complex *max = new alglib::complex[16];
    int index[32];
    find_max(pIFFT, max, index, size_signal, length_frame);

    //********************** НОРМИРОВАНИЕ ***************************
    //******************** ПОИСК САМОГО МАКСИМАЛЬНОГО ЭЛЕМЕНТА ******
    alglib::complex *norm_max_corr = new alglib::complex[16];
    norm_max(max, index, norm_max_corr);

    //************ ВЫДЕЛЕНИЕ ПАКЕТА *********************************
    std::complex <float> *package = new std::complex  <float> [length_data];
    pars_package(norm_max_corr, index, com, package, length_frame);

    //********************** ДЕМОДУЛЯТОР ****************************
    int bits[length_data/4];
    demodulation(package, bits, length_data);

    //******************** ПАРСИНГ **********************************
    int coordinate_data[168];
    for(int i = 0; i < 168; i++)
    {
        coordinate_data[i] = bits[32+i];
    }

    //****************** ПЕРЕВОД В 10ую СИСТЕМУ *********************
    translation(coordinate_data);

}
