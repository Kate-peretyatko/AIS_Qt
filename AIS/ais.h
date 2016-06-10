#ifndef AIS_H
#define AIS_H

#include <complex>
#include "fasttransforms.h"

//*********** ОТКРЫТИЕ ФАЙЛА ************************************
void open_read_file(float data[], char file_name[]);
//**************** СКЛЕИВАНИЕ В КОМПЛЕКСНЫЕ ЧИСЛА ***************
void complex_sig(float data[], std::complex <float> *com, int size_signal);
//***************** ПРИВЕДЕНИЕ К ДИФФ. ВИДУ *********************
void diff_complex_sig(std::complex <float> *com, std::complex <float> *diff_com,
                      int size_signal);
//***************** ДЕЛЕНИЕ КОМПЛЕКСНОЙ ПЛОСКОСТИ ***************
void twd(std::complex <float> *twiddles);
//****************** БЫСТРАЯ КОРЕЕЛЯЦИЯ *************************
void correlation(std::complex <float> synchro[], std::complex <float> *com,
                 alglib::complex *pIFFT, int length_frame, int size_signal);
//*********************** ПОИСК МАКСИМУМА В КАЖДОМ ФРЕЙМЕ *******
void find_max(alglib::complex *pIFFT, alglib::complex *max, int index[],
              int size_signal, int length_frame);
//********************** НОРМИРОВАНИЕ ***************************
//******************** ПОИСК САМОГО МАКСИМАЛЬНОГО ЭЛЕМЕНТА ******
void norm_max(alglib::complex *max, int index[], alglib::complex *norm_max_corr);
//************ ВЫДЕЛЕНИЕ ПАКЕТА *********************************
void pars_package(alglib::complex *norm_max_corr, int index[], std::complex <float> *com,
                  std::complex <float> *package, int length_data);
//********************** ДЕМОДУЛЯТОР ****************************
void demodulation(std::complex <float> *package, int bits[], int length_data);
//****************** ПЕРЕВОД В 10ую СИСТЕМУ *********************
void translation(int coordinate_data[]);

//****************** МОДУЛЬ АИС *********************************
void module_ais(int size_singnal, int length_frame, int length_data, char file_name[]);

#endif // AIS_H
