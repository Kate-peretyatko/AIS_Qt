#include "stdafx.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fasttransforms.h"


using namespace alglib;
using namespace std;

int fft(int argc, char **argv)
{
    //
    // first we demonstrate forward FFT:
    // [1i,1i,1i,1i] is converted to [4i, 0, 0, 0]

    complex_1d_array z = "[ 0.3827+0.9239i,  0.7071+0.7071i,  0.9239+0.3827i,   1.0000 ]";
    complex_1d_array y = "[ 1.5330e+005 +1.3392e+006i, 2.2845e+006 -2.1809e+006i, -1.9066e+006 -1.3138e+006i, -1.1442e+005 +1.4571e+006i ]";
    fftc1d(z);
    fftc1d(y);
    printf("%s\n", z.tostring(3).c_str()); // EXPECTED: [4i,0,0,0]
    printf("%s\n",y.tostring(3).c_str());
    //
    // now we convert [4i, 0, 0, 0] back to [1i,1i,1i,1i]
    // with backward FFT
    //
    fftc1dinv(z);
    printf("%s\n", z.tostring(3).c_str()); // EXPECTED: [1i,1i,1i,1i]
    return 0;
}

