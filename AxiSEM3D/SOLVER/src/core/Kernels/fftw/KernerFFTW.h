//time to freq fftw for kernels 

#pragma once

#include <fftw3.h>
#include <vector>
#include "eigenc.h"

#ifdef _USE_DOUBLE
    typedef fftw_plan PlanFFTW;
    #define complexFFTW reinterpret_cast<fftw_complex*>
    #define planR2CFFTW fftw_plan_many_dft_r2c
    #define planC2RFFTW fftw_plan_many_dft_c2r
    #define destroyFFTW fftw_destroy_plan
    #define execFFTW fftw_execute
#else
    typedef fftwf_plan PlanFFTW;
    #define complexFFTW reinterpret_cast<fftwf_complex*>
    #define planR2CFFTW fftwf_plan_many_dft_r2c
    #define planC2RFFTW fftwf_plan_many_dft_c2r
    #define destroyFFTW fftwf_destroy_plan
    #define execFFTW fftwf_execute
#endif

class KernerFFTW {


};
