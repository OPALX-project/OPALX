#ifndef WAKE_H_
#define WAKE_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <string>
#include <fftw3.h>
#include "simpson.h"


using namespace std;

void calcWakeSerie();
void CalcWakeFFT(ofstream* WakeStream, ofstream* FFTStream, int Lbunch, string material, string direction, string Mode, double a_Wake);
void fft(double* in , fftw_complex  *out, int N);
void CalcWake(int Lbunch, string material, string direction, string mode, double a_Wake, double *wake);



#endif /*WAKE_H_*/
