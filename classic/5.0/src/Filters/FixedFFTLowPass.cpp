#include "Filters/FixedFFTLowPass.h"
#include "fftw3.h"
#include "Physics/Physics.h"

FixedFFTLowPassFilter::FixedFFTLowPassFilter(const int &N): 
  number_frequencies_m(N) 
{ }

void FixedFFTLowPassFilter::apply(vector<double> &LineDensity)
{
  const int M = LineDensity.size();

  fftw_complex *FourCoefs = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* M);
  fftw_complex *PsudoComplexValues = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* M);
  fftw_plan p;

  for (int i = 0; i < M; ++ i) {
    PsudoComplexValues[i][0] = LineDensity[i];
    PsudoComplexValues[i][1] = 0.0;
  }

  p = fftw_plan_dft_1d(M, PsudoComplexValues, FourCoefs, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  FourCoefs[0][0] /= M;
  for (int i = 1; i < number_frequencies_m + 1; ++ i) {
    FourCoefs[i][0] /= M;
    FourCoefs[i][1] /= M;
    FourCoefs[M - i][0] /= M;
    FourCoefs[M - i][1] /= M;
  }

  for (int i = number_frequencies_m + 1; i < M - number_frequencies_m; ++ i) {
    FourCoefs[i][0] = 0.0;
    FourCoefs[i][1] = 0.0;
  }

  p = fftw_plan_dft_1d(M, FourCoefs, PsudoComplexValues, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  for (int i = 0; i < M; ++ i) {
    LineDensity[i] = PsudoComplexValues[i][0];
  }

  fftw_destroy_plan(p);
  fftw_free(FourCoefs);
  fftw_free(PsudoComplexValues);
}

void FixedFFTLowPassFilter::calc_derivative(vector<double> &LineDensity, const double &h)
{
  const int M = LineDensity.size();
  const double gff = 2.* Physics::pi / ((M-1) * h * M);
  double temp;

  fftw_complex *FourCoefs = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* M);
  fftw_complex *PsudoComplexValues = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* M);
  fftw_plan p;

  for (int i = 0; i < M; ++ i) {
    PsudoComplexValues[i][0] = LineDensity[i];
    PsudoComplexValues[i][1] = 0.0;
  }

  p = fftw_plan_dft_1d(M, PsudoComplexValues, FourCoefs, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  FourCoefs[0][0] = 0.0;
  FourCoefs[0][1] = 0.0;
  for (int i = 1; i < number_frequencies_m + 1; ++ i) {
    temp = FourCoefs[i][0];
    FourCoefs[i][0] = -FourCoefs[i][1] * gff * i;
    FourCoefs[i][1] = temp * gff * i;
    temp = FourCoefs[M - i][0];
    FourCoefs[M - i][0] = FourCoefs[M - i][1] * gff * i;
    FourCoefs[M - i][1] = -temp * gff * i;
  }

  for (int i = number_frequencies_m + 1; i < M - number_frequencies_m; ++ i) {
    FourCoefs[i][0] = 0.0;
    FourCoefs[i][1] = 0.0;
  }

  p = fftw_plan_dft_1d(M, FourCoefs, PsudoComplexValues, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  for (int i = 0; i < M; ++ i) {
    LineDensity[i] = PsudoComplexValues[i][0];
  }

  fftw_destroy_plan(p);
  fftw_free(FourCoefs);
  fftw_free(PsudoComplexValues);

}

