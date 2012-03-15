#include "Filters/RelativeFFTLowPass.h"
#include "fftw3.h"
#include "Physics/Physics.h"

RelativeFFTLowPassFilter::RelativeFFTLowPassFilter(const double &threshold): 
  threshold_m(threshold) 
{ }

void RelativeFFTLowPassFilter::apply(vector<double> &LineDensity)
{
  const int M = LineDensity.size();
  double max_four_coef;

  fftw_complex *FourCoefs = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* M);
  fftw_complex *PsudoComplexValues = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* M);
  fftw_plan p;

  for (int i = 0; i < M; ++ i) {
    PsudoComplexValues[i][0] = LineDensity[i];
    PsudoComplexValues[i][1] = 0.0;
  }

  p = fftw_plan_dft_1d(M, PsudoComplexValues, FourCoefs, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  for (int i = 0; i < M; ++ i) {
    FourCoefs[i][0] /= M;
    FourCoefs[i][1] /= M;
    if (fabs(FourCoefs[i][0]) > max_four_coef) {
      max_four_coef = fabs(FourCoefs[i][0]);
    }
    if (fabs(FourCoefs[i][1]) > max_four_coef) {
      max_four_coef = fabs(FourCoefs[i][1]);
    }
  }

  for (int i = 0; i < M; ++ i) {
    if (fabs(FourCoefs[i][0])/max_four_coef < threshold_m) {
      FourCoefs[i][0] = 0.0;
    }
    if (fabs(FourCoefs[i][1])/max_four_coef < threshold_m) {
      FourCoefs[i][1] = 0.0;
    }
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

void RelativeFFTLowPassFilter::calc_derivative(vector<double> &LineDensity, const double &h)
{
  const int M = LineDensity.size();
  const double gff = 2.* Physics::pi / ((M-1) * h * M);
  double max_four_coef;
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

  for (int i = 1; i < M; ++ i) {
    FourCoefs[i][0] /= M;
    FourCoefs[i][1] /= M;
    if (fabs(FourCoefs[i][0]) > max_four_coef) {
      max_four_coef = fabs(FourCoefs[i][0]);
    }
    if (fabs(FourCoefs[i][1]) > max_four_coef) {
      max_four_coef = fabs(FourCoefs[i][1]);
    }
  }

  FourCoefs[0][0] = 0.0;
  FourCoefs[0][1] = 0.0;
  for (int i = 1; i < M; ++ i) {
    temp = FourCoefs[i][0];
    if (fabs(FourCoefs[i][1])/max_four_coef < threshold_m) {
      FourCoefs[i][0] = 0.0;
    } else {
      FourCoefs[i][0] = -FourCoefs[i][1] * gff * i;
    }
    if (fabs(FourCoefs[i][0])/max_four_coef < threshold_m) {
      FourCoefs[i][1] = 0.0;
    } else {
      FourCoefs[i][1] = temp * gff * i;
    }
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
