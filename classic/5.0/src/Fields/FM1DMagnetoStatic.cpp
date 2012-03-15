#include <fstream>
#include <ios>

#include "Fields/FM1DMagnetoStatic.hh"
#include "Physics/Physics.h"
#include "fftw3.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

#define LONGITUDINALFIRST

FM1DMagnetoStatic::FM1DMagnetoStatic(string aFilename)
  :Filename_m(aFilename),
   realFourCoefs_m(NULL),
   imagFourCoefs_m(NULL)
{
  int tmpInt;
  string tmpString;
  double tmpDouble;

  Type = T1DMagnetoStatic;
  ifstream file(Filename_m.c_str());

  if (file.good())
    {
      file >> tmpString >> accuracy_m;
      file >> zbegin_m >> zend_m >> num_gridpz_m;
      file >> rbegin_m >> rend_m >> tmpInt;

      rbegin_m /= 100.;
      rend_m /= 100.;
      zbegin_m /= 100.;
      zend_m /= 100.;

      num_gridpz_m++;

      length_m = num_gridpz_m * (zend_m - zbegin_m) / (num_gridpz_m - 1);

      file.close();
    }
  else
    {
      Inform msg("FM1DMS ");
      msg << "* ************** W A R N I N G *****************************************************" << endl;
      msg << "* NO FILE \"" << Filename_m << "\" found!" << endl;
      msg << "* **********************************************************************************" << endl;
    }
}

FM1DMagnetoStatic::~FM1DMagnetoStatic()
{
  if (realFourCoefs_m != NULL)
    {
      delete[] realFourCoefs_m;
      delete[] imagFourCoefs_m;
    }
}

void FM1DMagnetoStatic::readMap()
{
  if (realFourCoefs_m == NULL)
    {
      Inform msg("FM1DMS ");
      ifstream in(Filename_m.c_str());

      int tmpInt;
      string tmpString;
      double tmpDouble;

      int num_gridpzp;
      double *RealValues;
      fftw_complex* FourCoefs;
      fftw_plan p;
      double ez, ezp, ezpp, ezppp, kz;
      double somefactor, somefactor_base;

      in >> tmpString >> accuracy_m;
      in >> tmpDouble >> tmpDouble >> tmpInt;
      num_gridpzp = (int)floor(++num_gridpz_m/2.) + 1;
      in >> tmpDouble >> tmpDouble >> tmpInt;

      realFourCoefs_m = new double[accuracy_m];
      imagFourCoefs_m = new double[accuracy_m - 1];

      RealValues = (double*) fftw_malloc(sizeof(double) * num_gridpz_m);
      FourCoefs  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_gridpzp);

      for (int i = 0; i < num_gridpz_m; i++)
        in >> RealValues[i];

      p = fftw_plan_dft_r2c_1d(num_gridpz_m, RealValues, FourCoefs, FFTW_ESTIMATE);
      fftw_execute(p);

      realFourCoefs_m[0] = FourCoefs[0][0] / num_gridpz_m;
      for (int i = 1; i < accuracy_m; i++)
        {
          realFourCoefs_m[i] = 2./num_gridpz_m * FourCoefs[i][0];
          imagFourCoefs_m[i - 1] = -2./num_gridpz_m * FourCoefs[i][1];
        }

      fftw_destroy_plan(p);
      fftw_free(FourCoefs);
      fftw_free(RealValues);

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* read in fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;

    }
}

void FM1DMagnetoStatic::freeMap()
{
  if (realFourCoefs_m != NULL)
    {
      Inform msg("FM1DMS ");

      delete[] realFourCoefs_m;
      delete[] imagFourCoefs_m;

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* freed fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;

    }
}

bool FM1DMagnetoStatic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
//   if (realFourCoefs_m == NULL) readMap();

  const double RR2 = R(0)*R(0) + R(1)*R(1);

  const double kz = two_pi * (R(2) - zbegin_m)/length_m;

  double ez = realFourCoefs_m[0];
  double ezp = 0.0;
  double ezpp = 0.0;
  double ezppp = 0.0;
  double somefactor_base, somefactor;

  for (int l = 1; l < accuracy_m ; l++)
    {
      somefactor_base = l * two_pi / length_m;
      somefactor = 1.0;
      ez    +=              ( realFourCoefs_m[l] * cos(kz * l) + imagFourCoefs_m[l] * sin(kz * l));somefactor *= somefactor_base;
      ezp   += somefactor * (-realFourCoefs_m[l] * sin(kz * l) + imagFourCoefs_m[l] * cos(kz * l));somefactor *= somefactor_base;
      ezpp  += somefactor * (-realFourCoefs_m[l] * cos(kz * l) - imagFourCoefs_m[l] * sin(kz * l));somefactor *= somefactor_base;
      ezppp += somefactor * ( realFourCoefs_m[l] * sin(kz * l) - imagFourCoefs_m[l] * cos(kz * l));
    }

  const double BfieldR = -ezp/2. + ezppp * RR2 / 16.;

//   E(0) +=  0.0
//   E(1) +=  0.0;
//   E(2) +=  0.0;
  B(0) += -BfieldR * R(0);
  B(1) += -BfieldR * R(1);
  B(2) += ez - ezpp * RR2 / 4.;

  return false;
}

bool FM1DMagnetoStatic::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{

}

void FM1DMagnetoStatic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
  zBegin = zbegin_m;
  zEnd = zend_m;
  rBegin = rbegin_m;
  rEnd = rend_m;
}

void FM1DMagnetoStatic::swap()
{}

void FM1DMagnetoStatic::getInfo(Inform *msg)
{
  (*msg) << Filename_m << " (1D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m" << endl;
}

void FM1DMagnetoStatic::rescale(double factor)
{
  if (factor != scaleFactor_m)
    {
      zbegin_m *= factor/scaleFactor_m;
      zend_m *= factor/scaleFactor_m;
      rbegin_m *= factor/scaleFactor_m;
      rend_m *= factor/scaleFactor_m;
      length_m *= factor/scaleFactor_m;
      scaleFactor_m = factor;
    }
}

double FM1DMagnetoStatic::getFrequency() const
{
  return 0.0;
}

void FM1DMagnetoStatic::setFrequency(double freq)
{}
