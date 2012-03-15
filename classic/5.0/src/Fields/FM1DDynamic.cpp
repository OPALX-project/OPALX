#include <fstream>
#include <ios>

#include "Fields/FM1DDynamic.hh"
#include "Physics/Physics.h"
#include "fftw3.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

#define LONGITUDINALFIRST

FM1DDynamic::FM1DDynamic(string aFilename)
  :Filename_m(aFilename),
   realFourCoefs_m(NULL),
   imagFourCoefs_m(NULL)
{
  Inform msg("FM1DD ");

  Type = T1DDynamic;

  ifstream file(Filename_m.c_str());

  if (file.good())
    {
      int tmpInt;
      string tmpString;
      double tmpDouble;

      file >> tmpString >> accuracy_m;
      file >> zbegin_m >> zend_m >> num_gridpz_m;
      file >> frequency_m;
      file >> rbegin_m >> rend_m >> tmpInt;

      frequency_m *= two_pi * 1e6;
      xlrep_m = frequency_m / c;

      rbegin_m /= 100.;
      rend_m /= 100.;
      zbegin_m /= 100.;
      zend_m /= 100.;

      length_m = zend_m - zbegin_m;

      num_gridpz_m++;

      tmpInt = 0;
      while(!file.eof())
        {
          file >> tmpDouble;
          tmpInt++;
        }

      if (tmpInt - 1 != num_gridpz_m)
        {
          msg << "* ************** WARNING ***********************************************************" << endl;
          msg << "* " << Filename_m << ":" << endl;
          msg << "* NUMBER OF LINES " << tmpInt - 1 << " DOES NOT CORRESPOND TO THE NUMBER OF GRIDPOINTS, " << num_gridpz_m << endl;
          msg << "* TAKING THE LATTER!" << endl;
          msg << "* **********************************************************************************" << endl;
          num_gridpz_m = tmpInt - 1;
        }

      file.close();
    }
  else
    {
      msg << "* ************** W A R N I N G *****************************************************" << endl;
      msg << "* NO FILE \"" << Filename_m << "\" found!" << endl;
      msg << "* **********************************************************************************" << endl;
    }
}

FM1DDynamic::~FM1DDynamic()
{
  if (realFourCoefs_m != NULL)
    {
      delete[] realFourCoefs_m;
      delete[] imagFourCoefs_m;
    }
}

void FM1DDynamic::readMap()
{
  if (realFourCoefs_m == NULL)
    {
      Inform msg("FM1DD ");
      ifstream in(Filename_m.c_str());

      int tmpInt;
      string tmpString;
      double tmpDouble;

      int num_gridpzp;
      double *RealValues;
      fftw_complex* FourCoefs;
      fftw_plan p;
      double ez, ezp, ezpp, ezppp, kz, f, fp, xlrep;
      double somefactor, somefactor_base;

      in >> tmpString >> accuracy_m;
      in >> tmpDouble >> tmpDouble >> tmpInt;
      num_gridpzp = (int)floor(num_gridpz_m/2.) + 1;
      in >> tmpDouble;
      xlrep = tmpDouble * 1e6 * two_pi / c;
      in >> tmpDouble >> tmpDouble >> tmpInt;

      realFourCoefs_m = new double[accuracy_m];
      imagFourCoefs_m = new double[accuracy_m - 1];

      RealValues = (double*) fftw_malloc(sizeof(double) * num_gridpz_m);
      FourCoefs  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_gridpzp);

      for (int i = 0; i < num_gridpz_m; i++)
        in >> RealValues[i];
      in.close();

      p = fftw_plan_dft_r2c_1d(num_gridpz_m, RealValues, FourCoefs, FFTW_ESTIMATE);
      fftw_execute(p);

      realFourCoefs_m[0] = FourCoefs[0][0] / num_gridpz_m;
      for (int i = 1; i < accuracy_m; i++)
        {
          realFourCoefs_m[i] = 2./num_gridpz_m * FourCoefs[i][0];
          imagFourCoefs_m[i - 1] = 2./num_gridpz_m * FourCoefs[i][1];
        }

      fftw_destroy_plan(p);
      fftw_free(FourCoefs);
      fftw_free(RealValues);

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* read in fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;

    }
}

void FM1DDynamic::freeMap()
{
  if (realFourCoefs_m != NULL)
    {
      Inform msg("FM1DD ");

      delete[] realFourCoefs_m;
      delete[] imagFourCoefs_m;

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* freed fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;

    }
}

bool FM1DDynamic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
  //  if (realFourCoefs_m == NULL) readMap();

  const double RR2 = R(0)*R(0) + R(1)*R(1);

  const double kz = two_pi * (R(2) - zbegin_m)/length_m;

  double ez = realFourCoefs_m[0];
  double ezp = 0.0;
  double ezpp = 0.0;
  double ezppp = 0.0;
  double somefactor_base, somefactor;

  for (int l = 1; l < accuracy_m ; l++)
    {
      somefactor_base = l * two_pi / length_m;       // = \frac{dkz}{dz}
      somefactor = 1.0;
      const double coskzl = cos(kz * l);
      const double sinkzl = sin(kz * l);
      ez    +=              ( realFourCoefs_m[l] * coskzl + imagFourCoefs_m[l-1] * sinkzl); somefactor *= somefactor_base;
      ezp   += somefactor * (-realFourCoefs_m[l] * sinkzl + imagFourCoefs_m[l-1] * coskzl); somefactor *= somefactor_base;
      ezpp  += somefactor * (-realFourCoefs_m[l] * coskzl - imagFourCoefs_m[l-1] * sinkzl); somefactor *= somefactor_base;
      ezppp += somefactor * ( realFourCoefs_m[l] * sinkzl - imagFourCoefs_m[l-1] * coskzl);
    }
  const double f  = -(ezpp  + ez *  xlrep_m * xlrep_m)/16.;
  const double fp = -(ezppp + ezp * xlrep_m * xlrep_m)/16.;

  const double EfieldR = -(ezp/2. + fp * RR2);
  const double BfieldT =  (ez/2. + f * RR2) * xlrep_m / c;

  E(0) +=  EfieldR * R(0);
  E(1) +=  EfieldR * R(1);
  E(2) +=  ez + 4 * f * RR2;
  B(0) += -BfieldT * R(1);
  B(1) +=  BfieldT * R(0);
//   B(2) += 0.0;
  return false;
}

bool FM1DDynamic::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{

}

void FM1DDynamic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
  zBegin = zbegin_m;
  zEnd = zend_m;
  rBegin = rbegin_m;
  rEnd = rend_m;
}

void FM1DDynamic::swap()
{}

void FM1DDynamic::getInfo(Inform *msg)
{
  (*msg) << Filename_m << " (1D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m" << endl;
}

void FM1DDynamic::rescale(double factor)
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

double FM1DDynamic::getFrequency() const
{
  return frequency_m;
}

void FM1DDynamic::setFrequency(double freq)
{
  frequency_m = freq;
}
