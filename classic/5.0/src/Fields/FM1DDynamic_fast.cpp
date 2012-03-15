#include <fstream>
#include <ios>

#include "Fields/FM1DDynamic_fast.hh"
#include "Physics/Physics.h"
#include "fftw3.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

#define LONGITUDINALFIRST

FM1DDynamic_fast::FM1DDynamic_fast(string aFilename)
  :Filename_m(aFilename),
   FieldstrengthEz_m(NULL),
   FieldstrengthEr_m(NULL),
   FieldstrengthBt_m(NULL)
{
  Inform msg("FM1DD ");

  Type = T1DDynamic;
  ifstream file(Filename_m.c_str());

  if (file.good())
    {
      int tmpInt;
      string tmpString;
      double tmpDouble;

      file >> tmpString >> tmpInt;
      file >> zbegin_m >> zend_m >> num_gridpz_m;
      file >> frequency_m;
      file >> rbegin_m >> rend_m >> num_gridpr_m;

      frequency_m *= Physics::two_pi * 1e6;

      rbegin_m /= 100.;
      rend_m /= 100.;
      zbegin_m /= 100.;
      zend_m /= 100.;

      hr_m = (rend_m - rbegin_m)/num_gridpr_m;
      hz_m = (zend_m - zbegin_m)/num_gridpz_m;

      num_gridpr_m++;
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
          msg << "* NUMBER OF LINES " << tmpInt - 1 << " DOES NOT CORRESPOND TO THE NUMBER OF GRID POINTS, " << num_gridpz_m << endl;
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
      num_gridpz_m = num_gridpr_m = 0;
      zbegin_m = zend_m = 0.0;
      rbegin_m = rend_m = 0.0;
    }
}

FM1DDynamic_fast::~FM1DDynamic_fast()
{
  if (FieldstrengthEz_m != NULL)
    {
      delete[] FieldstrengthEz_m;
      delete[] FieldstrengthEr_m;
      delete[] FieldstrengthBt_m;
    }
}

void FM1DDynamic_fast::readMap()
{
  if (FieldstrengthEz_m == NULL)
    {
      Inform msg("FM1DD ");
      ifstream in(Filename_m.c_str());
      if (in.good() && num_gridpz_m > 0 && num_gridpr_m > 0)
        {
          int tmpInt;
          string tmpString;
          double tmpDouble;

          int accuracy;

          int num_gridpzp = (int)floor(num_gridpz_m/2.) + 1;
          double *RealValues;
          fftw_complex* FourCoefs;
          fftw_plan p;
          double ez, ezp, ezpp, ezppp, kz, f, fp, xlrep;
          double somefactor, somefactor_base;
          double R;

          FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
          FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];
          FieldstrengthBt_m = new double[num_gridpz_m * num_gridpr_m];

          in >> tmpString >> accuracy;
          in >> tmpDouble >> tmpDouble >> tmpInt;
          in >> tmpDouble;
          xlrep = tmpDouble * 1e6 * two_pi / c;
          in >> tmpDouble >> tmpDouble >> tmpInt;

          RealValues = (double*) fftw_malloc(sizeof(double) * num_gridpz_m);
          FourCoefs  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_gridpzp);

          for (int i = 0; i < num_gridpz_m; i++)
            in >> RealValues[i];

          p = fftw_plan_dft_r2c_1d(num_gridpz_m, RealValues, FourCoefs, FFTW_ESTIMATE);
          fftw_execute(p);

          tmpDouble = two_pi / num_gridpz_m;
          FourCoefs[0][0] /= num_gridpz_m;
          for (int i = 1; i <= accuracy; i++)
            {
              FourCoefs[i][0] *= 2./num_gridpz_m;
              FourCoefs[i][1] *= -2./num_gridpz_m;
            }

          for (int i = 0; i < num_gridpz_m; i++)
            {
              kz = tmpDouble * i;
              ez = FourCoefs[0][0];
              ezp = 0.0;
              ezpp = 0.0;
              ezppp = 0.0;
              for (int l = 1; l <= accuracy ; l++)
                {
                  somefactor_base = l * tmpDouble / hz_m;
                  somefactor = 1.0;
                  const double coskzl = cos(kz * l);
                  const double sinkzl = sin(kz * l);
                  ez    +=              ( FourCoefs[l][0] * coskzl + FourCoefs[l][1] * sinkzl); somefactor *= somefactor_base;
                  ezp   += somefactor * (-FourCoefs[l][0] * sinkzl + FourCoefs[l][1] * coskzl); somefactor *= somefactor_base;
                  ezpp  += somefactor * (-FourCoefs[l][0] * coskzl - FourCoefs[l][1] * sinkzl); somefactor *= somefactor_base;
                  ezppp += somefactor * ( FourCoefs[l][0] * sinkzl - FourCoefs[l][1] * coskzl);
                }

              f  = -(ezpp  + ez *  xlrep * xlrep)/16.;
              fp = -(ezppp + ezp * xlrep * xlrep)/16.;
              for (int j = 0; j < num_gridpr_m; j++)
                {
                  R = j * hr_m;
#ifdef LONGITUDINALFIRST
                  FieldstrengthEz_m[i + j * num_gridpz_m] =   ez + 4 * f * R * R;
                  FieldstrengthEr_m[i + j * num_gridpz_m] = -(ezp/2. + fp * R * R);
                  FieldstrengthBt_m[i + j * num_gridpz_m] =  (ez/2. + f * R * R) * xlrep / c;
#else
                  FieldstrengthEz_m[i * num_gridpr_m + j] =   ez + 4 * f * R * R;
                  FieldstrengthEr_m[i * num_gridpr_m + j] = -(ezp/2. + fp * R * R);
                  FieldstrengthBt_m[i * num_gridpr_m + j] =  (ez/2. + f * R * R) * xlrep / c;
#endif
                }
            }

//           ofstream testout("test.out");
//           testout.precision(8);
//           for (int i = 0; i < num_gridpr_m; ++i)
//             {
//               int i = 10;
//               for (int j = 0; j < num_gridpz_m; ++j)
//                 {
//                   testout << setw(15)
//                           << FieldstrengthEz_m[j + i * num_gridpz_m] << "    "
//                           << FieldstrengthEr_m[j + i * num_gridpz_m] << "    "
//                           << FieldstrengthBt_m[j + i * num_gridpz_m] << endl;
//                 }
//             }
//           testout.close();
          fftw_destroy_plan(p);
          fftw_free(FourCoefs);
          fftw_free(RealValues);

          msg << "* *********** I N F O ***********************************************************" << endl;
          msg << "* read in fieldmap \"" << Filename_m  << "\""<< endl;
          msg << "* *******************************************************************************" << endl;

      //ALTERNATIVE: probably slightly faster but uses more memory
//       FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m * 2];
//       FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m * 2];

//       if (swap_m)
//         {
//           for (int i = 0; i < num_gridpz_m; i++)
//             {
//               in >> FieldstrengthEr_m[2 * i] >> FieldstrengthEz_m[2 * i];
//               for (int j = 1; j < num_gridpr_m - 1; j++)
//                 {
//                   in >> FieldstrengthEr_m[2 * (i + j * num_gridpz_m)] >> FieldstrengthEz_m[2 * (i + j * num_gridpz_m)];
//                   in >> FieldstrengthEr_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] >> FieldstrengthEz_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1];
//                 }
//               in >> FieldstrengthEr_m[2 * (i + num_gridpr_m * num_gridpz_m - num_gridpz_m)] >> FieldstrengthEz_m[2 * (i + num_gridpr_m * num_gridpz_m - num_gridpz_m)];
//             }
//         }
//       else
//         {
//           //j == 0
//           for (int i = 0; i < num_gridpz_m; i++)
//             in >> FieldstrengthEz_m[2 * i] >> FieldstrengthEr_m[2 * i];

//           for (int j = 1; j < num_gridpr_m - 1; j++)
//             for (int i = 0; i < num_gridpz_m; i++)
//               {
//                 in >> FieldstrengthEz_m[2 * (i + j * num_gridpz_m)] >> FieldstrengthEr_m[2 * (i + j * num_gridpz_m)];
//                 in >> FieldstrengthEz_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] >> FieldstrengthEr_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1];
//               }

//           //j == num_gridpr_m - 1
//           for (int i = 0; i < num_gridpz_m; i++)
//             in >> FieldstrengthEz_m[2 * (i + num_gridpr_m * num_gridpz_m - num_gridpz_m)] >> FieldstrengthEr_m[2 * (i + num_gridpr_m * num_gridpz_m - num_gridpz_m)];
//         }
          in.close();
        }
    }
}

void FM1DDynamic_fast::freeMap()
{
  if (FieldstrengthEz_m != NULL)
    {
      Inform msg("FM1DD ");
      delete[] FieldstrengthEz_m;
      delete[] FieldstrengthEr_m;
      delete[] FieldstrengthBt_m;

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* freed fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;
    }
}


bool FM1DDynamic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
//   if (FieldstrengthEz_m == NULL) readMap();

  const double RR = sqrt(R(0)*R(0) + R(1)*R(1));

  const int indexr = (int)floor(RR / hr_m);
  const double leverr = (RR / hr_m) - indexr;

  const int indexz = (int)floor((R(2)) / hz_m);
  const double leverz = (R(2) / hz_m) - indexz;

  if ((indexz < 0) || (indexz + 2 > num_gridpz_m) || (indexr < 0) || (indexr + 2 > num_gridpr_m)){
    //cerr << getName() << ".getFieldstrength(): out of boudaries (z,r) = (" << R(2) << "," << RR << ")" << endl ;
    return true;
  }

#ifdef LONGITUDINALFIRST
  const int index1 = indexz + indexr * num_gridpz_m;
  const int index2 = index1 + num_gridpz_m;
  const double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1]
                       + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index1 + 1]
                       + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index2]
                       + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];

  const double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1]
                       + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index1 + 1]
                       + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index2]
                       + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

  const double BfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBt_m[index1]
                       + leverz         * (1.0 - leverr) * FieldstrengthBt_m[index1 + 1]
                       + (1.0 - leverz) * leverr         * FieldstrengthBt_m[index2]
                       + leverz         * leverr         * FieldstrengthBt_m[index2 + 1];
#else
  const int index1 = indexz * num_gridpr_m + indexr;
  const int index2 = index1 + num_gridpr_m;
  const double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1]
                       + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index1 + 1]
                       + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index2]
                       + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];

  const double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1]
                       + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index1 + 1]
                       + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index2]
                       + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

  const double BfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBt_m[index1]
                       + (1.0 - leverz) * leverr         * FieldstrengthBt_m[index1 + 1]
                       + leverz         * (1.0 - leverr) * FieldstrengthBt_m[index2]
                       + leverz         * leverr         * FieldstrengthBt_m[index2 + 1];
#endif
//   //ALTERNATIVE: probably faster but uses more memory
//   int index = 2*(indexz + indexr * num_gridpz_m);
//   double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index]
//                  + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthEr_m[index + 3];

//   double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index]
//                  + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthEz_m[index + 3];

  if( RR != 0 ) {
    E(0) += EfieldR * R(0);
    E(1) += EfieldR * R(1);
    B(0) += -BfieldT * R(1);
    B(1) +=  BfieldT * R(0);
  }
  E(2) += EfieldZ;
//   B(2) += 0.0;
  return false;
}

bool FM1DDynamic_fast::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{

}


void FM1DDynamic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
  zBegin = zbegin_m;
  zEnd = zend_m;
  rBegin = rbegin_m;
  rEnd = rend_m;
}

void FM1DDynamic_fast::swap()
{}

void FM1DDynamic_fast::getInfo(Inform *msg)
{
  (*msg) << Filename_m << " (1D dynamic fast); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m" << endl;
}

void FM1DDynamic_fast::rescale(double factor)
{
  if (factor != scaleFactor_m)
    {
      zbegin_m *= factor/scaleFactor_m;
      zend_m *= factor/scaleFactor_m;
      rbegin_m *= factor/scaleFactor_m;
      rend_m *= factor/scaleFactor_m;
      hz_m *= factor/scaleFactor_m;
      hr_m *= factor/scaleFactor_m;
      scaleFactor_m = factor;
    }
}

double FM1DDynamic_fast::getFrequency() const
{
  return frequency_m;
}

void FM1DDynamic_fast::setFrequency(double freq)
{
  frequency_m = freq;
}
