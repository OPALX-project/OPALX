#include <fstream>
#include <ios>

#include "Fields/FM1DMagnetoStatic_fast.hh"
#include "Physics/Physics.h"
#include "fftw3.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

#define LONGITUDINALFIRST

FM1DMagnetoStatic_fast::FM1DMagnetoStatic_fast(string aFilename)
  :Filename_m(aFilename),
   FieldstrengthBz_m(NULL),
   FieldstrengthBr_m(NULL)
{
  int tmpInt;
  string tmpString;
  double tmpDouble;

  Type = T1DMagnetoStatic;
  ifstream file(Filename_m.c_str());

  if (file.good())
    {
      file >> tmpString >> tmpInt;
      file >> zbegin_m >> zend_m >> num_gridpz_m;
      file >> rbegin_m >> rend_m >> num_gridpr_m;

      rbegin_m /= 100.;
      rend_m /= 100.;
      zbegin_m /= 100.;
      zend_m /= 100.;

      hr_m = (rend_m - rbegin_m)/num_gridpr_m;
      hz_m = (zend_m - zbegin_m)/num_gridpz_m;
      
      num_gridpr_m++;
      num_gridpz_m++;

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


FM1DMagnetoStatic_fast::~FM1DMagnetoStatic_fast()
{
  if (FieldstrengthBz_m != NULL)
    {
      delete[] FieldstrengthBz_m;
      delete[] FieldstrengthBr_m;
    }
}

void FM1DMagnetoStatic_fast::readMap()
{
  if (FieldstrengthBz_m == NULL)
    {
      Inform msg("FM1DMS ");
      ifstream in(Filename_m.c_str());

      int tmpInt;
      string tmpString;
      double tmpDouble;

      int accuracy;

      int num_gridpzp = (int)floor(num_gridpz_m/2.) + 1;
      double *RealValues;
      fftw_complex* FourCoefs;
      fftw_plan p;
      double ez, ezp, ezpp, ezppp, kz;
      double somefactor, somefactor_base;
      double R;

      FieldstrengthBz_m = new double[num_gridpz_m * num_gridpr_m];
      FieldstrengthBr_m = new double[num_gridpz_m * num_gridpr_m];

      in >> tmpString >> accuracy;
      in >> tmpDouble >> tmpDouble >> tmpInt;
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
          FourCoefs[i][1] *= 2./num_gridpz_m;
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
              ez    +=              ( FourCoefs[l][0] * cos(kz * l) + FourCoefs[l][1] * sin(kz * l));somefactor *= somefactor_base;
              ezp   += somefactor * (-FourCoefs[l][0] * sin(kz * l) + FourCoefs[l][1] * cos(kz * l));somefactor *= somefactor_base;
              ezpp  += somefactor * (-FourCoefs[l][0] * cos(kz * l) - FourCoefs[l][1] * sin(kz * l));somefactor *= somefactor_base;
              ezppp += somefactor * ( FourCoefs[l][0] * sin(kz * l) - FourCoefs[l][1] * cos(kz * l));
            }
          for (int j = 0; j < num_gridpr_m; j++)
            {
              R = j * hr_m;
#ifdef LONGITUDINALFIRST
              FieldstrengthBz_m[i + j * num_gridpz_m] =  ez - ezpp * R * R / 4.;
              FieldstrengthBr_m[i + j * num_gridpz_m] = -ezp/2. + ezppp * R * R / 16.;
#else
              FieldstrengthBz_m[i * num_gridpr_m + j] =   ez -ezpp * R * R / 4.;
              FieldstrengthBr_m[i * num_gridpr_m + j] = -ezp/2. + ezppp * R * R / 16.;
#endif
            }
        }

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
//       in.close();
    } 
}

void FM1DMagnetoStatic_fast::freeMap()
{
  if (FieldstrengthBz_m != NULL)
    {
      Inform msg("FM1DD ");
      
      delete[] FieldstrengthBz_m;
      delete[] FieldstrengthBr_m;

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* freed fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;

    }
}

bool FM1DMagnetoStatic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
//   if (FieldstrengthBz_m == NULL) readMap();

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
  const double BfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBr_m[index1] 
                       + leverz         * (1.0 - leverr) * FieldstrengthBr_m[index1 + 1]
                       + (1.0 - leverz) * leverr         * FieldstrengthBr_m[index2]
                       + leverz         * leverr         * FieldstrengthBr_m[index2 + 1];
  
  const double BfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBz_m[index1] 
                       + leverz         * (1.0 - leverr) * FieldstrengthBz_m[index1 + 1]
                       + (1.0 - leverz) * leverr         * FieldstrengthBz_m[index2]
                       + leverz         * leverr         * FieldstrengthBz_m[index2 + 1];

#else
  const int index1 = indexz * num_gridpr_m + indexr;
  const int index2 = index1 + num_gridpr_m;
  const double BfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBr_m[index1] 
                       + (1.0 - leverz) * leverr         * FieldstrengthBr_m[index1 + 1]
                       + leverz         * (1.0 - leverr) * FieldstrengthBr_m[index2]
                       + leverz         * leverr         * FieldstrengthBr_m[index2 + 1];
  
  const double BfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBz_m[index1] 
                       + (1.0 - leverz) * leverr         * FieldstrengthBz_m[index1 + 1]
                       + leverz         * (1.0 - leverr) * FieldstrengthBz_m[index2]
                       + leverz         * leverr         * FieldstrengthBz_m[index2 + 1];

#endif
//   //ALTERNATIVE: probably faster but uses more memory
//   int index = 2*(indexz + indexr * num_gridpz_m);
//   double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBr_m[index] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthBr_m[index + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthBr_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthBr_m[index + 3];
  
//   double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBz_m[index] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthBz_m[index + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthBz_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthBz_m[index + 3];

  E(0) += 0.0;
  E(1) += 0.0;
  E(2) += 0.0;
  B(0) += -BfieldR * R(0);
  B(1) += -BfieldR * R(1);
  B(2) += BfieldZ;
  return false;
}

bool FM1DMagnetoStatic_fast::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const 
{

}


void FM1DMagnetoStatic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
  zBegin = zbegin_m;
  zEnd = zend_m;
  rBegin = rbegin_m;
  rEnd = rend_m;
}

void FM1DMagnetoStatic_fast::swap()
{}

void FM1DMagnetoStatic_fast::getInfo(Inform *msg)
{
  (*msg) << Filename_m << " (1D magnetostatic fast); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m" << endl;
}

void FM1DMagnetoStatic_fast::rescale(double factor)
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

double FM1DMagnetoStatic_fast::getFrequency() const
{
  return 0.0;
}

void FM1DMagnetoStatic_fast::setFrequency(double freq)
{}
