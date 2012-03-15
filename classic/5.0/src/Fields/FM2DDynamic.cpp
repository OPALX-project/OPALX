#include <fstream>
#include <ios>

#include "Fields/FM2DDynamic.hh"
#include "Physics/Physics.h"
#define LONGITUDINALFIRST
#define MEMORYRESTRICTED

using namespace std;
using Physics::mu_0;

FM2DDynamic::FM2DDynamic(string aFilename)
  :Filename_m(aFilename),
   FieldstrengthEz_m(NULL),
   FieldstrengthEr_m(NULL),
   FieldstrengthHt_m(NULL)
{
  Inform msg("FM2DD ");
  int tmpInt;
  string tmpString;
  double tmpDouble1, tmpDouble2;

  Type = T2DDynamic;
  ifstream file(Filename_m.c_str());

  if (file.good())
    {
      file >> tmpString >> tmpString;
      if (tmpString == "ZX")
        {
          swap_m = true;
          file >> rbegin_m >> rend_m >> num_gridpr_m;
          file >> frequency_m;
          file >> zbegin_m >> zend_m >> num_gridpz_m;
        }
      else if (tmpString == "XZ")
        {
          swap_m = false;
          file >> zbegin_m >> zend_m >> num_gridpz_m;
          file >> frequency_m;
          file >> rbegin_m >> rend_m >> num_gridpr_m;
        }
      else
        cerr << "unknown orientation of 2D dynamic fieldmap" << endl;

      frequency_m *= Physics::two_pi * 1e6;

      rbegin_m /= 100.0;
      rend_m /= 100.0;
      zbegin_m /= 100.0;
      zend_m /= 100.0;

      hr_m = (rend_m - rbegin_m) / num_gridpr_m;
      hz_m = (zend_m - zbegin_m) / num_gridpz_m;

      num_gridpr_m++;
      num_gridpz_m++;

      tmpInt = 0;
      while(!file.eof())
        {
    	  file >> tmpDouble1 >> tmpDouble1 >> tmpDouble1 >> tmpDouble1;
          tmpInt++;
        }

      if (tmpInt - 1 != (num_gridpz_m*num_gridpr_m))
        {
          msg << "* ************** WARNING ***********************************************************" << endl;
          msg << "* " << Filename_m << ":" << endl;
          msg << "* NUMBER OF LINES, " << tmpInt - 1 << " DOES NOT CORRESPOND TO THE NUMBER OF GRID POINTS, " << num_gridpz_m * num_gridpr_m << endl;
          if (floor(tmpInt/num_gridpz_m) == tmpInt/num_gridpz_m)
            {
              num_gridpr_m = static_cast<int>(tmpInt/num_gridpz_m);
              msg << "* SETTING num_gridpr_m TO " << num_gridpr_m - 1<< endl;
            }
          else if (floor(tmpInt/num_gridpr_m) == tmpInt/num_gridpr_m)
            {
              num_gridpz_m = static_cast<int>(tmpInt/num_gridpr_m);
              msg << "* SETTING num_gridpr_m TO " << num_gridpz_m - 1<< endl;
            }
          else if (floor(tmpInt/(num_gridpz_m - 1)) == tmpInt/(num_gridpz_m - 1))
            {
              num_gridpz_m--;
              num_gridpr_m = static_cast<int>(tmpInt/num_gridpz_m);
              msg << "* SETTING num_gridpz_m TO " << num_gridpz_m - 1<< " AND num_gridpr_m TO " << num_gridpr_m - 1<< endl;
            }
          else if (floor(tmpInt/(num_gridpr_m - 1)) == tmpInt/(num_gridpr_m - 1))
            {
              num_gridpr_m--;
              num_gridpz_m = static_cast<int>(tmpInt/num_gridpr_m);
              msg << "* SETTING num_gridpz_m TO " << num_gridpz_m - 1<< " AND num_gridpr_m TO " << num_gridpr_m - 1<< endl;
            }
          else
            {
              msg << "* DELETING THE ELEMENT!" << endl;
              zend_m = zbegin_m;
            }
          msg << "* **********************************************************************************" << endl;
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


FM2DDynamic::~FM2DDynamic()
{
  if (FieldstrengthEz_m != NULL)
    {
      delete[] FieldstrengthEz_m;
      delete[] FieldstrengthEr_m;
      delete[] FieldstrengthHt_m;
    }
}

void FM2DDynamic::readMap()
{
  if (FieldstrengthEz_m == NULL)
    {
      Inform msg("FM2DD ");
      ifstream in(Filename_m.c_str());
      int tmpInt;
      string tmpString;
      double tmpDouble, Ezmax = 0.0;

      in >> tmpString >> tmpString;
      in >> tmpDouble >> tmpDouble >> tmpInt;
      in >> tmpDouble;
      in >> tmpDouble >> tmpDouble >> tmpInt;

      FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
      FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];
      FieldstrengthHt_m = new double[num_gridpz_m * num_gridpr_m];

      if (swap_m)
        for (int i = 0; i < num_gridpz_m; i++)
          {
            for (int j = 0; j < num_gridpr_m; j++)
              in >> FieldstrengthEr_m[i + j * num_gridpz_m] >> FieldstrengthEz_m[i + j * num_gridpz_m] >> FieldstrengthHt_m[i + j * num_gridpz_m] >> tmpDouble;
            if (fabs(FieldstrengthEz_m[i]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[i]);
          }
      else
        {
          for (int j = 0; j < num_gridpr_m; j++)
            for (int i = 0; i < num_gridpz_m; i++)
              in >> FieldstrengthEz_m[i + j * num_gridpz_m] >> FieldstrengthEr_m[i + j * num_gridpz_m] >> tmpDouble >> FieldstrengthHt_m[i + j * num_gridpz_m];

          for (int i = 0; i < num_gridpz_m; i++)
            if (fabs(FieldstrengthEz_m[i]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[i]);
        }

      for (int i = 0; i < num_gridpr_m * num_gridpz_m; i++)
        {
          FieldstrengthEz_m[i] /= Ezmax;
          FieldstrengthEr_m[i] /= Ezmax;
          FieldstrengthHt_m[i] /= Ezmax;
        }
      in.close();

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* read in fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;

    }
}

void FM2DDynamic::freeMap()
{
  if (FieldstrengthEz_m != NULL)
    {
      Inform msg("FM2DD");
      delete[] FieldstrengthEz_m;
      delete[] FieldstrengthEr_m;
      delete[] FieldstrengthHt_m;

      msg << "* ************** I N F O ********************************************************" << endl;
      msg << "* freed fieldmap \"" << Filename_m << "\"" << endl;
      msg << "* **********************************************************************************" << endl;
    }
}

bool FM2DDynamic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
  const double RR = sqrt(R(0)*R(0) + R(1)*R(1));

  const int indexr = (int)floor(RR / hr_m);
  const double leverr = RR / hr_m - indexr;

  const int indexz = (int)floor(R(2) / hz_m);
  const double leverz = R(2) / hz_m - indexz;

  if ((indexz < 0) || (indexz + 2 > num_gridpz_m))
    return false;
  if (indexr + 2 > num_gridpr_m)
    return true;

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

  const double HfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthHt_m[index1]
                       + leverz         * (1.0 - leverr) * FieldstrengthHt_m[index1 + 1]
                       + (1.0 - leverz) * leverr         * FieldstrengthHt_m[index2]
                       + leverz         * leverr         * FieldstrengthHt_m[index2 + 1];

  if( RR > 1e-10 ) {
    E(0) += EfieldR * R(0)/RR;
    E(1) += EfieldR * R(1)/RR;
    B(0) += -mu_0 * HfieldT * R(1) / RR * 1e-6;
    B(1) +=  mu_0 * HfieldT * R(0) / RR * 1e-6;
  }
  E(2) += EfieldZ;
  return false;
}

bool FM2DDynamic::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{

}

void FM2DDynamic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
  zBegin = zbegin_m;
  zEnd = zend_m;
  rBegin = rbegin_m;
  rEnd = rend_m;
}

void FM2DDynamic::swap()
{
  if (swap_m) swap_m = false;
  else swap_m = true;
}

void FM2DDynamic::getInfo(Inform *msg)
{
  (*msg) << Filename_m << " (2D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

void FM2DDynamic::rescale(double factor)
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

double FM2DDynamic::getFrequency() const
{
  return frequency_m;
}

void FM2DDynamic::setFrequency(double freq)
{
  frequency_m = freq;
}
