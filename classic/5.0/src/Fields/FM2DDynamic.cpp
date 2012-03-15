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

// #ifdef LONGITUDINALFIRST     
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

// #else
//       if (swap_m)
//         for (int i = 0; i < num_gridpz_m; i++)
//           {
//             for (int j = 0; j < num_gridpr_m; j++)
//               in >> FieldstrengthEr_m[i * num_gridpr_m + j] >> FieldstrengthEz_m[i * num_gridpr_m + j] >> FieldstrengthHt_m[i * num_gridpr_m + j] >> tmpDouble;
//             if (fabs(FieldstrengthEz_m[i * num_gridpr_m]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[i * num_gridpr_m]);
//           }
//       else
//         {
//           for (int j = 0; j < num_gridpr_m; j++)
//             for (int i = 0; i < num_gridpz_m; i++)
//               in >> FieldstrengthEz_m[i * num_gridpr_m + j] >> FieldstrengthEr_m[i * num_gridpr_m + j] >> tmpDouble >> FieldstrengthHt_m[i * num_gridpr_m + j];
          
//           for (int i = 0; i < num_gridpz_m; i++)
//             if (fabs(FieldstrengthEz_m[i * num_gridpr_m]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[i * num_gridpr_m]);
//         }
// #endif
      for (int i = 0; i < num_gridpr_m * num_gridpz_m; i++)
        {
          FieldstrengthEz_m[i] /= Ezmax;
          FieldstrengthEr_m[i] /= Ezmax;
          FieldstrengthHt_m[i] /= Ezmax;
        }
//       //ALTERNATIVE: probably slightly faster but uses more memory
//       FieldstrengthEz_m = new double[num_gridpz_m * (num_gridpr_m - 1) * 2];
//       FieldstrengthEr_m = new double[num_gridpz_m * (num_gridpr_m - 1) * 2];
//       FieldstrengthHt_m = new double[num_gridpz_m * (num_gridpr_m - 1) * 2];

// #ifdef LONGITUDINALFIRST     
//       if (swap_m)
//         {
//           for (int i = 0; i < num_gridpz_m; i++)
//             {
//               in >> FieldstrengthEr_m[2 * i] 
//                  >> FieldstrengthEz_m[2 * i] 
//                  >> FieldstrengthHt_m[2 * i] 
//                  >> tmpDouble;
//               if (fabs(FieldstrengthEz_m[2 * i]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[2 * i]);

//               for (int j = 1; j < num_gridpr_m - 1; j++)
//                 {
//                   in >> FieldstrengthEr_m[2 * (i + j * num_gridpz_m)] 
//                      >> FieldstrengthEz_m[2 * (i + j * num_gridpz_m)] 
//                      >> FieldstrengthHt_m[2 * (i + j * num_gridpz_m)] 
//                      >> tmpDouble;
//                   FieldstrengthEr_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] = FieldstrengthEr_m[2 * (i + j * num_gridpz_m)];
//                   FieldstrengthEz_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] = FieldstrengthEz_m[2 * (i + j * num_gridpz_m)];
//                   FieldstrengthHt_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] = FieldstrengthHt_m[2 * (i + j * num_gridpz_m)];
//                 }                
//               in >> FieldstrengthEr_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1] 
//                  >> FieldstrengthEz_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1] 
//                  >> FieldstrengthHt_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1] 
//                  >> tmpDouble;
//             }
//         }
//       else
//         {
//           //j == 0
//           for (int i = 0; i < num_gridpz_m; i++)
//             {
//               in >> FieldstrengthEz_m[2 * i] 
//                  >> FieldstrengthEr_m[2 * i] 
//                  >> tmpDouble 
//                  >> FieldstrengthHt_m[2 * i];
//               if (fabs(FieldstrengthEz_m[2 * i]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[2 * i]);
//             }

//           for (int j = 1; j < num_gridpr_m - 1; j++)
//             for (int i = 0; i < num_gridpz_m; i++)
//               {
//                 in >> FieldstrengthEz_m[2 * (i + j * num_gridpz_m)] 
//                    >> FieldstrengthEr_m[2 * (i + j * num_gridpz_m)] 
//                    >> tmpDouble 
//                    >> FieldstrengthHt_m[2 * (i + j * num_gridpz_m)];
//                 FieldstrengthEz_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] = FieldstrengthEz_m[2 * (i + j * num_gridpz_m)];
//                 FieldstrengthEr_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] = FieldstrengthEr_m[2 * (i + j * num_gridpz_m)];
//                 FieldstrengthHt_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] = FieldstrengthHt_m[2 * (i + j * num_gridpz_m)];
//               }                

//           //j == num_gridpr_m - 1 
//           for (int i = 0; i < num_gridpz_m; i++)
//             in >> FieldstrengthEz_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1] 
//                >> FieldstrengthEr_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1]
//                >> tmpDouble
//                >> FieldstrengthHt_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1];

//         }

//       for (int i = 0; i < 2 * (num_gridpr_m - 1) * num_gridpz_m; i++)
//         {
//           FieldstrengthEz_m[i] /= Ezmax;
//           FieldstrengthEr_m[i] /= Ezmax;
//           FieldstrengthHt_m[i] /= Ezmax;
//         }
// #else
//       FieldstrengthEz_m = new double[num_gridpr_m * (num_gridpz_m - 1) * 2];
//       FieldstrengthEr_m = new double[num_gridpr_m * (num_gridpz_m - 1) * 2];
//       FieldstrengthHt_m = new double[num_gridpr_m * (num_gridpz_m - 1) * 2];

//       if (swap_m)
//         {
//           //j == 0
//           for (int i = 0; i < num_gridpr_m; i++)
//             {
//               in >> FieldstrengthEr_m[2 * i] 
//                  >> FieldstrengthEz_m[2 * i] 
//                  >> FieldstrengthHt_m[2 * i] 
//                  >> tmpDouble;
//             }

//           for (int j = 1; j < num_gridpz_m - 1; j++)
//             for (int i = 0; i < num_gridpr_m; i++)
//               {
//                 in >> FieldstrengthEr_m[2 * (i + j * num_gridpr_m)] 
//                    >> FieldstrengthEz_m[2 * (i + j * num_gridpr_m)] 
//                    >> FieldstrengthHt_m[2 * (i + j * num_gridpr_m)] 
//                    >> tmpDouble;
//                 FieldstrengthEz_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthEz_m[2 * (i + j * num_gridpr_m)];
//                 FieldstrengthEr_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthEr_m[2 * (i + j * num_gridpr_m)];
//                 FieldstrengthHt_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthHt_m[2 * (i + j * num_gridpr_m)];
//               }                

//           //j == num_gridpr_m - 1 
//           for (int i = 0; i < num_gridpr_m; i++)
//             in >> FieldstrengthEr_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1] 
//                >> FieldstrengthEz_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1]
//                >> FieldstrengthHt_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1]
//                >> tmpDouble;

//           for (int i = 0; i < num_gridpz_m; i++)
//             if (fabs(FieldstrengthEz_m[2 * i * num_gridpr_m]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[2 * i * num_gridpr_m]);          

//        }
//       else
//         {
//           for (int i = 0; i < num_gridpr_m; i++)
//             {
//               in >> FieldstrengthEz_m[2 * i] 
//                  >> FieldstrengthEr_m[2 * i] 
//                  >> tmpDouble 
//                  >> FieldstrengthHt_m[2 * i];
              
//               for (int j = 1; j < num_gridpz_m - 1; j++)
//                 {
//                   in >> FieldstrengthEz_m[2 * (i + j * num_gridpr_m)] 
//                      >> FieldstrengthEr_m[2 * (i + j * num_gridpr_m)] 
//                      >> tmpDouble 
//                      >> FieldstrengthHt_m[2 * (i + j * num_gridpr_m)];
//                   FieldstrengthEr_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthEr_m[2 * (i + j * num_gridpr_m)];
//                   FieldstrengthEz_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthEz_m[2 * (i + j * num_gridpr_m)];
//                   FieldstrengthHt_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthHt_m[2 * (i + j * num_gridpr_m)];
//                 }                
//               in >> FieldstrengthEz_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1] 
//                  >> FieldstrengthEr_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1]
//                  >> tmpDouble 
//                  >> FieldstrengthHt_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1];
//             }
//           for (int i = 0; i < num_gridpz_m; i++)
//             if (fabs(FieldstrengthEz_m[2 * i * num_gridpr_m]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[2 * i * num_gridpr_m]);          
 
//         }

//       for (int i = 0; i < 2 * (num_gridpz_m - 1) * num_gridpr_m; i++)
//         {
//           FieldstrengthEz_m[i] /= Ezmax;
//           FieldstrengthEr_m[i] /= Ezmax;
//           FieldstrengthHt_m[i] /= Ezmax;
//         }
// #endif
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
  //  if (FieldstrengthEz_m == NULL) readMap();

  const double RR = sqrt(R(0)*R(0) + R(1)*R(1));

  const int indexr = (int)floor(RR / hr_m);
  const double leverr = (RR / hr_m) - indexr;
  
  const int indexz = (int)floor((R(2)) / hz_m);
  const double leverz = (R(2) / hz_m) - indexz;

  if ((indexz < 0) || (indexz + 2 > num_gridpz_m))
    return false;
  if (indexr + 2 > num_gridpr_m)
    return true;

// #ifdef LONGITUDINALFIRST     
  const int index1 = indexz + indexr * num_gridpz_m;
  const int index2 = index1 + num_gridpz_m;

  if (index2 + 2 > num_gridpr_m * num_gridpz_m)
    {
      Inform msg("FM2DD");
      msg << "* ************** W A R N I N G *****************************************************" << endl;
      msg << "* index2 + 2 = " << index2 + 2 << " > num_gridpr_m * num_gridpz_m = " << num_gridpr_m << " * " << num_gridpz_m << " = " << num_gridpr_m * num_gridpz_m << endl;
      msg << "* **********************************************************************************" << endl;
      return false;
    }
  else
    {
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
      
// #else
//   int index1 = indexz * num_gridpr_m + indexr;
//   int index2 = index1 + num_gridpr_m;

//   double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index1 + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index2]
//                  + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];
  
//   double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index1 + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index2]
//                  + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

//   double HfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthHt_m[index1] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthHt_m[index1 + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthHt_m[index2]
//                  + leverz         * leverr         * FieldstrengthHt_m[index2 + 1];
// #endif
  //ALTERNATIVE: probably faster but uses more memory
// #ifdef LONGITUDINALFIRST
//   int index = 2*(indexz + indexr * num_gridpz_m);
//   double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthEr_m[index + 3];
  
//   double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthEz_m[index + 3];

//   double HfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthHt_m[index] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthHt_m[index + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthHt_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthHt_m[index + 3];
// #else
//   int index = 2*(indexr + indexz * num_gridpr_m);
//   double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index] 
//                  + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index + 1]
//                  + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthEr_m[index + 3];
  
//   double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index] 
//                  + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index + 1]
//                  + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthEz_m[index + 3];

//   double HfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthHt_m[index] 
//                  + leverz         * (1.0 - leverr) * FieldstrengthHt_m[index + 1]
//                  + (1.0 - leverz) * leverr         * FieldstrengthHt_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthHt_m[index + 3];
// #endif
      if( RR != 0 ) {
        E(0) += EfieldR * R(0)/RR;
        E(1) += EfieldR * R(1)/RR;
        B(0) += -mu_0 * HfieldT * R(1) / RR * 1e-6;
        B(1) +=  mu_0 * HfieldT * R(0) / RR * 1e-6;
      }
      E(2) += EfieldZ;
      return false;
    }
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
