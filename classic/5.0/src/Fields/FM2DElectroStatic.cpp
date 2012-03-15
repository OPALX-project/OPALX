#include <fstream>
#include <ios>

#include "Fields/FM2DElectroStatic.hh"

using namespace std;

FM2DElectroStatic::FM2DElectroStatic(string aFilename)
  :Filename_m(aFilename),
   FieldstrengthEz_m(NULL),
   FieldstrengthEr_m(NULL)
{
  Inform msg("FM2DES ");
  int tmpInt;
  string tmpString;
  double tmpDouble1, tmpDouble2;

  ifstream file(Filename_m.c_str());

  if (file.good())
    {
      file >> tmpString >> tmpString;
      if (tmpString == "ZX")
        {
          swap_m = true;
          file >> rbegin_m >> rend_m >> num_gridpr_m;
          file >> zbegin_m >> zend_m >> num_gridpz_m;
        }
      else if (tmpString == "XZ")
        {
          swap_m = false;
          file >> zbegin_m >> zend_m >> num_gridpz_m;
          file >> rbegin_m >> rend_m >> num_gridpr_m;
        }          
      else
        cerr << "unknown orientation of 2D electrostatic fieldmap" << endl;
      
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


void FM2DElectroStatic::readMap()
{
  if (FieldstrengthEz_m == NULL)
    {
      ifstream in(Filename_m.c_str());
      int tmpInt;
      string tmpString;
      double tmpDouble, Ezmax = 0.0;

      FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
      FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];

      in >> tmpString >> tmpString;
      in >> tmpDouble >> tmpDouble >> tmpInt;
      in >> tmpDouble >> tmpDouble >> tmpInt;
          
      
      if (swap_m)
        for (int i = 0; i < num_gridpz_m; i++)
          for (int j = 0; j < num_gridpr_m; j++)
            in >> FieldstrengthEr_m[i + j * num_gridpz_m] >> FieldstrengthEz_m[i + j * num_gridpz_m];
            
      else
        for (int j = 0; j < num_gridpr_m; j++)
          for (int i = 0; i < num_gridpz_m; i++)
            in >> FieldstrengthEz_m[i + j * num_gridpz_m] >> FieldstrengthEr_m[i + j * num_gridpz_m];

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

bool FM2DElectroStatic::getFieldstrength(Vector_t R, Vector_t &E, Vector_t &B)
{
  double RR = sqrt(R(0)*R(0) + R(1)*R(1));

  int indexr = (int)floor(RR / hr_m);
  double leverr = (RR / hr_m) - indexr;
  
  int indexz = (int)floor((R(2)) / hz_m);
  double leverz = (R(2) / hz_m) - indexz;

  if ((indexz < 0) || (indexz + 2 > num_gridpz_m) || (indexr < 0) || (indexr + 2 > num_gridpr_m)){
    //cerr << getName() << ".getFieldstrength(): out of boudaries (z,r) = (" << R(2) << "," << RR << ")" << endl ;
    return true;
  }

  int index1 = indexz + indexr * num_gridpz_m;
  int index2 = index1 + num_gridpz_m;
  double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1] 
                 + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index1 + 1]
                 + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index2]
                 + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];
  
  double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1] 
                 + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index1 + 1]
                 + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index2]
                 + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

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
    E(0) += EfieldR * R(0)/RR;
    E(1) += EfieldR * R(1)/RR;
  }
  E(2) += EfieldZ;
  return false;
}

void FM2DElectroStatic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
  zBegin = zbegin_m;
  zEnd = zend_m;
  rBegin = rbegin_m;
  rEnd = rend_m;
}

void FM2DElectroStatic::swap()
{
  if (swap_m) swap_m = false;
  else swap_m = true;
}

double FM2DElectroStatic::GetFrequency()
{
  return 0.0;
}

void FM2DElectroStatic::getInfo(Inform *msg)
{
  (*msg) << Filename_m << " (2D electrostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

void FM2DElectroStatic::rescale(double factor)
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

double FM2DElectroStatic::getFrequency() const
{
  return 0.0;
}

void FM2DElectroStatic::setFrequency(double freq)
{ ;}
