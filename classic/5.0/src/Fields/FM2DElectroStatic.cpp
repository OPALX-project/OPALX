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

  Type =  T2DElectroStatic;
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

      tmpInt = 0;
      while(!file.eof())
        {
          file >> tmpDouble1 >> tmpDouble2;
          tmpInt++;
        }
      
      if (tmpInt - 1 != (num_gridpz_m*num_gridpr_m))
        {
          msg << "* ************** WARNING ***********************************************************" << endl;
          msg << "* " << Filename_m << ":" << endl;
          msg << "* NUMBER OF LINES " << tmpInt - 1 << " DOES NOT CORRESPOND TO THE NUMBER OF GRID POINTS, " << num_gridpz_m * num_gridpr_m << endl;
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

FM2DElectroStatic::~FM2DElectroStatic()
{
  if (FieldstrengthEz_m != NULL)
    {
      delete[] FieldstrengthEz_m;
      delete[] FieldstrengthEr_m;
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

void FM2DElectroStatic::freeMap()
{
  if (FieldstrengthEz_m != NULL)
    {
      Inform msg("FM2DES ");
      delete[] FieldstrengthEz_m;
      delete[] FieldstrengthEr_m;

      msg << "* ************** W A R N I N G *****************************************************" << endl;
      msg << "* freed fieldmap \"" << Filename_m << "\"" << endl;
      msg << "* **********************************************************************************" << endl;
    }
}

bool FM2DElectroStatic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
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

bool FM2DElectroStatic::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const 
{

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
