#include <fstream>
#include <ios>

#include "Fields/FM2DMagnetoStatic.hh"

using namespace std;

FM2DMagnetoStatic::FM2DMagnetoStatic(string aFilename)
  :Filename_m(aFilename),
   FieldstrengthBz_m(NULL),
   FieldstrengthBr_m(NULL)
{
  Inform msg("FM2DMS ");
  int tmpInt;
  string tmpString;
  double tmpDouble1, tmpDouble2;

  Type = T2DMagnetoStatic;
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
        cerr << "unknown orientation of 2D magnetostatic fieldmap" << endl;

      
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

FM2DMagnetoStatic::~FM2DMagnetoStatic()
{
  if (FieldstrengthBz_m != NULL)
    {
      delete[] FieldstrengthBz_m;
      delete[] FieldstrengthBr_m;
    }
}

void FM2DMagnetoStatic::readMap()
{
  if (FieldstrengthBz_m == NULL)
    {
      Inform msg("FM2DMS ");
      ifstream in(Filename_m.c_str());
      int tmpInt;
      string tmpString;
      double tmpDouble, Bzmax = 0.0;

      in >> tmpString >> tmpString;
      in >> tmpDouble >> tmpDouble >> tmpInt;
      in >> tmpDouble >> tmpDouble >> tmpInt;

      FieldstrengthBz_m = new double[num_gridpz_m * num_gridpr_m];
      FieldstrengthBr_m = new double[num_gridpz_m * num_gridpr_m];
      
// #ifdef LONGITUDINALFIRST     
      if (swap_m)
        for (int i = 0; i < num_gridpz_m; i++)
          {
            for (int j = 0; j < num_gridpr_m; j++)
              in >> FieldstrengthBr_m[i + j * num_gridpz_m] >> FieldstrengthBz_m[i + j * num_gridpz_m];
            
            if (fabs(FieldstrengthBz_m[i]) > Bzmax) Bzmax = fabs(FieldstrengthBz_m[i]);
          }
      else
        {
          for (int j = 0; j < num_gridpr_m; j++)
            for (int i = 0; i < num_gridpz_m; i++)
              in >> FieldstrengthBz_m[i + j * num_gridpz_m] >> FieldstrengthBr_m[i + j * num_gridpz_m];
          for (int i = 0; i < num_gridpz_m; i++)
            if (fabs(FieldstrengthBz_m[i]) > Bzmax) Bzmax = fabs(FieldstrengthBz_m[i]);
        }
// #else
//       if (swap_m)
//         for (int i = 0; i < num_gridpz_m; i++)
//           {
//             for (int j = 0; j < num_gridpr_m; j++)
//               in >> FieldstrengthBr_m[i * num_gridpr_m + j] >> FieldstrengthBz_m[i * num_gridpr_m + j];
            
//             if (fabs(FieldstrengthBz_m[i * num_gridpr_m]) > Bzmax) Bzmax = fabs(FieldstrengthBz_m[i * num_gridpr_m]);
//           }
//       else
//         {
//           for (int j = 0; j < num_gridpr_m; j++)
//             for (int i = 0; i < num_gridpz_m; i++)
//               in >> FieldstrengthBz_m[i * num_gridpr_m + j] >> FieldstrengthBr_m[i * num_gridpr_m + j];
//           for (int i = 0; i < num_gridpz_m; i++)
//             if (fabs(FieldstrengthBz_m[i * num_gridpr_m]) > Bzmax) Bzmax = fabs(FieldstrengthBz_m[i * num_gridpr_m]);
//         }
// #endif          


//       //ALTERNATIVE: probably slightly faster but uses more memory
// #ifdef LONGITUDINALFIRST     
//       FieldstrengthBz_m = new double[num_gridpz_m * (num_gridpr_m - 1) * 2];
//       FieldstrengthBr_m = new double[num_gridpz_m * (num_gridpr_m - 1) * 2];
      
//       if (swap_m)
//         {
//           for (int i = 0; i < num_gridpz_m; i++)
//             {
//               in >> FieldstrengthBr_m[2 * i] >> FieldstrengthBz_m[2 * i];
//               if (fabs(FieldstrengthBz_m[2 * i]) > Bzmax) Bzmax = fabs(FieldstrengthBz_m[2 * i]);
//               for (int j = 1; j < num_gridpr_m - 1; j++)
//                 {
//                   in >> FieldstrengthBr_m[2 * (i + j * num_gridpz_m)] >> FieldstrengthBz_m[2 * (i + j * num_gridpz_m)];
//                   FieldstrengthBr_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] = FieldstrengthBr_m[2 * (i + j * num_gridpz_m)];
//                   FieldstrengthBz_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] = FieldstrengthBz_m[2 * (i + j * num_gridpz_m)];
//                 }                
//               in >> FieldstrengthBr_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1] >> FieldstrengthBz_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1];
//             }
//         }
//       else
//         {
//           //j == 0
//           for (int i = 0; i < num_gridpz_m; i++)
//             {
//               in >> FieldstrengthBz_m[2 * i] >> FieldstrengthBr_m[2 * i];
//               if (fabs(FieldstrengthBz_m[2 * i]) > Bzmax) Bzmax = fabs(FieldstrengthBz_m[2 * i]);
//             }

//           for (int j = 1; j < num_gridpr_m - 1; j++)
//             for (int i = 0; i < num_gridpz_m; i++)
//               {
//                 in >> FieldstrengthBz_m[2 * (i + j * num_gridpz_m)] >> FieldstrengthBr_m[2 * (i + j * num_gridpz_m)];
//                 FieldstrengthBz_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] =  FieldstrengthBz_m[2 * (i + j * num_gridpz_m)];
//                 FieldstrengthBr_m[2 * (i + j * num_gridpz_m - num_gridpz_m) + 1] = FieldstrengthBr_m[2 * (i + j * num_gridpz_m)];
//               }                

//           //j == num_gridpr_m - 1 
//           for (int i = 0; i < num_gridpz_m; i++)
//             in >> FieldstrengthBz_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1] >> FieldstrengthBr_m[2 * (i + (num_gridpr_m - 2) * num_gridpz_m) + 1];
//         }

//       for (int i = 0; i < 2 * (num_gridpr_m - 1) * num_gridpz_m; i++)
//         {
//           FieldstrengthBz_m[i] /= Bzmax;
//           FieldstrengthBr_m[i] /= Bzmax;
//         }
// #else
//       FieldstrengthBz_m = new double[num_gridpr_m * (num_gridpz_m - 1) * 2];
//       FieldstrengthBr_m = new double[num_gridpr_m * (num_gridpz_m - 1) * 2];

//       if (swap_m)
//         {
//           for (int i = 0; i < num_gridpr_m; i++)
//             {
//               in >> FieldstrengthBr_m[2 * i] >> FieldstrengthBz_m[2 * i];
//             }
//           for (int j = 1; j < num_gridpz_m - 1; j++)
//             for (int i = 0; i < num_gridpr_m; i++)
//             {
//               in >> FieldstrengthBr_m[2 * (i + j * num_gridpr_m)] >> FieldstrengthBz_m[2 * (i + j * num_gridpr_m)];
//               FieldstrengthBr_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthBr_m[2 * (i + j * num_gridpr_m)];
//               FieldstrengthBz_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthBz_m[2 * (i + j * num_gridpr_m)];
//             }                
//           for (int i = 0; i < num_gridpr_m; i++)
//             {
//               in >> FieldstrengthBr_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1] >> FieldstrengthBz_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1];
//             }
//         }
//       else
//         {
//           //j == 0
//           for (int i = 0; i < num_gridpr_m; i++)
//             {
//               in >> FieldstrengthBz_m[2 * i] >> FieldstrengthBr_m[2 * i];

//               for (int j = 1; j < num_gridpz_m - 1; j++)
//                 {
//                   in >> FieldstrengthBz_m[2 * (i + j * num_gridpr_m)] >> FieldstrengthBr_m[2 * (i + j * num_gridpr_m)];
//                   FieldstrengthBz_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthBz_m[2 * (i + j * num_gridpr_m)];
//                   FieldstrengthBr_m[2 * (i + j * num_gridpr_m - num_gridpr_m) + 1] = FieldstrengthBr_m[2 * (i + j * num_gridpr_m)];
//                 }                
//               in >> FieldstrengthBz_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1] >> FieldstrengthBr_m[2 * (i + (num_gridpz_m - 2) * num_gridpr_m) + 1];
//             }
//         }
//       for (int i = 0; i < num_gridpz_m; i++)
//         if (fabs(FieldstrengthBz_m[2 * i * num_gridpr_m]) > Bzmax) Bzmax = fabs(FieldstrengthBz_m[2 * i * num_gridpr_m]);

//       for (int i = 0; i < 2 * (num_gridpz_m - 1) * num_gridpr_m; i++)
//         {
//           FieldstrengthBz_m[i] /= Bzmax;
//           FieldstrengthBr_m[i] /= Bzmax;
//         }

// #endif  
      in.close();
      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* read in fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;
    } 
}

void FM2DMagnetoStatic::freeMap()
{
  if (FieldstrengthBz_m != NULL)
    {
      Inform msg("FM2DMS ");
      delete[] FieldstrengthBz_m;
      delete[] FieldstrengthBr_m;

      FieldstrengthBz_m = NULL;
      FieldstrengthBr_m = NULL;

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* freed fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;
    }
}

bool FM2DMagnetoStatic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
  //  if (FieldstrengthBz_m == NULL) readMap();

  const double RR = sqrt(R(0)*R(0) + R(1)*R(1));

  const int indexr = abs((int)floor(RR / hr_m));
  const double leverr = (RR / hr_m) - indexr;
  
  const int indexz = abs((int)floor((R(2)) / hz_m));
  const double leverz = (R(2) / hz_m) - indexz;

  if ((indexz < 0) || (indexz + 2 > num_gridpz_m) )
    return false;

  if(indexr + 2 > num_gridpr_m)
    return true;

// #ifdef LONGITUDINALFIRST     
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
// #else
//   int index1 = indexz * num_gridpr_m+ indexr;
//   int index2 = index1 + num_gridpr_m;
//   double BfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBr_m[index1] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthBr_m[index1 + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthBr_m[index2]
//                  + leverz         * leverr         * FieldstrengthBr_m[index2 + 1];
  
//   double BfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBz_m[index1] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthBz_m[index1 + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthBz_m[index2]
//                  + leverz         * leverr         * FieldstrengthBz_m[index2 + 1];

// #endif
//   ALTERNATIVE: probably faster but uses more memory
// #ifdef LONGITUDINALFIRST
//   int index = 2*(indexz + indexr * num_gridpz_m);
//   double BfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBr_m[index] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthBr_m[index + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthBr_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthBr_m[index + 3];
  
//   double BfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBz_m[index] 
//                  + (1.0 - leverz) * leverr         * FieldstrengthBz_m[index + 1]
//                  + leverz         * (1.0 - leverr) * FieldstrengthBz_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthBz_m[index + 3];
// #else
//   int index = 2*(indexr + indexz * num_gridpr_m);
//   double BfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBr_m[index] 
//                  + leverz         * (1.0 - leverr) * FieldstrengthBr_m[index + 1]
//                  + (1.0 - leverz) * leverr         * FieldstrengthBr_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthBr_m[index + 3];
  
//   double BfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBz_m[index] 
//                  + leverz         * (1.0 - leverr) * FieldstrengthBz_m[index + 1]
//                  + (1.0 - leverz) * leverr         * FieldstrengthBz_m[index + 2]
//                  + leverz         * leverr         * FieldstrengthBz_m[index + 3];
// #endif

  if( RR != 0 ) {
    B(0) += BfieldR * R(0)/RR;
    B(1) += BfieldR * R(1)/RR;
  }
  B(2) += BfieldZ;
  return false;
}

bool FM2DMagnetoStatic::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const 
{

  double BfieldR, BfieldZ;

  const double RR = sqrt(R(0)*R(0) + R(1)*R(1));

  const int indexr = (int)floor(RR / hr_m);
  const double leverr = (RR / hr_m) - indexr;
  
  const int indexz = (int)floor((R(2)) / hz_m);
  const double leverz = (R(2) / hz_m) - indexz;

  if ((indexz < 0) || (indexz + 2 > num_gridpz_m) )
    return false;

  if(indexr + 2 > num_gridpr_m)
    return true;

  const int index1 = indexz + indexr * num_gridpz_m;
  const int index2 = index1 + num_gridpz_m;
  if (dir == DZ)
    {
      if (indexz > 0)
        {
          if (indexz < num_gridpz_m - 1)
            {
              BfieldR = (1.0 - leverr) * ((FieldstrengthBr_m[index1-1] - 2. * FieldstrengthBr_m[index1] + FieldstrengthBr_m[index1+1]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                                       + (FieldstrengthBr_m[index1+1] - FieldstrengthBr_m[index1-1]) / (2.*hz_m))
                                   + leverr         * ((FieldstrengthBr_m[index2-1] - 2. * FieldstrengthBr_m[index2] + FieldstrengthBr_m[index2+1]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                                       + (FieldstrengthBr_m[index2+1] - FieldstrengthBr_m[index2-1]) / (2.*hz_m));
              
              BfieldZ = (1.0 - leverr) * ((FieldstrengthBz_m[index1-1] - 2. * FieldstrengthBz_m[index1] + FieldstrengthBz_m[index1+1]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                                       + (FieldstrengthBz_m[index1+1] - FieldstrengthBz_m[index1-1]) / (2.*hz_m))
                                   + leverr         * ((FieldstrengthBz_m[index2-1] - 2. * FieldstrengthBz_m[index2] + FieldstrengthBz_m[index2+1]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                                       + (FieldstrengthBz_m[index2+1] - FieldstrengthBz_m[index2-1]) / (2.*hz_m));
            }
          else
            {
              BfieldR = (1.0 - leverr) * ((FieldstrengthBr_m[index1] - 2. * FieldstrengthBr_m[index1-1] + FieldstrengthBr_m[index1-2]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                          - (3. * FieldstrengthBr_m[index1] - 4. * FieldstrengthBr_m[index1-1] + FieldstrengthBr_m[index1-2]) / (2.*hz_m))
                      + leverr         * ((FieldstrengthBr_m[index2] - 2. * FieldstrengthBr_m[index2-1] + FieldstrengthBr_m[index2-2]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                          - (3. * FieldstrengthBr_m[index2]  - 4. * FieldstrengthBr_m[index1-1] + FieldstrengthBr_m[index2-2]) / (2.*hz_m));
          
              BfieldZ = (1.0 - leverr) * ((FieldstrengthBz_m[index1] - 2. * FieldstrengthBz_m[index1-1] + FieldstrengthBz_m[index1-2]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                          - (3. * FieldstrengthBz_m[index1]  - 4. * FieldstrengthBr_m[index1-1] + FieldstrengthBz_m[index1-2]) / (2.*hz_m))
                      + leverr         * ((FieldstrengthBz_m[index2] - 2. * FieldstrengthBz_m[index2-1] + FieldstrengthBz_m[index2-2]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                          - (3. * FieldstrengthBz_m[index2]  - 4. * FieldstrengthBr_m[index1-1] + FieldstrengthBz_m[index2-2]) / (2.*hz_m));
            }
        }
      else
        {
          BfieldR = (1.0 - leverr) * ((FieldstrengthBr_m[index1] - 2. * FieldstrengthBr_m[index1+1] + FieldstrengthBr_m[index1+2]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                      - (3. * FieldstrengthBr_m[index1] - 4. * FieldstrengthBr_m[index1+1] + FieldstrengthBr_m[index1+2]) / (2.*hz_m))
                  + leverr         * ((FieldstrengthBr_m[index2] - 2. * FieldstrengthBr_m[index2+1] + FieldstrengthBr_m[index2+2]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                      - (3. * FieldstrengthBr_m[index2]  - 4. * FieldstrengthBr_m[index1+1] + FieldstrengthBr_m[index2+2]) / (2.*hz_m));
          
          BfieldZ = (1.0 - leverr) * ((FieldstrengthBz_m[index1] - 2. * FieldstrengthBz_m[index1+1] + FieldstrengthBz_m[index1+2]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                      - (3. * FieldstrengthBz_m[index1]  - 4. * FieldstrengthBr_m[index1+1] + FieldstrengthBz_m[index1+2]) / (2.*hz_m))
                  + leverr         * ((FieldstrengthBz_m[index2] - 2. * FieldstrengthBz_m[index2+1] + FieldstrengthBz_m[index2+2]) * (R(2) - indexz * hz_m) / (hz_m*hz_m)
                                      - (3. * FieldstrengthBz_m[index2]  - 4. * FieldstrengthBr_m[index1+1] + FieldstrengthBz_m[index2+2]) / (2.*hz_m));
        }
    }
  else
    {
      if (indexr > 0)
        {
          const int index_1 = index1 - num_gridpz_m;
          if (indexr < num_gridpr_m - 1)
            {
              BfieldR = (1.0 - leverz) * ((FieldstrengthBr_m[index_1+1] - 2. * FieldstrengthBr_m[index1+1] + FieldstrengthBr_m[index2+1]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                          + (FieldstrengthBr_m[index2+1] - FieldstrengthBr_m[index_1+1]) / (2.*hr_m))
                      + leverz         * ((FieldstrengthBr_m[index_1] - 2. * FieldstrengthBr_m[index1] + FieldstrengthBr_m[index2]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                          + (FieldstrengthBr_m[index2] - FieldstrengthBr_m[index_1]) / (2.*hr_m));
              
              BfieldZ = (1.0 - leverz) * ((FieldstrengthBz_m[index_1+1] - 2. * FieldstrengthBz_m[index1+1] + FieldstrengthBz_m[index2+1]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                          + (FieldstrengthBz_m[index2+1] - FieldstrengthBz_m[index_1+1]) / (2.*hr_m))
                      + leverz         * ((FieldstrengthBz_m[index_1] - 2. * FieldstrengthBz_m[index1] + FieldstrengthBz_m[index2]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                          + (FieldstrengthBz_m[index2] - FieldstrengthBz_m[index_1]) / (2.*hr_m));
            }
          else
            {
              const int index_2 = index_1 - num_gridpz_m;
              BfieldR = (1.0 - leverz) * ((FieldstrengthBr_m[index1+1] - 2. * FieldstrengthBr_m[index_1+1] + FieldstrengthBr_m[index_2+1]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                          - (3. * FieldstrengthBr_m[index1+1] - 4. * FieldstrengthBr_m[index_1+1] + FieldstrengthBr_m[index_2+1]) / (2.*hr_m))
                      + leverz         * ((FieldstrengthBr_m[index1] - 2. * FieldstrengthBr_m[index_1] + FieldstrengthBr_m[index_2]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                          - (3. * FieldstrengthBr_m[index1]  - 4. * FieldstrengthBr_m[index_1] + FieldstrengthBr_m[index_2]) / (2.*hr_m));
          
              BfieldZ = (1.0 - leverz) * ((FieldstrengthBz_m[index1+1] - 2. * FieldstrengthBz_m[index_1+1] + FieldstrengthBz_m[index_2+1]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                          - (3. * FieldstrengthBz_m[index1+1]  - 4. * FieldstrengthBr_m[index_1+1] + FieldstrengthBz_m[index_2+1]) / (2.*hr_m))
                      + leverz         * ((FieldstrengthBz_m[index1] - 2. * FieldstrengthBz_m[index_1] + FieldstrengthBz_m[index_2]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                          - (3. * FieldstrengthBz_m[index1]  - 4. * FieldstrengthBr_m[index_1] + FieldstrengthBz_m[index_2]) / (2.*hr_m));
            }
        }
      else
        {
          const int index3 = index2 + num_gridpz_m;
          BfieldR = (1.0 - leverz) * ((FieldstrengthBr_m[index1+1] - 2. * FieldstrengthBr_m[index2+1] + FieldstrengthBr_m[index3+1]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                      - (3. * FieldstrengthBr_m[index1+1] - 4. * FieldstrengthBr_m[index2+1] + FieldstrengthBr_m[index3+1]) / (2.*hr_m))
                  + leverz         * ((FieldstrengthBr_m[index1] - 2. * FieldstrengthBr_m[index2] + FieldstrengthBr_m[index3]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                      - (3. * FieldstrengthBr_m[index1]  - 4. * FieldstrengthBr_m[index2] + FieldstrengthBr_m[index3]) / (2.*hr_m));
          
          BfieldZ = (1.0 - leverz) * ((FieldstrengthBz_m[index1+1] - 2. * FieldstrengthBz_m[index2+1] + FieldstrengthBz_m[index3+1]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                      - (3. * FieldstrengthBz_m[index1+1]  - 4. * FieldstrengthBr_m[index2+1] + FieldstrengthBz_m[index3+1]) / (2.*hr_m))
                  + leverz         * ((FieldstrengthBz_m[index1] - 2. * FieldstrengthBz_m[index2] + FieldstrengthBz_m[index3]) * (RR - indexr * hr_m) / (hr_m*hr_m)
                                      - (3. * FieldstrengthBz_m[index1]  - 4. * FieldstrengthBr_m[index2] + FieldstrengthBz_m[index3]) / (2.*hr_m));
        }

    }

  if( RR > 1.e-12 ) {
    B(0) += BfieldR * R(0)/RR;
    B(1) += BfieldR * R(1)/RR;
  }
  B(2) += BfieldZ;
  return false;
}

void FM2DMagnetoStatic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
  zBegin = zbegin_m;
  zEnd = zend_m;
  rBegin = rbegin_m;
  rEnd = rend_m;
}

void FM2DMagnetoStatic::swap()
{
  if (swap_m) swap_m = false;
  else swap_m = true;
}

void FM2DMagnetoStatic::getInfo(Inform *msg)
{
  (*msg) << Filename_m << " (2D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" <<endl;
}

void FM2DMagnetoStatic::rescale(double factor)
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

double FM2DMagnetoStatic::getFrequency() const
{
  return 0.0;
}

void FM2DMagnetoStatic::setFrequency(double freq)
{ ;}
