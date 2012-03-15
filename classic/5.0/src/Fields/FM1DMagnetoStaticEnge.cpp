#include <fstream>
#include <ios>

#include "Fields/FM1DMagnetoStaticEnge.hh"
#include "Physics/Physics.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

FM1DMagnetoStaticEnge::FM1DMagnetoStaticEnge(string aFilename)
  :Filename_m(aFilename),
   EngeCoefs_entry_m(NULL),
   EngeCoefs_exit_m(NULL)
{
  int tmpInt;
  string tmpString;
  double tmpDouble;

  Type = T1DMagnetoStaticEnge;
  ifstream file(Filename_m.c_str());

  if (file.good())
    {
      file >> tmpString >> polynomialOrder_m;
      file >> zbegin_m >> zend_m >> tmpInt;

      zbegin_m /= 100.;
      zend_m /= 100.;

      length_m = zend_m - zbegin_m;

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

FM1DMagnetoStaticEnge::~FM1DMagnetoStaticEnge()
{
  if (EngeCoefs_entry_m != NULL)
    {
      delete[] EngeCoefs_entry_m;
      delete[] EngeCoefs_exit_m;
    }
}

void FM1DMagnetoStaticEnge::readMap()
{
  if (EngeCoefs_entry_m == NULL)
    {
      double tolerance = 1e-8;

      Inform msg("FM1DMS ");
      ifstream in(Filename_m.c_str());

      int tmpInt;
      string tmpString;
      double tmpDouble;
      double minValue = 99999999.99, maxValue = -99999999.99;
      int num_gridpz, num_gridpzp;
      int num_gridp_fringe_entry, num_gridp_fringe_exit;
      int num_gridp_before_fringe_entry, num_gridp_before_fringe_exit;
      double origin_shift_entry, origin_shift_exit;
      double *RealValues;
      double *rightHandSide;
      double *leastSquareMatrix;
      double dZ;

      in >> tmpString >> polynomialOrder_m;
      in >> tmpDouble >> tmpDouble >> num_gridpz;
      dZ = length_m / num_gridpz;

      EngeCoefs_entry_m = new double[polynomialOrder_m + 1];
      EngeCoefs_exit_m = new double[polynomialOrder_m + 1];

      RealValues = new double[num_gridpz + 1];

      for (int i = 0; i < num_gridpz + 1; ++i)
        {
          in >> RealValues[i];
          if (RealValues[i] > maxValue) maxValue = RealValues[i];
          else if (RealValues[i] < minValue) minValue = RealValues[i];
        }

      // normalise the values //
      for (int i = 0; i < num_gridpz + 1; ++i)
        RealValues[i] = (RealValues[i] - minValue)/(maxValue - minValue);

      // find begin of entry fringe field //
      int i = 0;
      while (i < num_gridpz + 1 && RealValues[i] < tolerance) ++i;
      num_gridp_before_fringe_entry = i - 1;
      field_sections_m[0] = dZ * (i - 1);

      // find end of entry fringe field //
      while (i < num_gridpz + 1 && RealValues[i] < 1. - tolerance) ++i;
      num_gridp_fringe_entry = i - 1 - num_gridp_before_fringe_entry;
      fringe_length_entry_m = dZ * num_gridp_fringe_entry;
      field_sections_m[2] = dZ * (i-1);

      // find begin of exit fringe field //
      while (i < num_gridpz + 1 && RealValues[i] >= 1. - tolerance) ++i;
      num_gridp_before_fringe_exit = i - 1;
      field_sections_m[3] = dZ * (i-1);

      while (i < num_gridpz + 1 && RealValues[i] > tolerance) ++i;
      num_gridp_fringe_exit = i - 1 - num_gridp_before_fringe_exit;
      fringe_length_exit_m = dZ * num_gridp_fringe_exit;
      field_sections_m[5] = dZ * (i-1);
      num_gridp_before_fringe_exit = i - 1;

      origin_shift_entry = field_sections_m[1] - field_sections_m[0];
      // set the origin of the polynomials
      // we assume that zbegin is negative and that the origin of the fieldmap coincides 
      // with the origin of the polynomial of the entry fringe field
      // furthermore we assume that the origin of the exit fringe field is as far away 
      // from the end of exit fringe field as the origin of the entry fringe field is away 
      // from the start of the entry fringe field
      field_sections_m[1] = -zbegin_m;
      field_sections_m[4] = field_sections_m[5] - origin_shift_entry;
      origin_shift_exit = field_sections_m[4] - field_sections_m[5];


      if (num_gridp_fringe_entry > num_gridp_fringe_exit)
        {
          leastSquareMatrix = new double[(polynomialOrder_m + 1) * num_gridp_fringe_entry];
          rightHandSide = new double[num_gridp_fringe_entry];
        }
      else
        {
          leastSquareMatrix = new double[(polynomialOrder_m + 1) * num_gridp_fringe_exit];
          rightHandSide = new double[num_gridp_fringe_exit];
        }

      for (int i = 0; i < num_gridp_fringe_entry; ++i)
        {
          double powerOfZ = 1.;
          double Z = dZ * i - origin_shift_entry;
          rightHandSide[i] = log(1./RealValues[num_gridp_before_fringe_entry + i] - 1.);
          for (int j = 0; j < polynomialOrder_m + 1; ++j)
            {
              leastSquareMatrix[i*(polynomialOrder_m + 1) + j] = powerOfZ;
              powerOfZ *= Z;
            }
        }

      QRDecomposition::solve(leastSquareMatrix, EngeCoefs_entry_m, rightHandSide, num_gridp_fringe_entry, polynomialOrder_m + 1);

      for (int i = 0; i < num_gridp_fringe_exit; ++i)
        {
          double powerOfZ = 1.;
          double Z = dZ * i - origin_shift_exit;
          rightHandSide[i] = log(1./RealValues[num_gridp_before_fringe_exit - i] - 1.);
          for (int j = 0; j < polynomialOrder_m + 1; ++j)
            {
              leastSquareMatrix[i*(polynomialOrder_m + 1) + j] = powerOfZ;
              powerOfZ *= Z;
            }
        }

      QRDecomposition::solve(leastSquareMatrix, EngeCoefs_exit_m, rightHandSide, num_gridp_fringe_exit, polynomialOrder_m + 1);

      delete[] RealValues;
      delete[] leastSquareMatrix;
      delete[] rightHandSide;

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* read in fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;

    } 
}

void FM1DMagnetoStaticEnge::freeMap()
{
  if (EngeCoefs_entry_m != NULL)
    {
      Inform msg("FM1DMS ");
      
      delete[] EngeCoefs_entry_m;
      delete[] EngeCoefs_exit_m;

      msg << "* *********** I N F O ***********************************************************" << endl;
      msg << "* freed fieldmap \"" << Filename_m  << "\""<< endl;
      msg << "* *******************************************************************************" << endl;

    }
}

bool FM1DMagnetoStaticEnge::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{ 
  if (R(2) >= field_sections_m[2] && R(2) <= field_sections_m[3])
    {
          B(1) = 1.0;
    }
  else if (R(2) > field_sections_m[0] && R(2) < field_sections_m[2])
    {
      double S, dSdz, d2Sdz2 = 0.0;
      double expS, f, dfdz, d2fdz2;

      const double &z = R(2);
          
      S = EngeCoefs_entry_m[polynomialOrder_m] * z;
      S += EngeCoefs_entry_m[polynomialOrder_m - 1];
      dSdz = polynomialOrder_m * EngeCoefs_entry_m[polynomialOrder_m + 1];
          
      for (int i = polynomialOrder_m - 2; i >= 0; i--)
        {
          S = S * z + EngeCoefs_entry_m[i];
          dSdz = dSdz * z + (i+1)*EngeCoefs_entry_m[i+1];
          d2Sdz2 = d2Sdz2 * z + (i+2)*(i+1)*EngeCoefs_entry_m[i+2];
        }
     
      f = 1.0 / (1.0 + exp(S));
      dfdz = - f*f * exp(S) * dSdz; // first derivative of f
      d2fdz2 = - f*f * exp(S) * d2Sdz2 + dSdz * dfdz + 2 * dfdz*dfdz / f; // second derivative of f
          
      B(0) = f;
      B(1) += dfdz;
      B(2) += d2fdz2;
    }
  else if (R(2) > field_sections_m[3] && R(2) < field_sections_m[5])
    {
      double S, dSdz, d2Sdz2 = 0.0;
      double expS, f, dfdz, d2fdz2;

      double z = zend_m - R(2);
            
      S = EngeCoefs_exit_m[polynomialOrder_m] * z;
      S += EngeCoefs_exit_m[polynomialOrder_m - 1];
      dSdz = polynomialOrder_m * EngeCoefs_exit_m[polynomialOrder_m + 1];
          
      for (int i = polynomialOrder_m - 2; i >= 0; i--)
        {
          S = S * z + EngeCoefs_exit_m[i];
          dSdz = dSdz * z + (i+1)*EngeCoefs_exit_m[i+1];
          d2Sdz2 = d2Sdz2 * z + (i+2)*(i+1)*EngeCoefs_exit_m[i+2];
        }
     
      f = 1.0 / (1.0 + exp(S));
      dfdz = - f*f * exp(S) * dSdz; // first derivative of f
      d2fdz2 = - f*f * exp(S) * d2Sdz2 + dSdz * dfdz + 2 * dfdz*dfdz / f; // second derivative of f
          
      B(0) = f;
      B(1) += dfdz;
      B(2) += d2fdz2;

    }
  return true;

}

void FM1DMagnetoStaticEnge::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
  zBegin = zbegin_m;
  zEnd = zend_m;
}

void FM1DMagnetoStaticEnge::swap()
{}

void FM1DMagnetoStaticEnge::getInfo(Inform *msg)
{
  (*msg) << Filename_m << " (1D magnetostatic [Enge]); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m" << endl;
}

void FM1DMagnetoStaticEnge::rescale(double factor)
{
  if (factor != scaleFactor_m)
    {
      zbegin_m *= factor/scaleFactor_m;
      zend_m *= factor/scaleFactor_m;
      length_m *= factor/scaleFactor_m;
      scaleFactor_m = factor;
    }
}

double FM1DMagnetoStaticEnge::getFrequency() const
{
  return 0.0;
}

void FM1DMagnetoStaticEnge::setFrequency(double freq)
{}

namespace QRDecomposition
{
  void solve(double* Matrix, double* Solution, double* rightHandSide, const int &M, const int &N)
  {
    double sinphi;
    double cosphi;
    double tempValue;
    double len;
    double *R = new double[M * N];
    double *tempVector = new double[M];
    double *residuum = new double[M];
    
    for (int i = 0; i < M; ++i)
      {
        for (int j = 0; j < N; ++j)
          R[i * N + j] = Matrix[i * N + j];
        tempVector[i] = rightHandSide[i];
      }

    /* using Givens rotations */
    for (int i = 0; i < N; ++i)
      {
        for (int j = i+1; j < M; ++j)
          {
            len = sqrt(R[j * N + i] * R[j * N + i] + R[i * (N + 1)] * R[i * (N + 1)]);
            sinphi = -R[j * N + i] / len;
            cosphi = R[i * (N + 1)] / len;
            
            for (int k = 0; k < N; ++k)
              {
                tempValue = cosphi * R[ i * N + k] - sinphi * R[ j * N + k];
                R[j * N + k] = sinphi * R[ i * N + k] + cosphi * R[ j * N + k];
                R[i * N + k] = tempValue;
              }
          }
      }

    /* one step of iterative refinement */

//     cout << "A^T*b" << endl;
    for (int i = 0; i < N; ++i)      /* A^T*b */
      {
        tempValue = 0.0;
        for (int j = 0; j < M; ++j)
          {
            tempValue += Matrix[j * N + i] * rightHandSide[j];
          }
        Solution[i] = tempValue;
      }
//     cout << "R^-TA^T*b" << endl;
    for (int i = 0; i < N; ++i)     /* R^-T*A^T*b */
      {
        tempValue = 0.0;
        for (int j = 0; j < i; ++j)
          tempValue += R[j * N + i] * residuum[j];
        residuum[i] = (Solution[i] - tempValue) / R[i * (N + 1)];
      }
//     cout << "R^-1R^-TA^T*b" << endl;
    for (int i = N - 1; i >= 0; --i) /* R^-1*R^-T*A^T*b */
      {
        tempValue = 0.0;
        for (int j = N - 1; j > i; --j)
          tempValue += R[i * N + j] * Solution[j];
        Solution[i] = (residuum[i] - tempValue) / R[i * (N + 1)];
      }
//     cout << "b - A*x" << endl;
    for (int i = 0; i < M; ++i)
      {
        tempValue = 0.0;
        for (int j = 0; j < N; ++j)
          tempValue += Matrix[i * N + j] * Solution[j];
        residuum[i] = rightHandSide[i] - tempValue;
      }
//     cout << "A^T*r" << endl;
    for (int i = 0; i < N; ++i)
      {
        tempValue = 0.0;
        for (int j = 0; j < M; ++j)
          tempValue += Matrix[j * N + i] * residuum[j];
        tempVector[i] = tempValue;
      }
//     cout << "R^-TA^T*r" << endl;
    for (int i = 0; i < N; ++i)
      {
        tempValue = 0.0;
        for (int j = 0; j < i; ++j)
          tempValue += R[j * N + i] * residuum[j];
        residuum[i] = (tempVector[i] - tempValue)/R[i * (N + 1)];
      }
//     cout << "R^-1R^-TA^T*r" << endl;
    for (int i = N - 1; i >= 0; --i)
      {
        tempValue = 0.0;
        for (int j = N - 1; j > i; --j)
          tempValue += R[i * N + j] * tempVector[j];
        tempVector[i] = (residuum[i] - tempValue) / R[i * (N + 1)];
        Solution[i] += tempVector[i];
      }

    for (int i = 0; i < N; ++i)
       cout << Solution[i] << endl;
    cout << endl;

    delete[] residuum;
    delete[] tempVector;
    delete[] R;
    
  }
  
}
