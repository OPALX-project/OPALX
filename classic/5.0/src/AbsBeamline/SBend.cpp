// ------------------------------------------------------------------------
// $RCSfile: SBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: SBend
//   Defines the abstract interface for a sector bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------
 
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include <iostream>
#include <fstream>


// Class SBend
// ------------------------------------------------------------------------

SBend::SBend():
  Component()
{
  setElType(isDipole);
}


SBend::SBend(const SBend &right):
  Component(right)
{
  setElType(isDipole);
}


SBend::SBend(const string &name):
  Component(name)
{
  setElType(isDipole);
}


SBend::~SBend()
{}


void SBend::accept(BeamlineVisitor &visitor) const
{
  visitor.visitSBend(*this);
}


double SBend::getNormalComponent(int n) const
{
  return getField().getNormalComponent(n);
}


double SBend::getSkewComponent(int n) const
{
  return getField().getSkewComponent(n);
}


void SBend::setNormalComponent(int n, double v)
{
  getField().setNormalComponent(n, v);
}


void SBend::setSkewComponent(int n, double v)
{
  getField().setSkewComponent(n, v);
}

bool SBend::apply(const int &i, const double &t, double E[], double B[])
{
  Vector_t Ev(0,0,0), Bv(0,0,0);
  if (apply(RefPartBunch_m->R[i],t,Ev,Bv)) return true;
  
  E[0] = Ev(0); E[1] = Ev(1); E[2] = Ev(2);
  B[0] = Bv(0); B[1] = Bv(1); B[2] = Bv(2);
  
  return false;
}

bool SBend::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  return apply(RefPartBunch_m->R[i],t,E,B);
}

bool SBend::apply( const Vector_t &R, const double &t, Vector_t &E, Vector_t &B)
{
  if (R(2) > z_entry_m && R(2) < z_exit_m)
    {
      if (R(2) > z_entry_m + fringe_length_entry_m && R(2) < z_exit_m - fringe_length_exit_m)
        {
          B(0) += amplitude_m * sin(B_xy_angle_m);
          B(1) += amplitude_m * cos(B_xy_angle_m);
        }else{
          double z_relative = - (R(2) - z_entry_m - enge_shift_entry_m) / gap_width_m;
          double S, dSdz, d2Sdz2 = 0.0;
          double expS, f, dfdz, d2fdz2;

          if (R(2) < z_entry_m + fringe_length_entry_m)
            {
              z_relative = - (R(2) - z_entry_m - enge_shift_entry_m) / gap_width_m;

            }else{
              z_relative = - (z_exit_m - R(2) - enge_shift_exit_m) / gap_width_m;

            }
            
          S = enge_coefficients_m[enge_polynomial_order_m - 1] * z_relative;
          S += enge_coefficients_m[enge_polynomial_order_m - 2];
          dSdz = (enge_polynomial_order_m - 1) * enge_coefficients_m[enge_polynomial_order_m];
          
          for (int i = enge_polynomial_order_m - 3; i >= 0; i--)
            {
              S = S * z_relative + enge_coefficients_m[i];
              dSdz = dSdz * z_relative + (i+1)*enge_coefficients_m[i+1];
              d2Sdz2 = d2Sdz2 * z_relative + (i+2)*(i+1)*enge_coefficients_m[i+2];
            }
     
          f = 1.0 / (1.0 + exp(S));
          dfdz = - f*f * exp(S) * dSdz; // first derivative of f
          d2fdz2 = - f*f * exp(S) * d2Sdz2 + dSdz * dfdz + 2 * dfdz*dfdz / f; // second derivative of f
          
          B(0) += amplitude_m * (dfdz - d2fdz2/2 * R(0)*sin(B_xy_angle_m) * R(0)*sin(B_xy_angle_m));
          B(1) += amplitude_m * (dfdz - d2fdz2/2 * R(1)*cos(B_xy_angle_m) * R(1)*cos(B_xy_angle_m));
          B(2) += amplitude_m * dfdz * R(1);

        }
   
    }

  return false;
      
}

void SBend::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  Inform msg("SBend ");

  RefPartBunch_m = bunch;

  ifstream in(filename_m.c_str());

  if (!in.good())
    {
      msg << "* ************** W A R N I N G *****************************************************" << endl;
      msg << "* NO FILE \"" << filename_m << "\" found!" << endl;
      msg << "* **********************************************************************************" << endl;
    }
  
  double trash;
  in >> z_entry_m >> fringe_length_entry_m >> trash >> enge_shift_entry_m;
  in >> z_exit_m >> fringe_length_exit_m >> trash >> enge_shift_exit_m;
  in >> enge_polynomial_order_m;

  if (enge_coefficients_m) delete[] enge_coefficients_m;

  enge_coefficients_m = new double[enge_polynomial_order_m];

  for (int i = 0; i < enge_polynomial_order_m; i++)
    in >> enge_coefficients_m[i];
    
  endField = startField + z_exit_m/100.0;
  startField += z_entry_m/100.0;

  fringe_length_entry_m /= 100.0;
  fringe_length_exit_m /= 100.0;
  enge_shift_entry_m /= 100.0;
  enge_shift_exit_m /= 100.0;
}

void SBend::finalise()
{}

void SBend::rescaleFieldMap(const double &scaleFactor)
{}

bool SBend::bends() const 
{ return true; }


void SBend::setAmplitudem(double vPeak)
{
  amplitude_m = vPeak;
}

void SBend::setGapWidth(double gapwidth)
{
  gap_width_m = gapwidth;
}

void SBend::setFieldMapFN(string fmapfn)
{
  filename_m = fmapfn;
}

string SBend::getFieldMapFN() const
{
  return filename_m;
}

void SBend::getDimensions(double &zBegin, double &zEnd) const
{

}

double SBend::getR() const
{
  return R_m;
}


const string& SBend::getType() const
{
    static const string type("SBend");
    return type;
}

