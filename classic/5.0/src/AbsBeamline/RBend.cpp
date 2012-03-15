// ------------------------------------------------------------------------
// $RCSfile: RBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend
//   Defines the abstract interface for a rectangular bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RBend.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include <iostream>
#include <fstream>


// Class RBend
// ------------------------------------------------------------------------

RBend::RBend():
  Component()
{}


RBend::RBend(const RBend &right):
  Component(right),
  filename_m(right.filename_m),
  amplitude_m(right.amplitude_m),
  cos_B_xy_angle_m(right.cos_B_xy_angle_m),
  sin_B_xy_angle_m(right.sin_B_xy_angle_m),
  s_entry_m(right.s_entry_m),
  s_exit_m(right.s_exit_m),
  sin_face_angle_entry_m(right.sin_face_angle_entry_m),
  sin_face_angle_exit_m(right.sin_face_angle_exit_m),
  cos_face_angle_entry_m(right.cos_face_angle_entry_m),
  cos_face_angle_exit_m(right.cos_face_angle_exit_m),
  tan_face_angle_entry_m(right.tan_face_angle_entry_m),
  tan_face_angle_exit_m(right.tan_face_angle_exit_m),
  fringe_length_entry_m(right.fringe_length_entry_m),
  fringe_length_exit_m(right.fringe_length_exit_m),
  gap_width_m(right.gap_width_m),
  enge_coefficients_m(NULL),
  enge_polynomial_order_m(0),
  enge_shift_entry_m(right.enge_shift_entry_m),
  enge_shift_exit_m(right.enge_shift_exit_m),
  wf_func_m(right.wf_func_m)
{}


RBend::RBend(const string &name):
  Component(name)
{
  cos_B_xy_angle_m = 1.0;
  sin_B_xy_angle_m = 0.0;
}


RBend::~RBend()
{}


void RBend::accept(BeamlineVisitor &visitor) const
{
  visitor.visitRBend(*this);
}


double RBend::getNormalComponent(int n) const
{
  return getField().getNormalComponent(n);
}


double RBend::getSkewComponent(int n) const
{
  return getField().getSkewComponent(n);
}


void RBend::setNormalComponent(int n, double v)
{
  getField().setNormalComponent(n, v);
}


void RBend::setSkewComponent(int n, double v)
{
  getField().setSkewComponent(n, v);
}

bool RBend::apply(const int &i, const double &t, double E[], double B[])
{
  Vector_t Ev(0,0,0), Bv(0,0,0);
  if (apply(RefPartBunch_m->R[i],t,Ev,Bv)) return true;
      
  E[0] = Ev(0); E[1] = Ev(1); E[2] = Ev(2);
  B[0] = Bv(0); B[1] = Bv(1); B[2] = Bv(2);
      
  return false;
}

bool RBend::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  const Vector_t &R = RefPartBunch_m->R[i];
  const Vector_t &P = RefPartBunch_m->P[i];
  double fraction = 1.0;

  const double eff_s_entry = tan_face_angle_entry_m * (R(0) * cos_B_xy_angle_m + R(1) * sin_B_xy_angle_m) + s_entry_m;
  const double eff_s_exit = tan_face_angle_exit_m * (R(0) * cos_B_xy_angle_m + R(1) * sin_B_xy_angle_m) + s_exit_m;

  if (R(2) > eff_s_entry && R(2) < eff_s_exit)
    {
      if (R(2) > eff_s_entry + fringe_length_entry_m && R(2) < eff_s_exit - fringe_length_exit_m)
        {
          double recpgamma = 1.0 / sqrt(1. + dot(P,P));
          double dz = P(2) * recpgamma * Physics::c * RefPartBunch_m->dt[i] / 2.0;
          if (fringe_length_entry_m == 0.0 && R(2) - dz < eff_s_entry )
            fraction = (R(2) - eff_s_entry) / dz;
          else if (fringe_length_exit_m == 0.0 && R(2) + dz > eff_s_exit)
            fraction = (eff_s_exit - R(2)) / dz;

          B(0) += fraction * amplitude_m * sin_B_xy_angle_m;
          B(1) += fraction * amplitude_m * cos_B_xy_angle_m;
        }else{
          double s_relative = - (R(2) - s_entry_m - enge_shift_entry_m) / gap_width_m;
          double S, dSdz, d2Sdz2 = 0.0;
          double expS, f, dfdz, d2fdz2;
          double sin_face_angle, cos_face_angle;
          const double y_tilde = R(0)*sin_B_xy_angle_m + R(1)*cos_B_xy_angle_m;

          if (R(2) < eff_s_entry + fringe_length_entry_m)
            {
              s_relative = - ((R(2) - eff_s_entry)*cos_face_angle_entry_m - enge_shift_entry_m) / gap_width_m;
              sin_face_angle = sin_face_angle_entry_m;
              cos_face_angle = cos_face_angle_entry_m;

            }else{
              s_relative = - ((s_exit_m - R(2))*cos_face_angle_exit_m - enge_shift_exit_m) / gap_width_m;
              sin_face_angle = sin_face_angle_exit_m;
              cos_face_angle = cos_face_angle_exit_m;

            }
            
          S = enge_coefficients_m[enge_polynomial_order_m] * s_relative;
          S += enge_coefficients_m[enge_polynomial_order_m - 1];
          dSdz = enge_polynomial_order_m * enge_coefficients_m[enge_polynomial_order_m + 1];
          
          for (int i = enge_polynomial_order_m - 2; i >= 0; i--)
            {
              S = S * s_relative + enge_coefficients_m[i];
              dSdz = dSdz * s_relative + (i+1)*enge_coefficients_m[i+1];
              d2Sdz2 = d2Sdz2 * s_relative + (i+2)*(i+1)*enge_coefficients_m[i+2];
            }
     
          f = 1.0 / (1.0 + exp(S));
          dfdz = - f*f * exp(S) * dSdz; // first derivative of f
          d2fdz2 = - f*f * exp(S) * d2Sdz2 + dSdz * dfdz + 2 * dfdz*dfdz / f; // second derivative of f
          
          B(0) += amplitude_m * (-(f - d2fdz2 * y_tilde * y_tilde) * sin_B_xy_angle_m + dfdz * y_tilde * sin_face_angle * cos_B_xy_angle_m);
          B(1) += amplitude_m * ((f - d2fdz2 * y_tilde * y_tilde) * cos_B_xy_angle_m + dfdz * y_tilde * sin_face_angle * sin_B_xy_angle_m);
          //(dfdz - d2fdz2/2 * R(1)*cos_B_xy_angle_m * R(1)*cos_B_xy_angle_m) * cos_face_angle;
          B(2) += amplitude_m * dfdz * y_tilde * cos_face_angle;

        }
   
    }

  return false;
}
  
bool RBend::apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B)
{
  const double eff_s_entry = tan_face_angle_entry_m * (R(0) * cos_B_xy_angle_m + R(1) * sin_B_xy_angle_m) + s_entry_m;
  const double eff_s_exit = tan_face_angle_exit_m * (R(0) * cos_B_xy_angle_m + R(1) * sin_B_xy_angle_m) + s_exit_m;

  if (R(2) > eff_s_entry && R(2) < eff_s_exit)
    {
      if (R(2) > eff_s_entry + fringe_length_entry_m && R(2) < eff_s_exit - fringe_length_exit_m)
        {
          B(0) += amplitude_m * sin_B_xy_angle_m;
          B(1) += amplitude_m * cos_B_xy_angle_m;
        }else{
          double s_relative = - (R(2) - s_entry_m - enge_shift_entry_m) / gap_width_m;
          double S, dSdz, d2Sdz2 = 0.0;
          double expS, f, dfdz, d2fdz2;
          double sin_face_angle, cos_face_angle;
          const double y_tilde = R(0)*sin_B_xy_angle_m + R(1)*cos_B_xy_angle_m;

          if (R(2) < eff_s_entry + fringe_length_entry_m)
            {
              s_relative = - ((R(2) - eff_s_entry)*cos_face_angle_entry_m - enge_shift_entry_m) / gap_width_m;
              sin_face_angle = sin_face_angle_entry_m;
              cos_face_angle = cos_face_angle_entry_m;

            }else{
              s_relative = - ((s_exit_m - R(2))*cos_face_angle_exit_m - enge_shift_exit_m) / gap_width_m;
              sin_face_angle = sin_face_angle_exit_m;
              cos_face_angle = cos_face_angle_exit_m;

            }
            
          S = enge_coefficients_m[enge_polynomial_order_m] * s_relative;
          S += enge_coefficients_m[enge_polynomial_order_m - 1];
          dSdz = enge_polynomial_order_m * enge_coefficients_m[enge_polynomial_order_m + 1];
          
          for (int i = enge_polynomial_order_m - 2; i >= 0; i--)
            {
              S = S * s_relative + enge_coefficients_m[i];
              dSdz = dSdz * s_relative + (i+1)*enge_coefficients_m[i+1];
              d2Sdz2 = d2Sdz2 * s_relative + (i+2)*(i+1)*enge_coefficients_m[i+2];
            }
     
          f = 1.0 / (1.0 + exp(S));
          dfdz = - f*f * exp(S) * dSdz; // first derivative of f
          d2fdz2 = - f*f * exp(S) * d2Sdz2 + dSdz * dfdz + 2 * dfdz*dfdz / f; // second derivative of f
          
          B(0) += amplitude_m * (-(f - d2fdz2 * y_tilde * y_tilde) * sin_B_xy_angle_m + dfdz * y_tilde * sin_face_angle * cos_B_xy_angle_m);
          B(1) += amplitude_m * ((f - d2fdz2 * y_tilde * y_tilde) * cos_B_xy_angle_m + dfdz * y_tilde * sin_face_angle * sin_B_xy_angle_m);
          //(dfdz - d2fdz2/2 * R(1)*cos_B_xy_angle_m * R(1)*cos_B_xy_angle_m) * cos_face_angle;
          B(2) += amplitude_m * dfdz * R(1) * cos_face_angle;

        }
   
    }

  return false;
      
}

void RBend::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  Inform msg("RBend ");

  RefPartBunch_m = bunch;
  ifstream in(filename_m.c_str());

  if (!in.good())
    {
      msg << "* ************** W A R N I N G *****************************************************" << endl;
      msg << "* NO FILE \"" << filename_m << "\" found!" << endl;
      msg << "* **********************************************************************************" << endl;
    }
  if (gap_width_m <= 0.0)
    {
      msg << "* ************** W A R N I N G *****************************************************" << endl;
      msg << "* GAP HEIGHT ZERO OR NOT SET! PLEASE CHECK ATTRIBUTE HGAP" << endl;
      msg << "* **********************************************************************************" << endl;
    }
  double tmpDouble;
  in >> s_entry_m >> fringe_length_entry_m >> tmpDouble >> enge_shift_entry_m;
  sin_face_angle_entry_m = sin(tmpDouble * Physics::pi / 180.0);
  cos_face_angle_entry_m = cos(tmpDouble * Physics::pi / 180.0);
  tan_face_angle_entry_m = tan(tmpDouble * Physics::pi / 180.0);
  in >> s_exit_m >> fringe_length_exit_m >> tmpDouble >> enge_shift_exit_m;
  sin_face_angle_exit_m = sin(tmpDouble * Physics::pi / 180.0);
  cos_face_angle_exit_m = cos(tmpDouble * Physics::pi / 180.0);
  tan_face_angle_exit_m = tan(tmpDouble * Physics::pi / 180.0);
  in >> enge_polynomial_order_m;

  if (enge_coefficients_m) delete[] enge_coefficients_m;

  enge_coefficients_m = new double[enge_polynomial_order_m + 1];

  for (int i = 0; i < enge_polynomial_order_m + 1; i++)
    in >> enge_coefficients_m[i];
    
  s_entry_m /= 100.0;
  s_exit_m /= 100.0;
  fringe_length_entry_m /= 100.0;
  fringe_length_exit_m /= 100.0;
  enge_shift_entry_m /= 100.0;
  enge_shift_exit_m /= 100.0;

  tmpDouble = startField;
  endField = startField + s_exit_m + fringe_length_exit_m;
  startField += s_entry_m - fringe_length_entry_m;

  s_entry_m += tmpDouble;
  s_exit_m += tmpDouble;
  online_m = true;

  /**
     Here we can initialize the wake function is present
     
  */
  if (wf_func_m) 
    msg << *wf_func_m << endl;
  else
    msg << "No wake function used" << endl;

}

void RBend::finalise()
{
  online_m = false;
}

void RBend::rescaleFieldMap(const double &scaleFactor)
{}

bool RBend::bends() const
{
  return true;
}

void RBend::setAmplitudem(double vPeak)
{
  amplitude_m = vPeak;
}

void RBend::setGapWidth(double gapwidth)
{
  gap_width_m = gapwidth;
}

void RBend::setFaceAngleEntry(double angle)
{
  sin_face_angle_entry_m = sin(angle);
  cos_face_angle_entry_m = cos(angle);
  tan_face_angle_entry_m = tan(angle);
}

void RBend::setFaceAngleExit(double angle)
{
  sin_face_angle_exit_m = sin(angle);
  cos_face_angle_exit_m = cos(angle);
  tan_face_angle_exit_m = tan(angle);
}

void RBend::setFieldMapFN(string fmapfn)
{
  filename_m = fmapfn;
}

string RBend::getFieldMapFN() const
{
  return filename_m;
}

void RBend::setWakeFunction(Wake *wf)
{
  wf_func_m = wf;
}

void RBend::setEngeCoefs(const vector<double> EngeCoefs)
{
  if (enge_polynomial_order_m != EngeCoefs.size())
    {
      enge_polynomial_order_m = EngeCoefs.size();
      if (enge_coefficients_m)
        delete[] enge_coefficients_m;
      enge_coefficients_m = new double[enge_polynomial_order_m];
    }
  int i = 0;
  for (int i = 0; i < EngeCoefs.size(); i++)
    enge_coefficients_m[i] = EngeCoefs[i];
}
