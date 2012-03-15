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
#include <iostream>
#include <fstream>


// Class RBend
// ------------------------------------------------------------------------

RBend::RBend():
  Component(),
  fast_m(false),
  filename_m(""),
  amplitude_m(0.0),
  field_orientation_m(1.0, 0.0),
  startField_m(0.0),
  length_m(0.0),
  fieldmap_m(NULL),
  sin_face_alpha_m(0.0),
  cos_face_alpha_m(1.0),
  tan_face_alpha_m(0.0),
  sin_face_beta_m(0.0),
  cos_face_beta_m(1.0),
  tan_face_beta_m(0.0),
  map_step_size_m(0.0),
  map_size_m(0),
  map_m(NULL),
  pusher_m(),
  design_energy_m(0.0)
{
  setElType(isDipole);
}


RBend::RBend(const RBend &right):
  Component(right),
  fast_m(right.fast_m),
  filename_m(right.filename_m),
  amplitude_m(right.amplitude_m),
  field_orientation_m(right.field_orientation_m),
  startField_m(right.startField_m),
  length_m(right.length_m),
  fieldmap_m(right.fieldmap_m),
  sin_face_alpha_m(right.sin_face_alpha_m),
  cos_face_alpha_m(right.cos_face_alpha_m),
  tan_face_alpha_m(right.tan_face_alpha_m),
  sin_face_beta_m(right.sin_face_beta_m),
  cos_face_beta_m(right.cos_face_beta_m),
  tan_face_beta_m(right.tan_face_beta_m),
  map_step_size_m(right.map_step_size_m),
  map_size_m(right.map_size_m),
  pusher_m(right.pusher_m),
  design_energy_m(right.design_energy_m)
{
  setElType(isDipole);
  if (map_size_m > 0) {
      map_m = new double[map_size_m + 1];
      for (int i = 0; i < map_size_m + 1; ++i)
        map_m[i] = right.map_m[i];
  } else {
      map_m = NULL;
  }
}


RBend::RBend(const string &name):
  Component(name),
  fast_m(false),
  filename_m(""),
  amplitude_m(0.0),
  field_orientation_m(1.0, 0.0),
  startField_m(0.0),
  length_m(0.0),
  fieldmap_m(NULL),
  sin_face_alpha_m(0.0),
  cos_face_alpha_m(1.0),
  tan_face_alpha_m(0.0),
  sin_face_beta_m(0.0),
  cos_face_beta_m(1.0),
  tan_face_beta_m(0.0),
  map_step_size_m(0.0),
  map_size_m(0),
  map_m(NULL),
  pusher_m(),
  design_energy_m(0.0)
{
  setElType(isDipole);
}


RBend::~RBend()
{
  if (map_m)
    delete[] map_m;
}


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
  if (apply(i,t,Ev,Bv)) return true;
      
  E[0] = Ev(0); E[1] = Ev(1); E[2] = Ev(2);
  B[0] = Bv(0); B[1] = Bv(1); B[2] = Bv(2);
      
  return false;
}

bool RBend::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  const Vector_t &X = RefPartBunch_m->X[i];
  Vector_t strength(0.0), info(0.0);

  fieldmap_m->getFieldstrength(X, strength, info);
  
  if (info(0) > 0.99) {
      
      double y_tilde = field_orientation_m(1) * X(0) + field_orientation_m(0) * X(1);
      B(0) += field_orientation_m(1) * amplitude_m * (strength(0) - strength(2)/2. * y_tilde * y_tilde);
      B(1) += field_orientation_m(0) * amplitude_m * (strength(0) - strength(2)/2. * y_tilde * y_tilde);
      if (info(1) > 0.99) {
          B(2) += amplitude_m * strength(1) * y_tilde;
      } else {
          B(2) -= amplitude_m * strength(1) * y_tilde;
      }
  } else if (fabs(info(0)) < 0.01) {
      B(0) += field_orientation_m(1) * amplitude_m;
      B(1) += field_orientation_m(0) * amplitude_m;
  }

  return false;
}
  
bool RBend::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B)
{
    int index = (int)floor((R(2) - startField_m) / map_step_size_m);
    if (index > 0 && index + 1 < map_size_m) {
        double lever = (R(2) - startField_m) / map_step_size_m - index;
        double z = (1. - lever) * map_m[index] + lever * map_m[index + 1];
        
        Vector_t X = Vector_t(R(0), R(1), z);
        Vector_t strength(0.0), info(0.0);
        fieldmap_m->getFieldstrength(X, strength, info);
        
        if (info(0) > 0.99) {
            double y_tilde = field_orientation_m(1) * X(0) + field_orientation_m(0) * X(1);
            B(0) += field_orientation_m(1) * amplitude_m * (strength(0) - strength(2)/2. * y_tilde * y_tilde);
            B(1) += field_orientation_m(0) * amplitude_m * (strength(0) - strength(2)/2. * y_tilde * y_tilde);
            if (info(1) > 0.99) {
                B(2) +=  amplitude_m * strength(1) * y_tilde;
            } else {
                B(2) -=  amplitude_m * strength(1) * y_tilde;
            }
        } else if (fabs(info(0)) < 0.01) {
            B(0) += field_orientation_m(1) * amplitude_m;
            B(1) += field_orientation_m(0) * amplitude_m;
        }
    }
    return false;
}

void RBend::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  using Physics::c;

  Inform msg("RBend ");

  double tmpDouble;
  double alpha, beta;
  double tolerance = 1e-10;

  RefPartBunch_m = bunch;
  pusher_m.initialise(bunch->getReference());

  *gmsg << "RBend " << getName() << " using file ";
  fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);
  if (fieldmap_m != NULL)
    {
      const double mass = bunch->getM();
      const double gamma = design_energy_m / mass + 1.;
      const double betagamma = sqrt(gamma * gamma - 1.);
      const double dt = bunch->getdT();
      const double charge = bunch->getQ();
      double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
      int j = 0;
      Vector_t tmp(0.0), Bfield(0.0), strength(0.0);
      Vector_t X;
      Vector_t P(-betagamma * sin_face_alpha_m, 0.0, betagamma * cos_face_alpha_m); // TODO: make it 3D
      bool EntryFringe_passed = false;
      double PathLengthEntryFringe = 0.0;  // in S coordinates. This value is different from zBegin due to the curvature!

      fieldmap_m->getInfo(gmsg);
      Fieldmap::readMap(filename_m);
      fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

      X = Vector_t(0.0, 0.0, 0.0);

      length_m = zEnd - zBegin;

      map_step_size_m = betagamma / gamma * c * dt;

      map_size_m = (int)floor(length_m / 2. * Physics::pi / map_step_size_m);
      map_m = new double[map_size_m + 1];
      map_m[0] = 0.0;
      while (map_m[j] < length_m && j < map_size_m)
        {
          strength = Vector_t(0.0);
          X /= Vector_t(c * dt);
          pusher_m.push(X, P, dt);
          X *= Vector_t(c * dt);
          fieldmap_m->getFieldstrength(X, strength, tmp);
          if (X(2) >= fabs(zBegin) && !EntryFringe_passed)
            {
                EntryFringe_passed = true;                       // this is the point where we enter we pass ELEMEDGE
                PathLengthEntryFringe = j * map_step_size_m;     // not the end of the entry fringe field as the name
                                                                 // suggests
            }
          Bfield(1) = amplitude_m * strength(0);
          tmp = Vector_t(0.0);
          X /= Vector_t(c * dt);
          pusher_m.kick(X, P, tmp, Bfield, dt);
          pusher_m.push(X, P, dt);
          X *= Vector_t(c * dt);

          map_m[++j] = X(2);
        }
      map_size_m = j;

      startField -= PathLengthEntryFringe;
      endField = startField + map_step_size_m * j;

      startField_m = startField;
      endField_m = endField;

      R_m = fabs(betagamma * mass / (c * amplitude_m));

      /**
         Here we can initialize the wake function is present
     
      */
      if (hasWake())
        {
          initWakefunction(*this);
          *gmsg << "RBend initialising wake function" << endl;
        }
    }
  else
    {
      endField = startField;
    }

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

void RBend::setFieldMapFN(string fmapfn)
{
  filename_m = fmapfn;
}

string RBend::getFieldMapFN() const
{
  return filename_m;
}

void RBend::getDimensions(double &zBegin, double &zEnd) const
{
  zBegin = startField_m;
  zEnd = endField_m;
}

double RBend::getR() const
{
  return R_m;
}


const string& RBend::getType() const
{
    static const string type("RBend");
    return type;
}

