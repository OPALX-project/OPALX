// ------------------------------------------------------------------------
// $RCSfile: TravelingWave.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TravelingWave
//   Defines the abstract interface for an accelerating structure.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include <iostream>
#include <fstream>
#include "fftw3.h"


// Class TravelingWave
// ------------------------------------------------------------------------

TravelingWave::TravelingWave():
  Component(),
  NumCells_m(0),
  lengthUnit_m(1.0),
  fast_m(false)
{}


TravelingWave::TravelingWave(const TravelingWave &right):
  Component(right),
  CoreFilename_m(right.CoreFilename_m),
//   EntryFilename_m(right.EntryFilename_m),
//   ExitFilename_m(right.ExitFilename_m),
  CoreFieldmap_m(right.CoreFieldmap_m),
//   EntryFringeField_m(right.EntryFringeField_m),
//   ExitFringeField_m(right.ExitFringeField_m),
  scale_m(right.scale_m),
  scaleCore_m(right.scaleCore_m),
  frequency_m(right.frequency_m),
  phase_m(right.phase_m),
  phaseCore1_m(right.phaseCore1_m),
  phaseCore2_m(right.phaseCore2_m),
  phaseExit_m(right.phaseExit_m),
  startField_m(right.startField_m),
  startCoreField_m(right.startCoreField_m),
  startExitField_m(right.startExitField_m),
  mappedStartExitField_m(right.mappedStartExitField_m),
  PeriodLength_m(right.PeriodLength_m),
  NumCells_m(right.NumCells_m),
  CellLength_m(right.CellLength_m),
  fast_m(right.fast_m),
  lengthUnit_m(right.lengthUnit_m),
  Mode_m(right.Mode_m)
{}


TravelingWave::TravelingWave(const string &name):
  Component(name)
{}


TravelingWave::~TravelingWave()
{
  Fieldmap::deleteFieldmap(CoreFilename_m);
//   Fieldmap::deleteFieldmap(EntryFilename_m);
//   Fieldmap::deleteFieldmap(ExitFilename_m);
}


void TravelingWave::accept(BeamlineVisitor &visitor) const
{
  visitor.visitTravelingWave(*this);
}

void TravelingWave::setFieldMapFN(string fn)
{
  CoreFilename_m = fn;
}

// void TravelingWave::setEntryFieldMapFN(string fn)
// {
//   EntryFilename_m = fn;
// }

// void TravelingWave::setExitFieldMapFN(string fn)
// {
//   ExitFilename_m = fn;
// }

string TravelingWave::getFieldMapFN() const
{
  return CoreFilename_m;
}

void TravelingWave::setAmplitudem(double vPeak)
{
  scale_m = vPeak;
}

void TravelingWave::setFrequencym(double freq)
{
  frequency_m = freq;
}

void TravelingWave::setPhasem(double phase)
{
  phase_m = phase;
}

void TravelingWave::setNumCells(int NumCells)
{
  NumCells_m = NumCells;
}

void TravelingWave::setFast(bool fast)
{
  fast_m = fast;
}


bool TravelingWave::getFast() const
{
  return fast_m;
}

bool TravelingWave::apply(const int &i, const double &t, double E[], double B[])
{
  Vector_t Ev(0,0,0), Bv(0,0,0);
  if (apply(RefPartBunch_m->R[i],Vector_t(0.0),t,Ev,Bv)) return true;

  E[0] = Ev(0); E[1] = Ev(1); E[2] = Ev(2);
  B[0] = Bv(0); B[1] = Bv(1); B[2] = Bv(2);

  return false;
}

bool TravelingWave::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  double tmpcos, tmpsin;
  Vector_t tmpR(RefPartBunch_m->R[i](0)-dx_m, RefPartBunch_m->R[i](1)-dy_m ,RefPartBunch_m->R[i](2) - startField_m - ds_m);
  Vector_t tmpE(0.0,0.0,0.0), tmpB(0.0,0.0,0.0);
  bool out_of_bounds = false;


  if (tmpR(2) < startCoreField_m)
    {
      tmpcos =  scale_m * cos( frequency_m * t + phase_m );
      tmpsin = -scale_m * sin( frequency_m * t + phase_m );

    }
  else if (tmpR(2) < startExitField_m)
    {
      Vector_t tmpE2(0.0,0.0,0.0), tmpB2(0.0,0.0,0.0);
      tmpR(2) -= startCoreField_m;
      const double z = tmpR(2);
      tmpR(2) = tmpR(2) - PeriodLength_m * floor(tmpR(2) / PeriodLength_m);
      tmpR(2) += startCoreField_m;

      tmpcos =  scaleCore_m * cos( frequency_m * t + phaseCore1_m);
      tmpsin = -scaleCore_m * sin( frequency_m * t + phaseCore1_m);
      out_of_bounds = CoreFieldmap_m->getFieldstrength(tmpR,tmpE,tmpB);
      E += tmpcos * tmpE;
      B += tmpsin * tmpB;

      tmpE = 0.0;
      tmpB = 0.0;

      tmpR(2) = z + CellLength_m;
      tmpR(2) = tmpR(2) - PeriodLength_m * floor(tmpR(2) / PeriodLength_m);
      tmpR(2) += startCoreField_m;

      tmpcos =  scaleCore_m * cos( frequency_m * t + phaseCore2_m);
      tmpsin = -scaleCore_m * sin( frequency_m * t + phaseCore2_m);

    }
  else
    {
      tmpcos =  scale_m * cos( frequency_m * t + phaseExit_m );
      tmpsin = -scale_m * sin( frequency_m * t + phaseExit_m );
      tmpR(2) -= mappedStartExitField_m;

    }

  out_of_bounds = out_of_bounds || CoreFieldmap_m->getFieldstrength(tmpR,tmpE,tmpB);
  E += tmpcos * tmpE;
  B += tmpsin * tmpB;

  return out_of_bounds;
}

bool TravelingWave::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B)
{
  double tmpcos, tmpsin;
  Vector_t tmpR(R(0)-dx_m, R(1)-dy_m ,R(2) - startField_m-ds_m);
  Vector_t tmpE(0.0,0.0,0.0), tmpB(0.0,0.0,0.0);
  bool out_of_bounds = false;


  if (tmpR(2) < startCoreField_m)
    {
      tmpcos =  scale_m * cos( frequency_m * t + phase_m );
      tmpsin = -scale_m * sin( frequency_m * t + phase_m );

    }
  else if (tmpR(2) < startExitField_m)
    {
      Vector_t tmpE2(0.0,0.0,0.0), tmpB2(0.0,0.0,0.0);
      tmpR(2) -= startCoreField_m;
      const double z = tmpR(2);
      tmpR(2) = tmpR(2) - PeriodLength_m * floor(tmpR(2) / PeriodLength_m);
      tmpR(2) += startCoreField_m;

      tmpcos =  scaleCore_m * cos( frequency_m * t + phaseCore1_m);
      tmpsin = -scaleCore_m * sin( frequency_m * t + phaseCore1_m);
      out_of_bounds = CoreFieldmap_m->getFieldstrength(tmpR,tmpE,tmpB);
      E += tmpcos * tmpE;
      B += tmpsin * tmpB;

      tmpE = 0.0;
      tmpB = 0.0;

      tmpR(2) = z + CellLength_m;
      tmpR(2) = tmpR(2) - PeriodLength_m * floor(tmpR(2) / PeriodLength_m);
      tmpR(2) += startCoreField_m;

      tmpcos =  scaleCore_m * cos( frequency_m * t + phaseCore2_m);
      tmpsin = -scaleCore_m * sin( frequency_m * t + phaseCore2_m);

    }
  else
    {
      tmpcos =  scale_m * cos( frequency_m * t + phaseExit_m );
      tmpsin = -scale_m * sin( frequency_m * t + phaseExit_m );
      tmpR(2) -= mappedStartExitField_m;

    }

  out_of_bounds = out_of_bounds || CoreFieldmap_m->getFieldstrength(tmpR,tmpE,tmpB);
  E += tmpcos * tmpE;
  B += tmpsin * tmpB;

  return out_of_bounds;

}

void TravelingWave::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  using Physics::pi;
  using Physics::two_pi;

  Inform msg("TravelingWave ");

  RefPartBunch_m = bunch;

  msg << getName() << " using file ";

  CoreFieldmap_m = Fieldmap::getFieldmap(CoreFilename_m,fast_m);
  if (CoreFieldmap_m != NULL)
    {
      CoreFieldmap_m->getInfo(&msg);
      if ((frequency_m - CoreFieldmap_m->getFrequency())/frequency_m > 0.0001)
        {
          msg << "************ WARNING ********************************************************" << endl;
          msg << " FREQUENCY IN INPUT FILE DIFFERENT THAN IN FIELD MAP;" << endl;
          msg << frequency_m << " <> " << CoreFieldmap_m->getFrequency() << "; TAKE ON THE LATTER" << endl;
          msg << "*****************************************************************************" << endl;
          frequency_m = CoreFieldmap_m->getFrequency();
        }

      if (dx_m > 1e-10 || dy_m > 1e-10 || ds_m > 1e-10)
        msg << "misaligned by dx = " << dx_m << ", dy = " << dy_m << ", dz = " << ds_m << endl;

      if (hasAttribute("MODE")){
        Mode_m = getAttribute("MODE");
      } else {
        Mode_m = 1./3.;
        msg << "* ************** W A R N I N G *****************************************************" << endl;
        msg << "* NO MODE GIVEN; 2\\pi/3 MODE ASSUMED;" << endl;
        msg << "* **********************************************************************************" << endl;
      }

      double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
      CoreFieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

      PeriodLength_m = (zEnd - zBegin) / 2.0;
      CellLength_m = PeriodLength_m * Mode_m;

      endField = startField + CellLength_m * NumCells_m + PeriodLength_m / 2.0;
      startField -= PeriodLength_m / 2.0;

      startField_m = startField;
      startCoreField_m = PeriodLength_m / 2.0;
      startExitField_m = startCoreField_m + NumCells_m * CellLength_m;
      mappedStartExitField_m = startExitField_m - 3.0 * PeriodLength_m / 2.0;

      scaleCore_m = scale_m / sin(2.0 * pi * Mode_m);
      phaseCore1_m = phase_m + pi * Mode_m / 2.0;
      phaseCore2_m = phase_m + pi * Mode_m * 1.5;
      phaseExit_m = phase_m + 2.0 * pi * (1.0 - NumCells_m * Mode_m);
      
      if (hasWake()) {
        initWakefunction(*this);
        *gmsg << "TravelingWave initialising wake function" << endl;
      }
    }
  else
    {
      endField = startField;
    }
}

void TravelingWave::finalise()
{}

void TravelingWave::rescaleFieldMap(const double &scaleFactor)
{
	//IFF: still needed?
  startField_m *= scaleFactor/lengthUnit_m;
  dx_m *= scaleFactor/lengthUnit_m;
  dy_m *= scaleFactor/lengthUnit_m;
  ds_m *= scaleFactor/lengthUnit_m;
  startCoreField_m *= scaleFactor/lengthUnit_m;
  startExitField_m *= scaleFactor/lengthUnit_m;
  mappedStartExitField_m *= scaleFactor/lengthUnit_m;
  PeriodLength_m *= scaleFactor/lengthUnit_m;
  CellLength_m *= scaleFactor/lengthUnit_m;
  CoreFieldmap_m->rescale(scaleFactor);
  lengthUnit_m = scaleFactor;

}

bool TravelingWave::bends() const
{
  return false;
}


void TravelingWave::goOnline()
{
  Fieldmap::readMap(CoreFilename_m);
  online_m = true;
}

void TravelingWave::goOffline()
{
  Fieldmap::freeMap(CoreFilename_m);
}

void TravelingWave::getDimensions(double &zBegin, double &zEnd) const 
{
  zBegin = startField_m;
  zEnd = startField_m + NumCells_m * CellLength_m;
}


const string& TravelingWave::getType() const
{
    static const string type("TravelingWave");
    return type;
}

