#ifndef CLASSIC_TravelingWave_HH
#define CLASSIC_TravelingWave_HH

// ------------------------------------------------------------------------
// $RCSfile: TravelingWave.h,v $
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


#include "AbsBeamline/Component.h"
#include "Fields/Fieldmap.hh"


// Class TravelingWave
// ------------------------------------------------------------------------
/// Interface for RF cavity.
//  Class TravelingWave defines the abstract interface for RF cavities.


class TravelingWave: public Component {

public:

  enum CavityType { SW, TW };
  /// Constructor with given name.
  explicit TravelingWave(const string &name);

  TravelingWave();
  TravelingWave(const TravelingWave &);
  virtual ~TravelingWave();

  /// Apply visitor to TravelingWave.
  virtual void accept(BeamlineVisitor &) const;

  /// Get RF amplitude.
  virtual double getAmplitude() const = 0;

  /// Get RF frequencey.
  virtual double getFrequency() const = 0;

  /// Get RF phase.
  virtual double getPhase() const = 0;

  /// Set the name of the field map
  void setFieldMapFN(string fmapfn);

/*   void setExitFieldMapFN(string fn); */

/*   void setExitFieldMapFN(string fn); */

  string getFieldMapFN() const;

  void setFast(bool fast);

  bool getFast() const;

  void setAmplitudem(double vPeak);

  void setFrequencym(double freq);

  void setPhasem(double phase);

  void setNumCells(int NumCells);

  virtual bool apply(const int &i, const double &t, double E[], double B[]);

  virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);

  virtual bool apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B);

  virtual void initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

  virtual void finalise();

  virtual void rescaleFieldMap(const double &scaleFactor);

  virtual bool bends() const;

  virtual void goOnline();

  virtual void goOffline();

  virtual string getType() { return "TravelingWave";}

  virtual void getDimensions(double &zBegin, double &zEnd) const;

 private:
  string CoreFilename_m;             /**< The name of the inputfile*/
/*   string EntryFilename_m; */
/*   string ExitFilename_m; */

  Fieldmap *CoreFieldmap_m;
/*   Fieldmap *EntryFringeField_m; */
/*   Fieldmap *ExitFringeField_m; */

  double scale_m;              /**< scale multiplier*/
  double scaleCore_m;

  double phase_m;              /**< phase shift of time varying field(degrees)*/
  double phaseCore1_m;
  double phaseCore2_m;
  double phaseExit_m;

  double frequency_m;          /**< Read in frequency of time varying field(MHz)*/

  double startField_m;
  double startCoreField_m;         /**< starting point of field(m)*/
  double startExitField_m;
  double mappedStartExitField_m;

  double PeriodLength_m;
  int NumCells_m;
  double CellLength_m;
  double Mode_m;

  double lengthUnit_m;

  bool fast_m;

  // Not implemented.
  void operator=(const TravelingWave &);
};

#endif // CLASSIC_TravelingWave_HH
