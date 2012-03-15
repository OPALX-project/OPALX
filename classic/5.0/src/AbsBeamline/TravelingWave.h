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

  void setMisalignment(double x, double y, double z);

  void getMisalignment(double &x, double &y, double &z) const;

  bool readFieldMap(double &startField, double &endField, double scaleFactor);

  /// Scale FieldMap with scaleFactor relative to current scaleFactor
  void rescaleFieldMap(double scaleFactor);

  void setAmplitudem(double vPeak);

  void setFrequencym(double freq);

  void setPhasem(double phase);
  
  void setNumCells(int NumCells);
  
  bool getFieldstrength(double R[], double t, double E[], double B[]) const;

  bool getFieldstrength(Vector_t R, double t, Vector_t &E, Vector_t &B) const;
  
 private:
  string CoreFilename_m;             /**< The name of the inputfile*/
/*   string EntryFilename_m; */
/*   string ExitFilename_m; */

  Fieldmap *CoreFieldmap_m;
/*   Fieldmap *EntryFringeField_m; */
/*   Fieldmap *ExitFringeField_m; */

  double scale_m;              /**< scale multiplier*/
  double scaleCore_m;

  double phase_m;              /**< phase shift of timevarying field(degrees)*/
  double phaseCore1_m;
  double phaseCore2_m;
  double phaseExit_m;

  double frequency_m;          /**< Read in frequency of timevarying field(MHz)*/

  double startField_m;
  double startCoreField_m;         /**< startingpoint of field(m)*/
  double startExitField_m;
  double mappedStartExitField_m;

  double PeriodLength_m;
  int NumCells_m;
  double CellLength_m;
  double Mode_m;

  double lengthUnit_m;

  double dx_m;
  double dy_m;
  double dz_m;

  bool fast_m;
  // Not implemented.
  void operator=(const TravelingWave &);
};

#endif // CLASSIC_TravelingWave_HH
