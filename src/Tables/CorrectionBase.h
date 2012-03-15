#ifndef OPAL_CorrectionBase_HH
#define OPAL_CorrectionBase_HH

// ------------------------------------------------------------------------
// $RCSfile: CorrectionBase.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Classes: CorrectionBase
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:22 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "Algorithms/PartData.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Beamlines/TBeamline.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"


#include <list>
using namespace std;

class ElementBase;


// Class CorrectionBase
// ------------------------------------------------------------------------
/// Abstract base class for all orbit correction commands.
//  Factors out all common behaviour for these algorithms.
//  A CorrectionBase contains a RBeamline<CorrectionBase::Row> which
//  holds pointers to the elements in the line, together with the ideal
//  transfer matrices from the beginning to the current position and the
//  orbit computed from the algorithm.  This information is used to find
//  the required corrections.  This class also contains a method to set
//  the orbit corrector strength once it has been found.

class CorrectionBase: public Action {

public:

  virtual ~CorrectionBase();

protected:

  /// The comon attributes for orbit correction commands.
  enum {
    // User-definable attributes:
    LINE,        // The beam line for the table,
    BEAM,        // The beam to be used.
    RANGE,       // The range in the lattice.
    SIZE
  };

  /// Structure for a row of the Twiss table.
  struct Row: public FlaggedElmPtr {

    Row(ElementBase *, int);
    explicit Row(const FlaggedElmPtr &);
    ~Row();

    /// The closed orbit after the element.
    FVector<double,6> orbit;

    /// The transfer matrix up to and including the element.
    FMatrix<double,6,6> matrix;

    /// The accumulated arc length.
    double arc;

    /// Flag telling when observation has been used.
    bool isUsed[2];
  };

  // The contained beamline type.
  typedef TBeamline<Row> TLine;


  /// Exemplar constructor.
  CorrectionBase(int size, const char *name, const char *help);

  /// Clone constructor.
  CorrectionBase(const string &name, CorrectionBase *parent);

  /// Add to kicker strength.
  //  Arguments:
  //  [ol]
  //  [li] The plane: (0 = x, 1 = y).
  //  [li] The position, given by an iterator into the Twiss table.
  //  [li] The kick change.
  //  [/ol]
  void addKick(int plane, Row &, double kick);

  /// List correctors before or after correction.
  void listCorrectors(bool list, int plane);

  /// List monitors before or after correction.
  void listMonitors(bool list, int plane);

  /// Set up the corrector and monitor tables.
  void setupTables();

  /// Routine to test for corrector or monitor.
  void test(ElementBase *);


  /// The flat beam line on which the correction is done.
  TLine itsLine;

  /// The particle reference data.
  PartData reference;

  /// The closed orbit guess.
  FVector<double,6> orbitGuess;

  /// Flags telling wether a corrector exists.
  bool isCorr[2];

  /// Flag telling wether a monitor exists.
  bool isMoni[2];

  // Types used for corrector and monitor tables.
  typedef list<Row *> LocalList;
  typedef LocalList::iterator LocalIter;

  // The table of correctors per plane.
  LocalList correctorTable[2];

  // The table of monitors per plane.
  LocalList monitorTable[2];

private:

  // Not implemented.
  CorrectionBase();
  CorrectionBase(const CorrectionBase &);
  void operator=(const CorrectionBase &);
};

#endif // OPAL_CorrectionBase_HH
