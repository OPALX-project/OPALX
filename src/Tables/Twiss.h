#ifndef OPAL_Twiss_HH
#define OPAL_Twiss_HH

// ------------------------------------------------------------------------
// $RCSfile: Twiss.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Twiss
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:22 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Table.h"
#include "AbstractObjects/Expressions.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Beamlines/TBeamline.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include <iosfwd>
#include <vector>

class AbstractMapper;
class Beam;
class DoomWriter;
class PartData;
class PlaceRep;
class RangeRep;


/// Class Twiss
// ------------------------------------------------------------------------
/// Abstract base class for table buffers holding lattice function.

class Twiss: public Table {

  friend class Insertion;
  friend class Period;

public:

  /// Structure for a row of the Twiss table.
  class Row: public FlaggedElmPtr {

    friend class Insertion;
    friend class Period;
    friend class Twiss;

  public:

    Row(ElementBase *, int);
    explicit Row(const FlaggedElmPtr &);
    ~Row();

    /// Closed orbit.
    const FVector<double,6> &getCO() const;

    /// Transfer matrix.
    const FMatrix<double,6,6> &getMatrix() const;

    /// Arc length.
    double getS() const;

    /// Phase for mode i.
    double getMUi(int i) const;

  private:

    /// The closed orbit after the element.
    FVector<double,6> orbit;

    /// The transfer matrix up to and including the element.
    FMatrix<double,6,6> matrix;

    /// The accumulated arc length.
    double arc;

    /// Phases for the three modes.
    double mu[3];
  };

  // The contained beamline type.
  typedef TBeamline<Row> TLine;

  virtual ~Twiss();

  /// Write the table to the DOOM data base.
  virtual void doomPut(DoomWriter &) const;

  /// Fill in summary record.
  virtual void doomSummary(DoomWriter &) const = 0;

  /// Access to first row.
  //  Version for constant table.
  TLine::const_iterator begin() const;

  /// Access to first row.
  //  Version for non-constant table.
  TLine::iterator begin();

  /// Access to last row.
  //  Version for constant table.
  TLine::const_iterator end() const;

  /// Access to last row.
  //  Version for non-constant table.
  TLine::iterator end();

  /// Check validity of the table definition.
  virtual void execute();

  /// Return a selected value in a selected row.
  virtual double getCell(const PlaceRep &row, const string &col);

  /// Return the default print columns.
  virtual CellArray getDefault() const;

  /// Return column [b]col[/b] of this table, limited by [b]range[/b].
  virtual std::vector<double>
  getColumn(const RangeRep &range, const string &col);

  /// Return current table row in iteration.
  const Row &getCurrent() const;

  /// Return emittance for mode 1.
  double getEX() const;

  /// Return emittance for mode 2.
  double getEY() const;

  /// Return emittance for mode 3.
  double getET() const;

  /// Return the length of the table.
  virtual double getLength();

  /// Return embedded CLASSIC beamline.
  virtual const Beamline *getLine() const;

  /// Return a table row, possible user-defined.
  virtual std::vector<double>
  getRow(const PlaceRep &, const std::vector<string> &);

  /// Check dependency.
  //  Return true, if this table depends on the named object.
  virtual bool isDependent(const string &name) const;

  /// Return column expression.
  //  Return an expression which denotes a column of the twiss table,
  //  identified by its name.
  virtual Expressions::PtrToScalar<double>
  makeColumnExpression(const string &colName) const;

  /// Check compatibility.
  //  True, if [b]rhs[/b] is derived from [b]Twiss[/b].
  virtual bool matches(Table *rhs) const;

  /// Print the body to this TWISS table.
  void printTableBody(std::ostream &, const CellArray &) const;

  /// Print standard information about the TWISS table.
  void printTableTitle(std::ostream &, const char *title) const;

  /// Write TFS body for this table.
  void tfsBody(std::ostream &, const CellArray &) const;

  /// List TFS descriptors for table summary.
  //  Writes the maximum and r.m.s. orbit and dispersion.
  void tfsSummary(std::ostream &) const;

  /// Write TFS descriptors for a TWISS table.
  //  Writes the line and beam info, as well as length and tunes.
  virtual void tfsTwissDescriptors(std::ostream &) const = 0;


  // Access to items in the table.
  // ----------------------------------------------------------------------
  /// Return initial curly A matrix.
  FMatrix<double,6,6> getCurlyA() const;

  /// Curly A map for given row.
  FMatrix<double,6,6> getCurlyA(const Row &) const;

  /// Accumulated transfer map.
  FMatrix<double,6,6> getMatrix(const Row &) const;

  /// Return initial closed orbit.
  FVector<double,6> getOrbit() const;

  /// Get orbit in given row.
  FVector<double,6> getOrbit(const Row &) const;

  /// Initial envelope (Sigma) matrix.
  FMatrix<double,6,6> getSigma() const;

  /// Envelope (Sigma) matrix for given row.
  FMatrix<double,6,6> getSigma(const Row &) const;


  /// Arc length for given row.
  double getS(const Row &, int=0, int=0) const;

  /// Three modes, "naive" Twiss functions.
  //  Index (0 ... 2) is mode.
  double getMUi(const Row &, int i1, int=0) const;
  double getBETi(const Row &, int i1, int=0) const;
  double getALFi(const Row &, int i1, int=0) const;

  /// Mais-Ripken beta functions.
  //  First index (0 ... 2) is plane, second index (0 ... 2) is mode.
  double getBETik(const Row &, int i1, int i2) const;

  /// Mais-Ripken alpha functions.
  //  First index (0 ... 2) is plane, second index (0 ... 2) is mode.
  double getALFik(const Row &, int i1, int i2) const;

  /// Mais-Ripken gamma functions.
  //  First index (0 ... 2) is plane, second index (0 ... 2) is mode.
  double getGAMik(const Row &, int i1, int i2) const;

  /// Closed orbit.
  //  Index (0 ... 5) is plane.
  double getCO(const Row &, int i1, int=0) const;

  /// Dispersion.
  //  Index (0 ... 5) is plane.
  double getDisp(const Row &, int i1, int=0) const;

  /// Eigenvectors.
  //  First index (0 ... 5) is plane, second index (0 ... 5) is column.
  double getEigen(const Row &, int i1, int i2) const;

  /// Sigma matrix.
  //  Both indices (0 ... 5) refer to planes.
  double getSigma(const Row &, int i1, int i2) const;

  /// Transfer matrix.
  //  First index (0 ... 5) is row, second index (0 ... 5) is column.
  double getMatrix(const Row &, int i1, int i2) const;

protected:

  /// The common attributes for all objects having the "TWISS" interface.
  //  Must be accessible to classes Insertion and Period.
  enum {
    // User-definable attributes:
    LINE,        // The beam line for the table,
    BEAM,        // The beam to be used.
    RANGE,       // The range in the lattice.
    ORDER,       // The order for the calculation.
    STATIC,      // The flag for suppressing recalculation.
    METHOD,      // the algorithm for filling
    REVBEAM,     // If true, beam runs backwrads.
    REVTRACK,    // If true, track lattice functions against the beam.

    // Read-only attributes:
    BETXMAX,     // Maximum beta functions
    BETYMAX,
    XCMAX,       // Maximum closed orbit excursion
    YCMAX,
    XCRMS,       // R.M.S. closed orbit excursion
    YCRMS,
    DXMAX,       // Maximum dispersion
    DYMAX,
    DXRMS,       // R.M.S. dispersion
    DYRMS,
    SIZE
  };

  /// Exemplar constructor.
  Twiss(int size, const char *name, const char *help);

  /// Clone constructor.
  Twiss(const string &name, Twiss *parent);


  /// Number of table columns.
  static const int numColumns = 11;

  /// The initial closed orbit.
  FVector<double,6> orbit;

  /// The initial curly A matrix.
  FMatrix<double,6,6> curly_A;

private:

  // Not implemented.
  Twiss(const Twiss &);
  void operator=(const Twiss &);

  // Return the table row.
  Row &findRow(const PlaceRep &row);

  // This method is called to fill the table after initialising.
  void put();

  // Access to current table row.
  mutable TLine::const_iterator current;

  // The table contents.
  TLine *itsTable;

  // Pointer to the filling algorithm.
  AbstractMapper *itsMapper;

  // The attached beam.
  const Beam *beam;

  // The particle reference data.
  const PartData *reference;

  // The truncation order (1 for LINEAR, 2 for all other).
  int order;

  // The direction flags.
  bool revBeam;   // true, if beam runs from right (s=C) to left (s=0).
  bool revTrack;  // true, if tracking against the beam.
  bool revPath;   // true, if tracking from right (s=C) to left (s=0).

  // The name of the analysed line.
  string itsLine;
};

#endif // OPAL_Twiss_HH
