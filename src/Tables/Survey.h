#ifndef OPAL_Survey_HH
#define OPAL_Survey_HH

// ------------------------------------------------------------------------
// $RCSfile: Survey.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Survey
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:22 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Table.h"
#include "AbstractObjects/Expressions.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Beamlines/TBeamline.h"
#include <iosfwd>
#include <vector>

class DoomWriter;
class PlaceRep;
class RangeRep;
class Surveyor;


// Class Survey
// ------------------------------------------------------------------------
/// The SURVEY command.

class Survey: public Table {

public:

  /// The class for one row of the survey table.
  class Row: public FlaggedElmPtr {

    friend class ::Survey;

  public:

    Row(ElementBase *, int);
    explicit Row(const FlaggedElmPtr &);
    ~Row();

    /// Return the accumulated geometry transform.
    const Euclid3D &getMap() const;

    /// Return the accumulated length.
    double getS() const;

  private:

    /// The accumulated geometry transform.
    Euclid3D euclid;

    /// The accumulated length.
    double s;
  };

  /// Exemplar constructor.
  Survey();

  virtual ~Survey();

  /// Make clone.
  virtual Survey *clone(const string &name);

  /// Store the table in the DOOM data base.
  virtual void doomPut(DoomWriter &) const;

  /// Check validity of survey definition.
  virtual void execute();

  /// Fill the buffer using the survey algorithm.
  virtual void fill();

  /// Return a selected value in a selected row.
  virtual double getCell(const PlaceRep &row, const string &col);

  /// Return the default print columns.
  virtual CellArray getDefault() const;

  /// Return column [b]col[/b] of this table, limited by [b]range[/b].
  virtual std::vector<double>
  getColumn(const RangeRep &range, const string &col);

  /// Return current row of table.
  const Row &getCurrent() const;


  /// Return the length of the table.
  virtual double getLength();

  /// Return embedded CLASSIC beamline.
  virtual const Beamline *getLine() const;

  /// Return a table row, possible user-defined.
  virtual std::vector<double>
  getRow(const PlaceRep &, const std::vector<string> &);


  /// Arc length for given row.
  double getS(const Row &, int=0, int=0) const;

  /// Position and orientation of local system.
  const Euclid3D &getMap(const Row &) const;

  /// X component of displacement.
  double getX(const Row &, int=0, int=0) const;

  /// Y component of displacement.
  double getY(const Row &, int=0, int=0) const;

  /// Z component of displacement.
  double getZ(const Row &, int=0, int=0) const;

  /// Rotation about X.
  double getPhi(const Row &, int=0, int=0) const;

  /// Rotation about Y.
  double getTheta(const Row &, int=0, int=0) const;

  /// Rotation about Z.
  double getPsi(const Row &, int=0, int=0) const;

  /// Local axis vectors.
  //  First index (1 ... 3) is coordinate, second index (1 ... 3) is vector.
  double getW(const Row &, int i1, int i2) const;


  /// Find dependency.
  //  Return true, if this table depends on the named object.
  virtual bool isDependent(const string &name) const;

  /// Return column.
  //  Return an expression which denotes the selected column,
  //  identified by its name.
  virtual Expressions::PtrToScalar<double>
  makeColumnExpression(const string &colName) const;

  /// Write TFS file for this table.
  virtual void makeTFS(std::ostream &, const CellArray &) const;

  /// Check compatibility.
  //  True, if [b]rhs[/b] is a survey table.
  virtual bool matches(Table *rhs) const;

  /// Print list for the table.
  virtual void printTable(std::ostream &, const CellArray &) const;

private:

  // Not implemented.
  Survey(const Survey &);
  void operator=(const Survey &);

  // Clone constructor.
  Survey(const string &name, Survey *parent);


  // Set the current table row for a given place specification.
  const Row &findRow(const PlaceRep &row);

  // The contained beamline type.
  typedef TBeamline<Row> TLine;

  // Access to current table row.
  mutable TLine::const_iterator current;

  // The table contents.
  TLine *itsTable;

  // The algorithm for filling the buffer.
  Surveyor *itsVisitor;

  // Current accumulated design length.
  double s;

  // The name of the surveyed line.
  string itsLine;

  /// Number of table columns.
  static const int numColumns = 7;
};

#endif // OPAL_Survey_HH
