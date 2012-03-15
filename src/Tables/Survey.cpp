// ------------------------------------------------------------------------
// $RCSfile: Survey.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Survey
//   The Survey class implements a buffer for survey tables.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/Survey.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/RangeRep.h"
#include "AbstractObjects/TableRowRep.h"
#include "Attributes/Attributes.h"
#include "Algorithms/Surveyor.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Physics/Physics.h"
#include "Tables/Flatten.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "doomex.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#if defined(__GNUC__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif


// Local structures.
// ------------------------------------------------------------------------

namespace {

  // Describe a column with name and function to get it.
  struct ColDesc {

    // The column name.
    const char *colName;

    // Method to return the selected column from the given row.
    double (Survey::*get)(const Survey::Row &, int, int) const;

    int printWidth;
    int printPrecision;

    int ind_1;
    int ind_2;
  };

  // The column entries table.
  const ColDesc allColumns[] = {
    { "S",     &Survey::getS,     14,  6, 0, 0 },

    { "X",     &Survey::getX,     10,  6, 0, 0 },
    { "Y",     &Survey::getY,     10,  6, 0, 0 },
    { "Z",     &Survey::getZ,     10,  6, 0, 0 },

    { "THETA", &Survey::getTheta, 14, 10, 0, 0 },
    { "PHI",   &Survey::getPhi,   14, 10, 0, 0 },
    { "PSI",   &Survey::getPsi,   14, 10, 0, 0 },

    { "W11",   &Survey::getW,     10,  6, 1, 1 },
    { "W21",   &Survey::getW,     10,  6, 2, 1 },
    { "W31",   &Survey::getW,     10,  6, 3, 1 },
    { "W12",   &Survey::getW,     10,  6, 1, 2 },
    { "W22",   &Survey::getW,     10,  6, 2, 2 },
    { "W32",   &Survey::getW,     10,  6, 3, 2 },
    { "W13",   &Survey::getW,     10,  6, 1, 3 },
    { "W23",   &Survey::getW,     10,  6, 2, 3 },
    { "W33",   &Survey::getW,     10,  6, 3, 3 },

    { 0,       0,                  0,  0, 0, 0 }
  };


  // The print column entries table.
  const ColDesc printColumns[] = {
    { "S",     &Survey::getS,     14,  6, 0, 0 },

    { "X",     &Survey::getX,     10,  6, 0, 0 },
    { "Y",     &Survey::getY,     10,  6, 0, 0 },
    { "Z",     &Survey::getZ,     10,  6, 0, 0 },

    { "THETA", &Survey::getTheta, 14, 10, 0, 0 },
    { "PHI",   &Survey::getPhi,   14, 10, 0, 0 },
    { "PSI",   &Survey::getPsi,   14, 10, 0, 0 },

    { 0,       0,                  0,  0, 0, 0 }
  };


  // The TFS column entries table.
  const ColDesc listTFScolumns[] = {
    { "S",   &Survey::getS, 14, 6, 0, 0 },

    { "X",   &Survey::getX, 10, 6, 0, 0 },
    { "Y",   &Survey::getY, 10, 6, 0, 0 },
    { "Z",   &Survey::getZ, 10, 6, 0, 0 },

    { "W11", &Survey::getW, 10, 6, 0, 0 },
    { "W21", &Survey::getW, 10, 6, 1, 0 },
    { "W31", &Survey::getW, 10, 6, 2, 0 },
    { "W12", &Survey::getW, 10, 6, 0, 1 },
    { "W22", &Survey::getW, 10, 6, 1, 1 },
    { "W32", &Survey::getW, 10, 6, 2, 1 },
    { "W13", &Survey::getW, 10, 6, 0, 2 },
    { "W23", &Survey::getW, 10, 6, 1, 2 },
    { "W33", &Survey::getW, 10, 6, 2, 2 },

    { 0,       0,            0, 0, 0, 0 }
  };

  // Find a column by name.
  const ColDesc *findCol(const Survey &table, const string &colName)
  {
    for (const ColDesc *col = allColumns; col->colName; ++col) {
      if (colName == col->colName) {
	return col;
      }
    }

    throw OpalException("Survey::findCol()",
		       "Survey table \"" + table.getOpalName() +
		       "\" has no column named \"" + colName + "\".");
  }


  // Local class Column.
  //   Returns the value for a given column in the current row of a given
  //   table.
  class Column: public Expressions::Scalar<double> {

  public:

    //: Constructor.
    //  Identify the table by its name [b]tab[/b], and the column by its
    //  name [b]col[/b] and the function [b]col[/b].
    //  The row is specified as the ``current'' row of the table.
    Column(const Survey &tab, const string &scol, const ColDesc &col);

    Column(const Column &);
    virtual ~Column();

    //: Make clone.
    virtual Expressions::Scalar<double> *clone() const;

    //: Evaluate.
    virtual double evaluate() const;

    //: Print expression.
    virtual std::ostream &print(std::ostream &os, int precedence=99) const;

  private:

    // Not implemented.
    Column();
    const Column &operator=(const Column &);

    // The Table referred.
    const Survey &itsTable;

    // Column name.
    string colName;

    // The function returning the column value from the given row.
    double (Survey::*get)(const Survey::Row &, int, int) const;

    // Indices for calling get().
    int ind_1, ind_2;
  };


  // Implementation.
  // ------------------------------------------------------------------------

  Column::Column(const Survey &tab,
		 const string &name,
		 const ColDesc &desc):
    itsTable(tab),
    colName(name),
    get(desc.get),
    ind_1(desc.ind_1),
    ind_2(desc.ind_2)
  {}


  Column::Column(const Column &rhs):
    itsTable(rhs.itsTable),
    colName(rhs.colName),
    get(rhs.get),
    ind_1(rhs.ind_1),
    ind_2(rhs.ind_2)
  {}


  Column::~Column()
  {}


  Expressions::Scalar<double> *Column::clone() const
  {
    return new Column(*this);
  }


  double Column::evaluate() const
  {
    // Apply the function to the current row of the table.
    return (itsTable.*get)(itsTable.getCurrent(), ind_1, ind_2);
  }


  std::ostream &Column::print(std::ostream &os, int) const
  {
    return os << colName;
  }
};


// Class Survey::Row
// ------------------------------------------------------------------------

Survey::Row::Row(ElementBase *elem, int occur):
  FlaggedElmPtr(elem)
{
  setCounter(occur);
}


Survey::Row::Row(const FlaggedElmPtr &rhs):
  FlaggedElmPtr(rhs)
{}


Survey::Row::~Row()
{}


double Survey::Row::getS() const
{
  return s;
}


const Euclid3D &Survey::Row::getMap() const
{
  return euclid;
}


// Class Survey
// ------------------------------------------------------------------------

// The attributes of class Survey.
namespace {
  enum {
    LINE,        // The lattice to be used.
    RANGE,       // The range in the lattice.
    STATIC,      // The flag for suppressing recalculation.
    INIT,        // If given, take initial values from another table.
    X0,          // The initial position  X.
    Y0,          // The initial position  Y.
    Z0,          // The initial position  Z.
    THETA0,      // The initial azimuth   THETA.
    PHI0,        // The initial elevation PHI.
    PSI0,        // The initial roll      PSI.
    L,           // Sum of design lengths.
    SIZE
  };
}


Survey::Survey():
  Table(SIZE, "SURVEY",
	"The \"SURVEY\" command defines a table of survey data which "
	"can be matched or tabulated."),
  itsTable(0), itsVisitor(0)
{
  itsAttr[LINE] = Attributes::makeString
    ("LINE", "Name of the lattice to be surveyed");
  itsAttr[RANGE] = Attributes::makeRange
    ("RANGE", "The range in the lattice");
  itsAttr[STATIC] = Attributes::makeBool
    ("STATIC", "If true, the table is not recalculated at each iteration");
  itsAttr[INIT] = Attributes::makeTableRow
    ("INIT", "If given, this table position is used to initialise");
  itsAttr[X0] = Attributes::makeReal
    ("X0", "Initial X position in m");
  itsAttr[Y0] = Attributes::makeReal
    ("Y0", "Initial Y position in m");
  itsAttr[Z0] = Attributes::makeReal
    ("Z0", "Initial Z position in m");
  itsAttr[THETA0] = Attributes::makeReal
    ("THETA0", "Initial azimuth angle in rad");
  itsAttr[PHI0] = Attributes::makeReal
    ("PHI0", "Initial pitch angle in rad");
  itsAttr[PSI0] = Attributes::makeReal
    ("PSI0", "Initial roll angle in rad");

  // READ ONLY.
  itsAttr[L] = Attributes::makeReal
    ("L", "Sum of design lengths in m");
  itsAttr[L].setReadOnly(true);
}


Survey::Survey(const string &name, Survey *parent):
  Table(name, parent), itsTable(new TLine(name)), itsVisitor(0)
{}


Survey::~Survey()
{
  if (itsTable)   delete itsTable;
  if (itsVisitor) delete itsVisitor;
}


Survey *Survey::clone(const string &name)
{
  return new Survey(name, this);
}


void Survey::doomPut(DoomWriter &writer) const
{
  // Write the definition for this table.
  Object::doomPut(writer);

  // Make sure the table is up to date.
  // Cast away const, to allow logically constant table to update.
  const_cast<Survey *>(this)->fill();

  // Make the table header.
  const string &tableName = getOpalName();
  {
    static const int headKeyList[] = { 1, TABLE_HEAD };
    DoomWriter headWriter(tableName, headKeyList);
    headWriter.setParentName(itsLine);
    headWriter.setTypeName("TABLE_HEADER");

    static const char *columnName[numColumns] = {
      "PHI", "PSI", "S", "THETA", "X", "Y", "Z"
	};
    static const int columnNumber[numColumns] = {
      5, 6, 0, 4, 1, 2, 3
	};
    for (int i = 0; i < numColumns; ++i) {
      headWriter.putInt(i, columnNumber[i]);
      headWriter.putString(i, columnName[i]);
    }

    int numRows = itsTable->size();
    headWriter.putInt(numColumns + 0, 1);          // 1 delta value (ignored).
    headWriter.putInt(numColumns + 1, numRows);    // number of rows.
    headWriter.putInt(numColumns + 3, 3);          // position at end.
    // Header object will be written when headWriter goes out of scope.
  }

  // Write the table data.
  {
    static const int rowKeyList[] = { 2, TABLE_BODY, 0 };
    DoomWriter rowWriter(tableName, rowKeyList);
    rowWriter.setParentName(itsLine);
    rowWriter.setTypeName("SUB_TABLE");
    int lineCount = 0;
    int realCount = 0;
    for (current = itsTable->begin(); current != itsTable->end(); ++current) {
      rowWriter.putInt(lineCount,      current->getCounter());
      rowWriter.putString(lineCount,   current->getElement()->getName());
      rowWriter.putReal(realCount + 0, getS(*current));
      rowWriter.putReal(realCount + 1, getX(*current));
      rowWriter.putReal(realCount + 2, getY(*current));
      rowWriter.putReal(realCount + 3, getZ(*current));
      rowWriter.putReal(realCount + 4, getPhi(*current));
      rowWriter.putReal(realCount + 5, getTheta(*current));
      rowWriter.putReal(realCount + 6, getPsi(*current));
      lineCount++;
      realCount += numColumns;
    }
    // Body object will be written when Rowwriter goes out of scope.
  }
}


void Survey::execute()
{
  // Find BeamSequence definition.
  itsLine = Attributes::getString(itsAttr[LINE]);
  BeamSequence *use = BeamSequence::find(itsLine);

  // Make sure all is up-to-date.
  OPAL.update();

  // Create flat list with space for data storage.
  RangeRep range = Attributes::getRange(itsAttr[RANGE]);
  Flatten<Row> flattener(*use->fetchLine(), *itsTable, range);
  flattener.execute();

  if (itsTable->empty() && Options::warn) {
    std::cerr << "\n### Warning ### Survey table \"" << getOpalName()
	      << "\" contains no elements.\n" << std::endl;
  } else {
    itsTable->front().setSelectionFlag(true);
    itsTable->back().setSelectionFlag(true);
  }

  // Fill the table.
  itsVisitor = new Surveyor(*itsTable, false);
  fill();

  // Set static flag.
  if (Attributes::getBool(itsAttr[STATIC])) dynamic = false;
  printTable(std::cout, getDefault());
}


void Survey::fill()
{
  if (refill) {
    // Reset the surveyor.
    double x, y, z, phi, theta, psi;
    if (itsAttr[INIT]) {
      TableRowRep rowrep = Attributes::getTableRow(itsAttr[INIT]);
      Table *init = Table::find(rowrep.getTabName());

      if (matches(init)) {
	PlaceRep pinit = rowrep.getPosition();
	Survey  *sinit = dynamic_cast<Survey *>(init);
	sinit->fill();
	const Row &row = sinit->findRow(pinit);
	itsVisitor->setMap(sinit->getMap(row));
      } else {
	throw OpalException("Survey::fill()", "Table \"" + rowrep.getTabName() +
			   "\" is not suitable for initialising survey \"" +
			   getOpalName() + "\".");
      }
    } else {
      x     = Attributes::getReal(itsAttr[X0]);
      y     = Attributes::getReal(itsAttr[Y0]);
      z     = Attributes::getReal(itsAttr[Z0]);
      phi   = Attributes::getReal(itsAttr[PHI0]);
      theta = Attributes::getReal(itsAttr[THETA0]);
      psi   = Attributes::getReal(itsAttr[PSI0]);

      // Use the OPAL-8 conventions.
      Rotation3D rot = Rotation3D::YRotation(theta) *
	Rotation3D::XRotation(- phi) * Rotation3D::ZRotation(psi);
      itsVisitor->setMap(Euclid3D(Vector3D(x, y, z), rot));
    }

    s = 0.0;

    for (TLine::iterator i = itsTable->begin(); i != itsTable->end(); ++i) {
      // Advance through element.
      i->accept(*itsVisitor);

      // Update accumulated length.
      ElementBase &elem = *i->getElement();
      if (! dynamic_cast<Beamline *>(&elem)) {
	s += elem.getElementLength();
      }

      // Copy row to table line.
      itsVisitor->getMap(i->euclid);
      i->s = s;
    }

    Euclid3D map;
    itsVisitor->getMap(map);
    map.getAll(x, y, z, phi, theta, psi);
    Attributes::setReal(itsAttr[L], s);
    refill = false;
  }
}


const Survey::Row &Survey::findRow(const PlaceRep &place)
{
  PlaceRep row(place);
  row.initialize();
  for (TLine::const_iterator i = itsTable->begin();
       i != itsTable->end(); ++i) {
    row.enter(*i);
    if (row.isActive()) return *i;
    row.leave(*i);
  }

#if defined(__GNUC__) && __GNUC__ < 3
  char buffer[128];
  std::ostrstream os(buffer, sizeof(buffer));
#else
  std::ostringstream os;
#endif
  os << row << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
  string name(buffer);
  throw OpalException("Survey::setRow()", "Row \"" + name +
		     "\" not found in survey table \"" + getOpalName() + "\".");
#else
  throw OpalException("Survey::setRow()", "Row \"" + os.str() +
		     "\" not found in survey table \"" + getOpalName() + "\".");
#endif
}


double Survey::getCell(const PlaceRep &place, const string &name)
{
  const Row &row = findRow(place);
  const ColDesc *col = findCol(*this, name);
  return (this->*(col->get))(row, col->ind_1, col->ind_2);
}


Table::CellArray Survey::getDefault() const
{
  CellArray columns;
  for (const ColDesc *col = printColumns; col->colName; ++col) {
    Expressions::PtrToScalar<double> expr =
      new Column(*this, col->colName, *col);
    columns.push_back(Cell(expr, col->printWidth, col->printPrecision));
  }
  return columns;
}


std::vector<double> Survey::getColumn(const RangeRep &rng, const string &name)
{
  // Find proper column function.
  const ColDesc *col = findCol(*this, name);
  RangeRep range(rng);
  range.initialize();
  std::vector<double> column;

  for (current = itsTable->begin(); current != itsTable->end(); ++current) {
    range.enter(*current);
    if (range.isActive()) {
      column.push_back((this->*(col->get))(*current, col->ind_1, col->ind_2));
    }
    range.leave(*current);
  }

  return column;
}


const Survey::Row &Survey::getCurrent() const
{
  return *current;
}


double Survey::getLength()
{
  return itsTable->getElementLength();
}


const Beamline *Survey::getLine() const
{
  return itsTable;
}


std::vector<double>
Survey::getRow(const PlaceRep &pos, const std::vector<string> &columns)
{
  const Row &row = findRow(pos);
  std::vector<double> result;

  if (columns.empty()) {
    // Standard column selection.
    for (const ColDesc *col = printColumns; col->colName != 0; ++col) {
      result.push_back((this->*(col->get))(row, col->ind_1, col->ind_2));
    }
  } else {
    // User column selection.
    for (std::vector<string>::const_iterator name = columns.begin();
	 name != columns.end(); ++name) {
      const ColDesc *col = findCol(*this, *name);
      result.push_back((this->*(col->get))(row, col->ind_1, col->ind_2));
    }
  }

  return result;
}


bool Survey::isDependent(const string &name) const
{
  // Test if name refers to LINE attribute.
  if (itsLine == name) return true;

  // Test if name occurs in table.
  for (TLine::const_iterator i = itsTable->begin(); i != itsTable->end(); ++i) {
    if (i->getElement()->getName() == name) return true;
  }

  // Otherwise replacement is not required.
  return false;
}


Expressions::PtrToScalar<double>
Survey::makeColumnExpression(const string &colName) const
{
  return new Column(*this, colName, *findCol(*this, colName));
}


void Survey::makeTFS(std::ostream &os, const CellArray &cells) const
{
  // Save the formatting flags.
  std::streamsize old_prec = os.precision(12);
  os.setf(std::ios::fixed, std::ios::floatfield);

  // Write table descriptors.
  os << "@ TYPE     %s  SURVEY\n";
  Table::tfsTableDescriptors(os);
  os << "@ LINE     %s  " << itsAttr[LINE] << '\n'
     << "@ LENGTH   %le " << itsTable->back().s << '\n';

  // Write column headers, names.
  os << "* NAME";
  for (CellArray::const_iterator col = cells.begin();
       col < cells.end(); ++col) {
#if defined(__GNUC__) && __GNUC__ < 3
    char buffer[256];
    std::ostrstream ss(buffer, 256);
#else
    std::ostringstream ss;
#endif
    col->itsExpr->print(ss, 0);
    ss << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
    os << "  " << buffer;
    string image(buffer); // This does nothing!
#else
    os << ss.str();
#endif
  }
  os << '\n';

  // Write column headers, formats.
  os << "$ %s";
  for (CellArray::size_type i = 0; i < cells.size(); ++i) {
    os << " %le";
  }
  os << '\n';

  // Write table body.
  for (current = itsTable->begin(); current != itsTable->end(); ++current) {
    if (current->getSelectionFlag()) {
      os << "  " << current->getElement()->getName();

      if (cells.empty()) {
	// Standard column selection.
	for (const ColDesc *col = listTFScolumns; col->colName != 0; ++col) {
	  os << std::setw(col->printWidth)
	     << std::setprecision(col->printPrecision)
	     << ' ' << (this->*(col->get))(*current, col->ind_1, col->ind_2);
	}
      } else {
	for (CellArray::const_iterator cell = cells.begin();
	     cell != cells.end(); ++cell) {
	  os << std::setw(cell->printWidth)
	     << std::setprecision(cell->printPrecision)
	     << ' ' << cell->itsExpr->evaluate();
	}
      }

      os << '\n';
    }
  }

  // Restore the formatting flags.
  os.precision(old_prec);
  os.setf(std::ios::fixed, std::ios::floatfield);
}


bool Survey::matches(Table *rhs) const
{
  return dynamic_cast<Survey *>(rhs) != 0;
}


void Survey::printTable(std::ostream &os, const CellArray &cells) const
{
  // Write table specific header.
  OPAL.printTitle(os);
  os << '\n';
  os << "Geometric layout for line: " << itsAttr[LINE]
     << ", range: " << itsAttr[RANGE] << "." << '\n';

  // Find line length.
  int lineLength = 16;
  for (CellArray::const_iterator cell = cells.begin();
       cell < cells.end(); ++cell) {
    lineLength += cell->printWidth;
  }
  os << string(lineLength, '-') << '\n';

  // Print table header.
  os << "Element         ";
  for (CellArray::const_iterator cell = cells.begin();
       cell < cells.end(); ++cell) {
#if defined(__GNUC__) && __GNUC__ < 3
    char buffer[256];
    std::ostrstream ss(buffer, 256);
#else
    std::ostringstream ss;
#endif
    cell->itsExpr->print(ss, 0);
    ss << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
    string image(buffer);
#else
    std::string image=ss.str();
#endif

    if (int(image.length()) < cell->printWidth) {
      // Right adjust the column header.
      os << string(cell->printWidth - image.length(), ' ') << image;
    } else {
      // Truncate the column header.
      os << ' ' << string(image, 0, cell->printWidth - 3) << "..";
    }
  }
  os << '\n';
  os << string(lineLength, '-') << '\n';

  // Save the formatting flags.
  std::streamsize old_prec = os.precision(6);
  os.setf(std::ios::fixed, std::ios::floatfield);

  // Write table body.
  for (current = itsTable->begin(); current != itsTable->end(); ++current) {
    if (current->getSelectionFlag()) {
      string name = current->getElement()->getName();
      if (int occur = current->getCounter()) {
#if defined(__GNUC__) && __GNUC__ < 3
	char buffer[128];
	std::ostrstream tos(buffer, 128);
#else
	std::ostringstream tos;
#endif
	tos << name << '[' << occur << ']' << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
	name = buffer;
#else
	name = tos.str();
#endif
      }

      if (name.length() > 16) {
	// Truncate the element name.
	os << string(name, 0, 13) << ".. ";
      } else {
	// Left adjust the element name.
	os << name << string(16 - name.length(), ' ');
      }

      for (CellArray::const_iterator cell = cells.begin();
	   cell != cells.end(); ++cell) {
	os << std::setw(cell->printWidth)
	   << std::setprecision(cell->printPrecision)
	   << cell->itsExpr->evaluate();
      }
      os << '\n';
    }
  }

  // Write table specific summary.
  os << string(lineLength, '-') << '\n';

  // Restore the formatting flags.
  os.precision(old_prec);
  os.setf(std::ios::fixed, std::ios::floatfield);
}


double Survey::getS(const Survey::Row &row, int, int) const
{
  return current->getS();
}


const Euclid3D &Survey::getMap(const Survey::Row &row) const
{
  return row.getMap();
}


double Survey::getX(const Survey::Row &row, int, int) const
{
  return row.getMap().getX();
}


double Survey::getY(const Survey::Row &row, int, int) const
{
  return row.getMap().getY();
}


double Survey::getZ(const Survey::Row &row, int, int) const
{
  return row.getMap().getZ();
}


double Survey::getPhi(const Survey::Row &row, int, int) const
{
  const Rotation3D &rot = row.getMap().getRotation();
  return atan2(rot(1,2), sqrt(rot(1,0)*rot(1,0) + rot(1,1)*rot(1,1)));
}


double Survey::getTheta(const Survey::Row &row, int, int) const
{
  const Rotation3D &rot = row.getMap().getRotation();
  double arg = abs(rot(0,2)) + abs(rot(2,2));
  return (arg > 1.0e-10) ? atan2(rot(0,2), rot(2,2)) : 0.0;
}


double Survey::getPsi(const Survey::Row &row, int, int) const
{
  const Rotation3D &rot = row.getMap().getRotation();
  double arg = abs(rot(1,0)) + abs(rot(1,1));
  return (arg > 1.0e-10) ? atan2(rot(1,0), rot(1,1)) : 0.0;
}


double Survey::getW(const Survey::Row &row, int ind_1, int ind_2) const
{
  return row.getMap().M(ind_1, ind_2);
}
