// ------------------------------------------------------------------------
// $RCSfile: TableRow.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TableRow
//   A class used to parse a reference to a table lines.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/TableRow.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"

using namespace Expressions;


// Class TableRow
// ------------------------------------------------------------------------

namespace Attributes {

  TableRow::TableRow(const string &name, const string &help):
    AttributeHandler(name, help, 0)
  {}


  TableRow::~TableRow()
  {}
 

  void TableRow::doomGet
  (Attribute &attr, const DoomReader &reader, int index) const
  {
    // MISSING.
  }


  void TableRow::doomPut
  (const Attribute &attr, DoomWriter &writer, int index) const
  {
    // MISSING.
  }


  const string &TableRow::getType() const
  {
    static const string type("table line");
    return type;
  }


  void TableRow::parse(Attribute &attr, Statement &stat, bool) const
  {
    attr.set(new SValue<TableRowRep>(parseTableRow(stat)));
  }

};
