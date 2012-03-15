// ------------------------------------------------------------------------
// $RCSfile: AColumn.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AColumn
//   Expression class. Generates a list of values from a TABLE() function.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Expressions/AColumn.h"
#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/Table.h"
#include "Expressions/SCell.h"
#include <iostream>
#include <vector>


// Class AColumn
// ------------------------------------------------------------------------

namespace Expressions {
    
  AColumn::AColumn(const AColumn &rhs):
    Array<double>(rhs), tab_name(rhs.tab_name),
    col_name(rhs.col_name), itsRange(rhs.itsRange)
    {}

    
  AColumn::AColumn
    (const string &tName, const string &cName, const RangeRep &rng):
    Array<double>(), tab_name(tName), col_name(cName), itsRange(rng)
    {}


  AColumn::~AColumn()
    {}


  Array<double> *AColumn::clone() const
    {
      return new AColumn(*this);
    }


  std::vector<double> AColumn::evaluate() const
    {
      Table *table = Table::find(tab_name);
      table->fill();
      return table->getColumn(itsRange, col_name);
    }


  std::ostream &AColumn::print(std::ostream &os, int) const
    {
      os << "COLUMN(" << tab_name << ',' << col_name << ','
	 << itsRange << ')';
      return os;
    }

}
