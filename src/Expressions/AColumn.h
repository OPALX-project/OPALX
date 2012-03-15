#ifndef OPAL_AColumn_HH
#define OPAL_AColumn_HH

// ------------------------------------------------------------------------
// $RCSfile: AColumn.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AColumn
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/RangeRep.h"
#include <iosfwd>

class Table;


namespace Expressions {

  // Class AColumn.
  // ----------------------------------------------------------------------
  /// An array expression defined as a table column.
  //  The result is defined as the values contained in that column,
  //  delimited by a range.
  
  class AColumn: public Array<double> {
    
  public:

    /// Constructor.
    //  Use table name [b]tName[/b], column name [b]cName[/b], and range
    //  representation [b]rng[/b].
    AColumn(const string &tName, const string &cName,
	    const RangeRep &rng);
    
    AColumn(const AColumn &);
    ~AColumn();
    
    /// Make clone.
    virtual Array<double> *clone() const;

    /// Evaluate.
    virtual std::vector<double> evaluate() const;
  
    /// Print expression.
    virtual std::ostream &print(std::ostream &os, int precedence=99) const;
  
  private:
  
    // Not implemented.
    AColumn();
    const AColumn &operator=(const AColumn &);
  
    // Table and column name.
    const string tab_name;
    const string col_name;

    // The range specification.
    RangeRep itsRange;
  };

}

#endif // OPAL_AColumn_HH
