#ifndef OPAL_TableRow_HH
#define OPAL_TableRow_HH

// ------------------------------------------------------------------------
// $RCSfile: TableRow.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TableRow
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/TableRowRep.h"

class TableRowRep;
class DoomReader;
class DoomWriter;


// Class TableRow
// ------------------------------------------------------------------------

namespace Attributes {

  /// Parser for an attribute of type table row reference.
  class TableRow: public AttributeHandler {

  public:

    /// Constructor.
    //  Assign attribute name and help string.
    TableRow(const string &name, const string &help);

    virtual ~TableRow();

    /// Read the attribute from the DOOM data base.
    virtual void doomGet(Attribute &, const DoomReader &, int) const;

    /// Write the attribute to the DOOM data base.
    virtual void doomPut(const Attribute &, DoomWriter &, int) const;

    /// Return attribute type string ``table line''.
    virtual const string &getType() const;

    /// Parse the attribute.
    virtual void parse(Attribute &, Statement &, bool) const;

  private:

    // Not implemented.
    TableRow(const TableRow &);
    void operator=(const TableRow &);
  };

};

#endif // OPAL_TableRow_HH
