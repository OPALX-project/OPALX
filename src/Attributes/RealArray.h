#ifndef OPAL_RealArray_HH
#define OPAL_RealArray_HH

// ------------------------------------------------------------------------
// $RCSfile: RealArray.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class RealArray:
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"


// Class RealArray
// ------------------------------------------------------------------------

namespace Attributes {

  /// Parser for an attribute of type real array.
  class RealArray: public AttributeHandler {

  public:

    /// Constructor.
    //  Assign attribute name and help string.
    RealArray(const string &name, const string &help);

    virtual ~RealArray();

    /// Read the attribute from the DOOM data base.
    virtual void doomGet(Attribute &, const DoomReader &, int) const;

    /// Write the attribute to the DOOM data base.
    virtual void doomPut(const Attribute &, DoomWriter &, int) const;

    /// Return attribute type string ``real array''.
    virtual const string &getType() const;

    /// Parse the attribute.
    virtual void parse(Attribute &, Statement &, bool) const;

    /// Parse a component of the array.
    //  Identified by its index.
    virtual void parseComponent(Attribute &, Statement &, bool, int) const;

  private:

    // Not implemented.
    RealArray();
    RealArray(const RealArray &);
    void operator=(const RealArray &);
  };

};

#endif // OPAL_RealArray_HH
