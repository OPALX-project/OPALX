#ifndef OPAL_Reference_HH
#define OPAL_Reference_HH

// ------------------------------------------------------------------------
// $RCSfile: Reference.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Reference:
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeHandler.h"

class Attribute;


// Class Reference
// ------------------------------------------------------------------------

namespace Attributes {

  /// Parser for an attribute of type attribute reference.
  //  The attribute referred to may be logical, real, or string.
  class Reference: public AttributeHandler {

  public:

    /// Constructor.
    //  Assign attribute name and help string.
    Reference(const string &name, const string &help);

    virtual ~Reference();

    /// Read the attribute from the DOOM data base.
    virtual void doomGet(Attribute &, const DoomReader &, int) const;

    /// Write the attribute to the DOOM data base.
    virtual void doomPut(const Attribute &, DoomWriter &, int) const;

    /// Return attribute type string ``reference''.
    virtual const string &getType() const;

    /// Parse the attribute.
    virtual void parse(Attribute &, Statement &, bool) const;

  private:

    // Not implemented.
    Reference(const Reference &);
    void operator=(const Reference &);
  };

};

#endif // OPAL_Reference_HH
