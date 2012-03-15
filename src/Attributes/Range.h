#ifndef OPAL_Range_HH
#define OPAL_Range_HH

// ------------------------------------------------------------------------
// $RCSfile: Range.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Range
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/RangeRep.h"

class Attribute;
class DoomReader;
class DoomWriter;


// Class Range
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type range definition.
    class Range: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        Range(const string &name, const string &help);

        virtual ~Range();

        /// Return attribute type ``range''.
        virtual const string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

        /// Read the attribute from the DOOM data base.
        virtual void doomGet(Attribute &, const DoomReader &, int) const;

        /// Write the attribute to the DOOM data base.
        virtual void doomPut(const Attribute &, DoomWriter &, int) const;

    private:

        // Not implemented.
        Range();
        Range(const Range &);
        void operator=(const Range &);
    };

};

#endif // OPAL_Range_HH
