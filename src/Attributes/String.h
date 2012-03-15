#ifndef OPAL_String_HH
#define OPAL_String_HH

// ------------------------------------------------------------------------
// $RCSfile: String.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: String
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"


// Class String
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type string.
    class String: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        String(const string &name, const string &help);

        virtual ~String();

        /// Read the attribute from the DOOM data base.
        virtual void doomGet(Attribute &, const DoomReader &, int) const;

        /// Write the attribute to the DOOM data base.
        virtual void doomPut(const Attribute &, DoomWriter &, int) const;

        /// Return attribute type string ``string''.
        virtual const string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        String();
        String(const String &);
        void operator=(const String &);
    };

};

#endif // OPAL_String_HH
