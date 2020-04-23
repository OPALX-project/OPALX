#ifndef OPAL_UString_HH
#define OPAL_UString_HH

// ------------------------------------------------------------------------
// $RCSfile: UString.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: UString
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"


// Class UString
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type string.
    class UString: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        UString(const std::string &name, const std::string &help);

        virtual ~UString();

        /// Return attribute type string ``string''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        UString();
        UString(const UString &);
        void operator=(const UString &);
    };

};

#endif // OPAL_UString_HH