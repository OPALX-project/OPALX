#ifndef OPAL_Place_HH
#define OPAL_Place_HH

// ------------------------------------------------------------------------
// $RCSfile: Place.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Place
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/PlaceRep.h"

class Attribute;
class DoomReader;
class DoomWriter;


// Class Place
// ------------------------------------------------------------------------

namespace Attributes {

    /// Parser for an attribute of type place reference.
    class Place: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        Place(const string &name, const string &help);

        virtual ~Place();

        /// Read the attribute from the DOOM data base.
        virtual void doomGet(Attribute &, const DoomReader &, int) const;

        /// Write the attribute to the DOOM data base.
        virtual void doomPut(const Attribute &, DoomWriter &, int) const;

        /// Return attribute type string ``place''.
        virtual const string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        Place();
        Place(const Place &);
        void operator=(const Place &);
    };

};

#endif // OPAL_Place_HH
