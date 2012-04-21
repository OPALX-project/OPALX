#ifndef OPAL_Envelope_HH
#define OPAL_Envelope_HH

// ------------------------------------------------------------------------
// $RCSfile: Envelope.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Envelope
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:45 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "AbstractObjects/Attribute.h"
#include <iosfwd>

class Twiss;


// Class Envelope
// ------------------------------------------------------------------------
/// The ENVELOPE command.

class Envelope: public Action {

public:

    /// Exemplar constructor.
    Envelope();

    virtual ~Envelope();

    /// Make clone.
    virtual Envelope *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Envelope(const Envelope &);
    void operator=(const Envelope &);

    // Clone constructor.
    Envelope(const string &name, Envelope *parent);

    // Do the listing.
    void format(std::ostream &, const Twiss *);

    /// Print Twiss table in envelope representation.
    void formatPrint(std::ostream &, const Twiss *) const;

};

#endif // OPAL_Envelope_HH
