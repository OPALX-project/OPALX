#ifndef OPAL_OpalInstrument_HH
#define OPAL_OpalInstrument_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalInstrument.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalInstrument
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"


// Class OpalInstrument
// ------------------------------------------------------------------------
/// The INSTRUMENT element.

class OpalInstrument: public OpalElement {

public:

    /// Exemplar constructor.
    OpalInstrument();

    virtual ~OpalInstrument();

    /// Make clone.
    virtual OpalInstrument *clone(const std::string &name);

    /// Update the embedded CLASSIC drift.
    virtual void update();

private:

    // Not implemented.
    OpalInstrument(const OpalInstrument &);
    void operator=(const OpalInstrument &);

    // Clone constructor.
    OpalInstrument(const std::string &name, OpalInstrument *parent);
};

#endif // OPAL_OpalInstrument_HH
