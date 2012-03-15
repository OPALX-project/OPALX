#ifndef OPAL_TrackStart_HH
#define OPAL_TrackStart_HH

// ------------------------------------------------------------------------
// $RCSfile: TrackStart.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackStart
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:47 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class TrackStart
// ------------------------------------------------------------------------
/// The START command.

class TrackStart: public Action {

public:

    /// Exemplar constructor.
    TrackStart();

    virtual ~TrackStart();

    /// Make clone.
    virtual TrackStart *clone(const string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    TrackStart(const TrackStart &);
    void operator=(const TrackStart &);

    // Clone constructor.
    TrackStart(const string &name, TrackStart *parent);
};

#endif // OPAL_TrackStart_HH
