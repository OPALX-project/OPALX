#ifndef OPAL_TrackEnd_HH
#define OPAL_TrackEnd_HH

// ------------------------------------------------------------------------
// $RCSfile: TrackEnd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackEnd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"

class Sequence;

// Class TrackEnd
// ------------------------------------------------------------------------
/// The ENDTRACK command.

class TrackEnd : public Action {
public:
    /// Exemplar constructor.
    TrackEnd();

    virtual ~TrackEnd();

    /// Execute the command.
    virtual void execute();

private:
    // Not implemented.
    TrackEnd(const TrackEnd&);
    void operator=(const TrackEnd&);
};

#endif  // OPAL_TrackEnd_HH
