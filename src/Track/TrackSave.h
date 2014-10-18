#ifndef OPAL_TrackSave_HH
#define OPAL_TrackSave_HH

// ------------------------------------------------------------------------
// $RCSfile: TrackSave.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackSave
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:47 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class TrackSave
// ------------------------------------------------------------------------
/// The TSAVE command.

class TrackSave: public Action {

public:

    /// Exemplar constructor.
    TrackSave();

    virtual ~TrackSave();

    /// Make clone.
    virtual TrackSave *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    TrackSave(const TrackSave &);
    void operator=(const TrackSave &);

    // Clone constructor.
    TrackSave(const std::string &name, TrackSave *parent);
};

#endif // OPAL_TrackSave_HH
