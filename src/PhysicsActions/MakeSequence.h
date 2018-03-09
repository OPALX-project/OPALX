#ifndef OPAL_MakeSequence_HH
#define OPAL_MakeSequence_HH

// ------------------------------------------------------------------------
// $RCSfile: MakeSequence.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MakeSequence
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:22:04 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"


// Class MakeSequence
// ------------------------------------------------------------------------
/// The MAKESEQ command.

class MakeSequence: public Action {

public:

    /// Exemplar constructor.
    MakeSequence();

    virtual ~MakeSequence();

    /// Make clone.
    virtual MakeSequence *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    MakeSequence(const MakeSequence &);
    void operator=(const MakeSequence &);

    // Clone constructor.
    MakeSequence(const std::string &name, MakeSequence *parent);
};

#endif // __MakeSequence_HH
