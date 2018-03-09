#ifndef OPAL_AlignRemover_HH
#define OPAL_AlignRemover_HH

// ------------------------------------------------------------------------
// $RCSfile: AlignRemover.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignRemover
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/AlignBase.h"
#include <fstream>

class AlignWrapper;


// Class AlignRemover
// ------------------------------------------------------------------------
/// A visitor used to remove all misalignments from a line.
//  Uses an AlignHandler to access all elements which have a misalignment,
//  and sets all misalignments to zero.

class AlignRemover: public AlignBase {

public:

    AlignRemover();
    virtual ~AlignRemover();

    /// Remove misalignement.
    virtual void misalignment(const AlignWrapper &, int occur);

private:

    // Not implemented.
    AlignRemover(const AlignRemover &);
    void operator=(const AlignRemover &);
};

#endif // OPAL_AlignRemover_HH
