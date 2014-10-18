#ifndef OPAL_AlignWriter_HH
#define OPAL_AlignWriter_HH

// ------------------------------------------------------------------------
// $RCSfile: AlignWriter.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignWriter
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/AlignBase.h"


// Class AlignWriter
// ------------------------------------------------------------------------
/// A DOOM writer for writing alignment errors.
//  Retrieves misalignments from a beam line and uses an AlignHandler to
//  store them in the doom data base.


class AlignWriter: public AlignBase {

public:

    /// Constructor.
    //  Store [b]name[/b] in the DOOM environment.
    AlignWriter(const std::string &name);

    virtual ~AlignWriter();

    /// Write misalignment.
    //  The misalignment is read from the given AlignWrapper.
    virtual void misalignment(const AlignWrapper &, int occur);

private:

    // Not implemented.
    AlignWriter();
    AlignWriter(const AlignWriter &);
    void operator=(const AlignWriter &);
};

#endif // OPAL_AlignWriter_HH
