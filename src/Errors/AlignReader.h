#ifndef OPAL_AlignReader_HH
#define OPAL_AlignReader_HH

// ------------------------------------------------------------------------
// $RCSfile: AlignReader.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignReader
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/AlignBase.h"


// Class AlignReader
// ------------------------------------------------------------------------
/// DOOM reader for reading alignment errors.
//  Retrieves misalignments from the DOOM data base and uses an AlignHandler
//  to store them in the beam line.

class AlignReader: public AlignBase {

public:

    /// Constructor.
    //  Store [b]name[/b] in the DOOM environment.
    AlignReader(const string &name);

    virtual ~AlignReader();

    /// Read misalignment.
    //  The misalignment is written to the given AlignWrapper.
    virtual void misalignment(const AlignWrapper &, int occur);

private:

    // Not implemented.
    AlignReader();
    AlignReader(const AlignReader &);
    void operator=(const AlignReader &);
};

#endif // OPAL_AlignReader_HH
