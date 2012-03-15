#ifndef OPAL_LineWriter_HH
#define OPAL_LineWriter_HH 1

// ------------------------------------------------------------------------
// $RCSfile: LineWriter.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LineWriter
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2001/08/13 15:16:16 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/DoomWriter.h"
#include "Algorithms/DefaultVisitor.h"

class Beamline;
class Drift;
class ElementBase;
class FlaggedElmPtr;


// Class LineWriter
// ------------------------------------------------------------------------
/// Write a beam line or sequence list to the DOOM data base.

class LineWriter: public DefaultVisitor {

public:

    LineWriter(const Beamline &, const string &key, const string &pName);
    ~LineWriter();

    /// Override drift operation.
    //  Accumulate length, but do not save anything.
    virtual void visitDrift(const Drift &);

    /// Override beamline exit.
    virtual void visitMapIntegrator(const MapIntegrator &);

    /// Apply visitor to FlaggedElmPtr.
    //  Makes sure the selection flag is set.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

protected:

    /// Apply default.
    //  Add the element to the DOOM list.
    virtual void applyDefault(const ElementBase &);

private:

    // The Doom object to be filled with the LineWriter list.
    // Its destructor will write the object.
    DoomWriter itsWriter;

    // The current attribute position in the DoomWriter block.
    int itsCount;

    // The current longitudinal position.
    double itsPosition;

    // Keep track of the occurrence counts.
    int occurCount;
};

#endif // OPAL_LineWriter_HH
