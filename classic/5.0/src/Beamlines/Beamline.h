#ifndef CLASSIC_Beamline_HH
#define CLASSIC_Beamline_HH

// ------------------------------------------------------------------------
// $RCSfile: Beamline.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Beamline
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ElementBase.h"

class CLRangeError;


// Class Beamline
// ------------------------------------------------------------------------
/// An abstract sequence of beam line components.
//  A beam line is built as a list of objects derived from  ElmPtr. Each
//  ElmPtr (``element pointer'') points to a  ElementBase, and may contain
//  additional data describing the position, like lattice functions etc.

class Beamline: public ElementBase {

public:

    /// Constructor with given name.
    explicit Beamline(const string &name);

    Beamline();
    Beamline(const Beamline &);
    virtual ~Beamline();

    /// Apply visitor to all elements of the line.
    //  If the parameter [b]reverse[/b] is true, theline is traversed in
    //  reverse direction.  If any error occurs, this method may throw an
    //  exception.
    virtual void iterate(BeamlineVisitor &, bool reverse) const = 0;

private:

    // Not implemented.
    void operator=(const Beamline &);
};

#endif // CLASSIC_Beamline_HH
