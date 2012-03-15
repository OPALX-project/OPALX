#ifndef OPAL_AlignHandler_HH
#define OPAL_AlignHandler_HH

// ------------------------------------------------------------------------
// $RCSfile: AlignHandler.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignHandler
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"
#include "Expressions/SFunction.h"

class AlignWrapper;
class Beamline;
class ElementBase;
class AlignBase;
class FlaggedElmPtr;


// Class AlignHandler
// ------------------------------------------------------------------------
/// Access misalignment errors.
//  A ``Visitor'' class giving access to the misalignment errors.  Its
//  constructor expects a beam line and and an object ``base'' derived
//  from  ``AlignBase''.  When executed, it calls ``base.misalignment()''
//  for each element which has a misalignment.  The latter may do whatever
//  is suitable to the element's misalignment.

class AlignHandler: public DefaultVisitor {

public:

    /// Constructor.
    //  Set up the visitor to use the given beam line [b]bl[/b] and to apply
    //  the command [b]base[/b] to its contained elements.  If [b]full[/b] is
    //  true, apply to all elements, otherwise only to selected elemens.
    AlignHandler(const Beamline &bl, AlignBase &base, bool full);

    virtual ~AlignHandler();

    /// Apply visitor to AlignWrapper.
    //  Access the misalignment.
    virtual void visitAlignWrapper(const AlignWrapper &);

    /// Apply visitor to FlaggedElmPtr.
    //  Makes sure the selection flag is set.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

private:

    // Not implemented.
    AlignHandler();
    AlignHandler(const AlignHandler &);
    void operator=(const AlignHandler &);

    // Override default method for elements.
    virtual void applyDefault(const ElementBase &);

    // Misalignment command.
    AlignBase &command;

    // The S-function handler associated with this handler.
    SFunction sHandler;

    // If true, the handler is applied to all positions.
    bool applyToAll;

    // Keep track of selection.
    bool isSelected;
    int occurCount;
};

#endif // OPAL_AlignHandler_HH
