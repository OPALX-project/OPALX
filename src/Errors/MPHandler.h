#ifndef OPAL_MPHandler_HH
#define OPAL_MPHandler_HH

// ------------------------------------------------------------------------
// $RCSfile: MPHandler.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: MPHandler
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"

class Beamline;
class ElementBase;
class MPBase;
class MultipoleWrapper;
class RBendWrapper;
class SBendWrapper;


// Class MPHandler
// ------------------------------------------------------------------------
/// Access field errors.
//  A ``Visitor'' class giving access to the field errors. Its constructor
//  expects a beam line and and an object ``base'' derived ``MPBase''.
//  When executed, it calls ``base.fieldError()'' for each element which
//  has a field error.  The latter may do whatever is suitable to the
//  element's field error.

class MPHandler: public DefaultVisitor {

public:

  /// Constructor.
  //  Set up the visitor to use the given beam line [b]bl[/b] and to apply
  //  the command [b]base[/b] to its contained elements.  If [b]full[/b] is
  //  true, apply to all elements, otherwise only to selected elemens.
  MPHandler(const Beamline &bl, MPBase &base, bool full);

  virtual ~MPHandler();

  /// Apply visitor to MultipoleWrapper.
  virtual void visitMultipoleWrapper(const MultipoleWrapper &);

  /// Apply visitor to RBendWrapper.
  virtual void visitRBendWrapper(const RBendWrapper &);

  /// Apply visitor to SBendWrapper.
  virtual void visitSBendWrapper(const SBendWrapper &);

  /// Apply visitor to FlaggedElmPtr.
  //  Makes sure the selection flag is set.
  virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

private:

  // Not implemented.
  MPHandler();
  MPHandler(const MPHandler &);
  void operator=(const MPHandler &);

  // Override visit method for all other elements.
  // Make sure that the selection flag is reset.
  virtual void applyDefault(const ElementBase &);
    
  // The command to generate errors.
  MPBase &command;

  // If true, the visitor applies to all positions.
  bool applyToAll;

  // Keep track of selection.
  bool isSelected;
  int occurCount;
};

#endif // OPAL_MPHandler_HH
