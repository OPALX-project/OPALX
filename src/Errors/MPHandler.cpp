// ------------------------------------------------------------------------
// $RCSfile: MPHandler.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: MPHandler
//   This Visitor class provides access to multipole errors.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/MPHandler.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/SBend.h"
#include "AbstractObjects/Table.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "ComponentWrappers/RBendWrapper.h"
#include "ComponentWrappers/SBendWrapper.h"
#include "Fields/BMultipoleField.h"
#include "Errors/MPBase.h"
#include "Utilities/Options.h"


// Class MPHandler
// ------------------------------------------------------------------------

MPHandler::MPHandler(const Beamline &beamline,
		     MPBase &cmd,
		     bool full):
  DefaultVisitor(beamline, false, false),
  command(cmd),
  applyToAll(full)
{}


MPHandler::~MPHandler()
{}
 

void MPHandler::visitMultipoleWrapper(const MultipoleWrapper &wrap)
{
  if (isSelected) {
    command.fieldError(wrap.getName(), occurCount,
		       wrap.getDesign().getField(),
		       wrap.errorField());
  }

  isSelected = false;
}


void MPHandler::visitRBendWrapper(const RBendWrapper &wrap)
{
  if (isSelected) {
    command.fieldError(wrap.getName(), occurCount,
		       wrap.getDesign().getField(),
		       wrap.errorField());
  }

  isSelected = false;
}


void MPHandler::visitSBendWrapper(const SBendWrapper &wrap)
{ 
  if (isSelected) {
    command.fieldError(wrap.getName(), occurCount,
		       wrap.getDesign().getField(),
		       wrap.errorField());
  }

  isSelected = false;
}


void MPHandler::visitFlaggedElmPtr(const FlaggedElmPtr &fep)
{
  // Set selection flag for use in the embedded element.
  isSelected = applyToAll || fep.getSelectionFlag();
  occurCount = fep.getCounter();

  // Forward call to element.
  DefaultVisitor::visitFlaggedElmPtr(fep);
}


void MPHandler::applyDefault(const ElementBase &)
{
  isSelected = false;
}
