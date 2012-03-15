// ------------------------------------------------------------------------
// $RCSfile: Wake.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Wake
//   The class for the OPAL WAKE command.
//
// $Date: 2003/08/11 22:09:00 $
// $Author: A. Adelmann $
//
// ------------------------------------------------------------------------

#include "Structure/Wake.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

using namespace Expressions;
using namespace Physics;


// Class Wake
// ------------------------------------------------------------------------

// The attributes of class Wake.
namespace {
  enum {
    // DESCRIPTION OF SINGLE PARTICLE:
    TYPE,       // The type of the wake
    NBIN,       // Number of bins for the line density
    SIZE
  };
}

Wake::Wake():
  Definition(SIZE, "WAKE",
	     "The \"WAKE\" statement defines data for a the wakefuction "
	     "on an element.")
{
  itsAttr[TYPE] = Attributes::makeString
    ("TYPE", "Specifies the wake function: 1D-CSR, LONG-SHORT-RANGE, TRANSV-SHORT-RANGE, LONG-TRANSV-SHORT-RANGE");
  
  itsAttr[NBIN] = Attributes::makeReal
    ("NBIN", "Number of bins for the line density calculation");
  
  Wake *defWake = clone("UNNAMED_WAKE");
  defWake->builtin = true;
  
  try {
    defWake->update();
    OPAL.define(defWake);
  } catch (...) {
    delete defWake;
  }
}


Wake::Wake(const string &name, Wake *parent):
  Definition(name, parent)
{}


Wake::~Wake()
{}


bool Wake::canReplaceBy(Object *object)
{
  // Can replace only by another WAKE.
  return dynamic_cast<Wake *>(object) != 0;
}


Wake *Wake::clone(const string &name)
{
  return new Wake(name, this);
}


void Wake::execute()
{
  update();
}


Wake *Wake::find(const string &name)
{
  Wake *wake = dynamic_cast<Wake*>(OPAL.find(name));

  if (wake == 0) {
    throw OpalException("Wake::find()", "Wake \"" + name + "\" not found.");
  }
  return wake;
}


int Wake::getNumberOfBins()
{
  return (int)Attributes::getReal(itsAttr[NBIN]);
}


void Wake::update()
{
  // Set default name.
  if (getOpalName().empty()) setOpalName("UNNAMED_WAKE");
}


void Wake::tfsDescriptors(std::ostream &os) const
{
  os << "@ WAKE     %s  " << getOpalName() << '\n'
     << "@ BINS     %le " << Attributes::getReal(itsAttr[NBIN]) << '\n';
}


Inform &Wake::print(Inform &os) const
{
  os << "* ************* W A K E ************************************************************ " << endl;
  os << "* WAKE        " << getOpalName() << '\n'
     << "* BINS        " << Attributes::getReal(itsAttr[NBIN]) << '\n';
  os << "* ********************************************************************************** " << endl;
}
