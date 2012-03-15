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
#include "Solvers/GreenWakeFunction.hh"
#include "Solvers/CSRWakeFunction.hh"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ElementBase.h"


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
    CONST_LENGTH,// True if the length of the Bunch is considered as constant
    CONDUCT,	// Conductivity, either AC or DC
    Z0,		//
    FORM,	// From of the tube 
    RADIUS,	// Radius of the tube
    SIGMA,
    TAU,
    SIZE
  };
}

Wake::Wake():
  Definition(SIZE, "WAKE",
             "The \"WAKE\" statement defines data for a the wakefuction "
             "on an element."),
  wf_m(0)
{
  itsAttr[TYPE] = Attributes::makeString
    ("TYPE", "Specifies the wake function: 1D-CSR, LONG-SHORT-RANGE, TRANSV-SHORT-RANGE, LONG-TRANSV-SHORT-RANGE");
  
  itsAttr[NBIN] = Attributes::makeReal
    ("NBIN", "Number of bins for the line density calculation");

  itsAttr[CONST_LENGTH] = Attributes::makeBool
    ("CONST_LENGTH", "True if the length of the Bunch is considered as constant");
    
  itsAttr[CONDUCT] = Attributes::makeString
    ("CONDUCT", "Cundivity: DC, AC");
    
  itsAttr[Z0] = Attributes::makeReal
    ("Z0", "Impedanz of the beam pipe ");
    
  itsAttr[FORM] = Attributes::makeString
    ("FORM", "The form of the  beam pipe: ROUND");
        
  itsAttr[RADIUS] = Attributes::makeReal
    ("RADIUS", "The radius of the beam pipe [m]");
 
  itsAttr[SIGMA] = Attributes::makeReal
    ("SIGMA", "Material constant dependant on the  beam pipe material");
 
  itsAttr[TAU] = Attributes::makeReal
    ("TAU", "Material constant dependant on the  beam pipe material");
    
  Wake *defWake = clone("UNNAMED_WAKE");
  defWake->builtin = true;
  
  try 
    {
      defWake->update();
      OPAL.define(defWake);
    } 
  catch (...) 
    {
      delete defWake;
    }
}


Wake::Wake(const string &name, Wake *parent):
  Definition(name, parent),
  wf_m(0)
{}


Wake::~Wake()
{
  if (wf_m)
    delete wf_m;
}


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

void Wake::initWakefunction(ElementBase &element)
{
  *gmsg << "* ************* W A K E ************************************************************"<< endl;
  *gmsg << "Wake::initWakefunction " << endl;
  *gmsg << "* **********************************************************************************"<< endl;

  itsElement_m = &element;
  if (Attributes::getString(itsAttr[TYPE]) == "1D-CSR") {
    wf_m = new CSRWakeFunction(reference, itsElement_m);
  }
  else if (Attributes::getString(itsAttr[TYPE]) == "LONG-SHORT-RANGE") {
    double acMode = 1;
    if (Attributes::getString(itsAttr[CONDUCT]).compare("DC")==0) {
      acMode = 2;
    }
    wf_m = new GreenWakeFunction(reference, (int)(Attributes::getReal(itsAttr[NBIN])), (Attributes::getReal(itsAttr[Z0])), (Attributes::getReal(itsAttr[RADIUS])), 
				 (Attributes::getReal(itsAttr[SIGMA])), Physics::c, acMode, 
				 (Attributes::getReal(itsAttr[TAU])), 1, (Attributes::getBool(itsAttr[CONST_LENGTH])));     
  }
  else if (Attributes::getString(itsAttr[TYPE]) == "TRANSV-SHORT-RANGE") {
    double acMode = 1;
    if (Attributes::getString(itsAttr[CONDUCT]).compare("DC")==0) {
      acMode = 2;
    }
   
    wf_m = new GreenWakeFunction(reference, (int)(Attributes::getReal(itsAttr[NBIN])), (Attributes::getReal(itsAttr[Z0])), (Attributes::getReal(itsAttr[RADIUS])), 
				 (Attributes::getReal(itsAttr[SIGMA])), Physics::c, acMode, 
				 (Attributes::getReal(itsAttr[TAU])), 0, (Attributes::getBool(itsAttr[CONST_LENGTH])));     
  }
  else if (Attributes::getString(itsAttr[TYPE]) == "LONG-TRANSV-SHORT-RANGE") {
    
  }
  else {
    wf_m = 0;
    INFOMSG("no wakefunction attached" << endl);
  }
}

Inform &Wake::print(Inform &os) const
{
  os << "* ************* W A K E ************************************************************ " << endl;
  os << "* WAKE         " << getOpalName() << '\n'
     << "* BINS         " << Attributes::getReal(itsAttr[NBIN]) << '\n'
     << "* CONST_LENGTH " << Attributes::getReal(itsAttr[CONST_LENGTH]) << '\n'
     << "* CONDUCT      " << Attributes::getReal(itsAttr[CONDUCT]) << '\n'
     << "* Z0           " << Attributes::getReal(itsAttr[Z0]) << '\n'
     << "* FORM         " << Attributes::getReal(itsAttr[FORM]) << '\n'
     << "* RADIUS       " << Attributes::getReal(itsAttr[RADIUS]) << '\n'
     << "* SIGMA        " << Attributes::getReal(itsAttr[SIGMA]) << '\n'
     << "* TAU          " << Attributes::getReal(itsAttr[TAU]) << '\n';     
  os << "* ********************************************************************************** " << endl;
}
