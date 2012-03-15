// ------------------------------------------------------------------------
// $RCSfile: Geometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Geometry
//   The class for the OPAL GEOMETRY command.
//
// $Date: 2003/08/11 22:09:00 $
// $Author: A. Adelmann $
//
// ------------------------------------------------------------------------

//FIXME: cleanup
#include "Structure/BoundaryGeometry.h"
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

// Class BoundaryGeometry
// ------------------------------------------------------------------------

// The attributes of class Geometry.
namespace {
  enum {
    FGEOM,     // file holding the geometry
    LENGTH,   // length of elliptic tube
    S,        // start of the geometry
    A,		    // major semi-axis of elliptic tube
    B,	      // minor semi-axis of ellitpic tube 
    TOPO,     // BOX, ELLIPTIC if FGEOM is selected topo is over-written 
    SIZE
  };
}

BoundaryGeometry::BoundaryGeometry() : Definition(SIZE, "GEOMETRY", "The \"GEOMETRY\" statement defines the beam pipe geometry.")
{
  itsAttr[FGEOM] = Attributes::makeString
    ("FGEOM", "Specifies the geometry file [h5fed]", "");

  itsAttr[TOPO] = Attributes::makeString
    ("TOPO", "BOX, ELLIPTIC if FGEOM is selected topo is over-written ", "ELLIPTIC");
  
  itsAttr[LENGTH] = Attributes::makeReal
    ("LENGTH", "Specifies the length of a tube shaped elliptic beam pipe [m]", 1.0);
  
  itsAttr[S] = Attributes::makeReal
    ("S", "Specifies the start of a tube shaped elliptic beam pipe [m]", 0.0);
    
  itsAttr[A] = Attributes::makeReal
    ("A", "Specifies the major semi-axis of a tube shaped elliptic beam pipe [m]", 0.025);
    
  itsAttr[B] = Attributes::makeReal
    ("B", "Specifies the minor semi-axis of a tube shaped elliptic beam pipe [m]", 0.025);
    
  BoundaryGeometry *defGeometry = clone("UNNAMED_GEOMETRY");
  defGeometry->builtin = true;
  
  try {
    defGeometry->update();
    OPAL.define(defGeometry);
  } 
  catch (...) {
    delete defGeometry;
  }
}


BoundaryGeometry::BoundaryGeometry(const string &name, BoundaryGeometry *parent):
  Definition(name, parent)
{
  //empty so far
}


BoundaryGeometry::~BoundaryGeometry()
{
  //empty so far
}


bool BoundaryGeometry::canReplaceBy(Object *object)
{
  // Can replace only by another GEOMETRY.
  return dynamic_cast<Geometry *>(object) != 0;
}


BoundaryGeometry *BoundaryGeometry::clone(const string &name)
{
  return new BoundaryGeometry(name, this);
}


void BoundaryGeometry::execute()
{
  update();
}


BoundaryGeometry *BoundaryGeometry::find(const string &name)
{
  BoundaryGeometry *geom = dynamic_cast<BoundaryGeometry*>(OPAL.find(name));

  if(geom == 0) 
    throw OpalException("BoundaryGeometry::find()", "Geometry \"" + name + "\" not found.");
  
  return geom;
}


string BoundaryGeometry::getFilename()
{
  return (string)Attributes::getString(itsAttr[FGEOM]);
}

string BoundaryGeometry::getTopology()
{
  return (string)Attributes::getString(itsAttr[TOPO]);
}

double BoundaryGeometry::getA()
{
  return (double)Attributes::getReal(itsAttr[A]);
}

double BoundaryGeometry::getB()
{
  return (double)Attributes::getReal(itsAttr[B]);
}

void BoundaryGeometry::update()
{
  // Set default name.
  if (getOpalName().empty()) setOpalName("UNNAMED_GEOMETRY");
}

/*
void Geometry::initGeometryfunction(ElementBase &element)
{
  *gmsg << "* ************* G E O M E T R Y ****************************************************"<< endl;
  *gmsg << "Geometry::initGeometryfunction " << endl;
  *gmsg << "* **********************************************************************************"<< endl;

  itsElement_m = &element;
  if (Attributes::getString(itsAttr[FILE]) != "") {
    wf_m = new CSRGeometryFunction(reference, itsElement_m);
  }
  else if (Attributes::getString(itsAttr[TYPE]) == "PARAMETRIZED") {
    double acMode = 1;
  }
  else {
    INFOMSG("no geometry attached" << endl);
  }
}
*/

Inform &BoundaryGeometry::print(Inform &os) const
{
  os << "* ************* W A K E ************************************************************ " << endl;
  os << "* GEOMETRY     " << getOpalName() << '\n'
     << "* FGEOM        " << Attributes::getString(itsAttr[FGEOM]) << '\n'
     << "* TOPO         " << Attributes::getString(itsAttr[TOPO]) << '\n'
     << "* LENGTH       " << Attributes::getReal(itsAttr[LENGTH]) << '\n'
     << "* S            " << Attributes::getReal(itsAttr[S]) << '\n'
     << "* A            " << Attributes::getReal(itsAttr[A]) << '\n'
     << "* B            " << Attributes::getReal(itsAttr[B]) << endl ;
  os << "* ********************************************************************************** " << endl;
}
