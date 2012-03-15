// ------------------------------------------------------------------------
// $RCSfile: FieldSolver.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FieldSolver
//   The class for the OPAL FIELDSOLVER command.
//
// ------------------------------------------------------------------------
//
// $Date: 2003/08/11 22:09:00 $
// $Author: ADA $
//
// ------------------------------------------------------------------------

#include "Structure/FieldSolver.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

using namespace Expressions;
using namespace Physics;


// Class FieldSolver
// ------------------------------------------------------------------------

// The attributes of class FieldSolver.
namespace {
  enum {
    FSTYPE,   // The field solver name
    // FOR FFT BASED SOLVER
    MX,         // mesh sixe in x
    MY,         // mesh sixe in y
    MT,         //  mesh sixe in z
    PARFFTX,    // parallelized grind in x
    PARFFTY,    // parallelized grind in y
    PARFFTT,    // parallelized grind in z
    BCFFTX,     // boundary condition in x
    BCFFTY,     // boundary condition in y
    BCFFTT,     // boundary condition in z
    GREENSF,    // holds greensfunction to be used
    BBOXINCR,   // how much the boundingbox is increased
    // FOR XXX BASED SOLVER
    SIZE
  };
}


FieldSolver::FieldSolver():
  Definition(SIZE, "FIELDSOLVER",
	     "The \"FIELDSOLVER\" statement defines data for a the field solver ")
{

  itsAttr[FSTYPE] = Attributes::makeString("FSTYPE", "Name of the attached field solver: FFT, xxxx ");

  itsAttr[MX] = Attributes::makeReal("MX", "Meshsize in x");
  itsAttr[MY] = Attributes::makeReal("MY", "Meshsize in y");
  itsAttr[MT] = Attributes::makeReal("MT", "Meshsize in z(t)");

  itsAttr[PARFFTX] = Attributes::makeBool("PARFFTX", "True, dimension 0 i.e x is parallelized",false);
  itsAttr[PARFFTY] = Attributes::makeBool("PARFFTY", "True, dimension 1 i.e y is parallelized",false);
  itsAttr[PARFFTT] = Attributes::makeBool("PARFFTT", "True, dimension 2 i.e z(t) is parallelized",true);

  itsAttr[BCFFTX] = Attributes::makeString("BCFFTX", "Boundary conditions in x: open, parallel ");
  itsAttr[BCFFTY] = Attributes::makeString("BCFFTY", "Boundary conditions in y: open, parallel ");
  itsAttr[BCFFTT] = Attributes::makeString("BCFFTT", "Boundary conditions in z(t): open, parallel");

  itsAttr[GREENSF]  = Attributes::makeString("GREENSF", "Which Greensfunction to be used [STANDARD | INTEGRATED]","INTEGRATED");
  itsAttr[BBOXINCR] = Attributes::makeReal("BBOXINCR", "Increase of bounding box in \% ",2.0);
}


FieldSolver::FieldSolver(const string &name, FieldSolver *parent):
  Definition(name, parent)
{}


FieldSolver::~FieldSolver()
{}

FieldSolver *FieldSolver::clone(const string &name)
{
  return new FieldSolver(name, this);
}

void FieldSolver::execute()
{
  update();
}

FieldSolver *FieldSolver::find(const string &name)
{
  FieldSolver *fs = dynamic_cast<FieldSolver*>(OPAL.find(name));

  if (fs == 0) {
    throw OpalException("FieldSolver::find()", "FieldSolver \"" + name + "\" not found.");
  }
  return fs;
}


double FieldSolver::getMX() const
{
  return Attributes::getReal(itsAttr[MX]);
}

double FieldSolver::getMY() const
{
  return Attributes::getReal(itsAttr[MY]);
}

double FieldSolver::getMT() const
{
  return Attributes::getReal(itsAttr[MT]);
}

void FieldSolver::setMX(double value)
{
  Attributes::setReal(itsAttr[MX], value);
}


void FieldSolver::setMY(double value)
{
  Attributes::setReal(itsAttr[MY], value);
}

void FieldSolver::setMT(double value)
{
  Attributes::setReal(itsAttr[MT], value);
}


void FieldSolver::update()
{

}

void FieldSolver::initCartesianFields() 
{
  
  e_dim_tag decomp[3] = {SERIAL,SERIAL,SERIAL};
  
  NDIndex<3> domain;
  domain[0] = Index( (int)getMX() + 1 );
  domain[1] = Index( (int)getMY() + 1 );
  domain[2] = Index( (int)getMT() + 1 );
  
  if (Attributes::getBool(itsAttr[PARFFTX]))
    decomp[0] = PARALLEL;
  if (Attributes::getBool(itsAttr[PARFFTY]))
    decomp[1] = PARALLEL;
  if (Attributes::getBool(itsAttr[PARFFTT]))
    decomp[2] = PARALLEL;
  
  if (Attributes::getString(itsAttr[FSTYPE]) == "FFTPERIODIC") {
    decomp[0] = decomp[1] = SERIAL;
    decomp[2] = PARALLEL;
  }
  // create prototype mesh and layout objects for this problem domain
  mesh_m   = new Mesh_t(domain);
  FL_m     = new FieldLayout_t(*mesh_m, decomp);
  PL_m     = new Layout_t(*FL_m, *mesh_m);
}

void FieldSolver::initSolver(PartBunch &b) 
{
  itsBunch_m = &b;
  
  if(Attributes::getString(itsAttr[FSTYPE])=="FFT") {
    solver_m = new FFTPoissonSolver(mesh_m,FL_m,Attributes::getString(itsAttr[GREENSF]));
  }
  else {
    solver_m = 0;
    INFOMSG("no solver attached" << endl); // solver_m = 0;
  }
  itsBunch_m->set_meshEnlargement(Attributes::getReal(itsAttr[BBOXINCR])/100.0);
}


bool FieldSolver::hasValidSolver() {
  return (solver_m != 0);
}

Inform &FieldSolver::print(Inform &os) const
{
  os << "* ************* F I E L D S O L V E R ********************************************** " << endl;
  os << "* FIELDSOLVER  " << getOpalName() << '\n'
     << "* TYPE         " << Attributes::getString(itsAttr[FSTYPE]) << '\n'
     << "* N-PROCESSORS " << Ippl::getNodes() << '\n'
     << "* MX           " << Attributes::getReal(itsAttr[MX])   << '\n'
     << "* MY           " << Attributes::getReal(itsAttr[MY])   << '\n'
     << "* MT           " << Attributes::getReal(itsAttr[MT])   << '\n'     
     << "* BBOXINCR     " << Attributes::getReal(itsAttr[BBOXINCR]) << endl
     << "* GRRENSF      " << Attributes::getString(itsAttr[GREENSF]) << endl;
  if (Attributes::getBool(itsAttr[PARFFTX]))
    os << "* XDIM is parallel  " << endl;
  else
    os << "* XDIM is serial  " << endl;

  if (Attributes::getBool(itsAttr[PARFFTY]))
    os << "* YDIM is parallel  " << endl;
  else
    os << "* YDIM is serial  " << endl;

  if (Attributes::getBool(itsAttr[PARFFTT]))
    os << "* Z(T)DIM is parallel  " << endl;
  else
    os << "* Z(T)DIM is serial  " << endl;
  
  INFOMSG(*mesh_m << endl);
  INFOMSG(*PL_m << endl);
  if (solver_m)
    os << *solver_m << endl;
  os << "* ********************************************************************************** " << endl;
}
