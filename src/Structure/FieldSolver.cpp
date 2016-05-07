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
#include "Solvers/FFTPoissonSolver.h"
#include "Solvers/FFTBoxPoissonSolver.h"
#include "Solvers/P3MPoissonSolver.h"
#ifdef HAVE_SAAMG_SOLVER
#include "Solvers/MGPoissonSolver.h"
#endif
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "BoundaryGeometry.h"
#include "AbstractObjects/Element.h"
#include "Algorithms/PartBunch.h"

using namespace Expressions;
using namespace Physics;

//TODO: o add a FIELD for DISCRETIZATION, MAXITERS, TOL...

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
        BCFFTX,     // boundary condition in x [FFT only]
        BCFFTY,     // boundary condition in y [FFT only]
        BCFFTT,     // boundary condition in z [FFT only]
        GREENSF,    // holds greensfunction to be used [FFT only]
        BBOXINCR,   // how much the boundingbox is increased
        GEOMETRY,   // geometry of boundary [SAAMG only]
        ITSOLVER,   // iterative solver [SAAMG only]
        INTERPL,    // interpolation used for boundary points [SAAMG only]
        TOL,        // tolerance of the SAAMG preconditioned solver [SAAMG only]
        MAXITERS,   // max number of iterations [SAAMG only]
        PRECMODE,   // preconditioner mode [SAAMG only]
        RC,         // cutoff radius for PP interactions
	ALPHA,      // Green’s function splitting parameter
	EPSILON,    // regularization for PP interaction
        // FOR XXX BASED SOLVER
        SIZE
    };
}


FieldSolver::FieldSolver():
    Definition(SIZE, "FIELDSOLVER",
               "The \"FIELDSOLVER\" statement defines data for a the field solver ") {

    itsAttr[FSTYPE] = Attributes::makeString("FSTYPE", "Name of the attached field solver: FFT, FFTPERIODIC, SAAMG, AMR, and NONE ");

    itsAttr[MX] = Attributes::makeReal("MX", "Meshsize in x");
    itsAttr[MY] = Attributes::makeReal("MY", "Meshsize in y");
    itsAttr[MT] = Attributes::makeReal("MT", "Meshsize in z(t)");

    itsAttr[PARFFTX] = Attributes::makeBool("PARFFTX", "True, dimension 0 i.e x is parallelized", false);
    itsAttr[PARFFTY] = Attributes::makeBool("PARFFTY", "True, dimension 1 i.e y is parallelized", false);
    itsAttr[PARFFTT] = Attributes::makeBool("PARFFTT", "True, dimension 2 i.e z(t) is parallelized", true);

    //FFT ONLY:
    itsAttr[BCFFTX] = Attributes::makeString("BCFFTX", "Boundary conditions in x: open, dirichlet (box) ");
    itsAttr[BCFFTY] = Attributes::makeString("BCFFTY", "Boundary conditions in y: open, dirichlet (box) ");
    itsAttr[BCFFTT] = Attributes::makeString("BCFFTT", "Boundary conditions in z(t): open, periodoc");

    itsAttr[GREENSF]  = Attributes::makeString("GREENSF", "Which Greensfunction to be used [STANDARD | INTEGRATED]", "INTEGRATED");
    itsAttr[BBOXINCR] = Attributes::makeReal("BBOXINCR", "Increase of bounding box in % ", 2.0);

    // P3M only:
    itsAttr[RC]  = Attributes::makeReal("RC", "cutoff radius for PP interactions",0.0);
    itsAttr[ALPHA]  = Attributes::makeReal("ALPHA", "Green’s function splitting parameter",0.0);
    itsAttr[EPSILON]  = Attributes::makeReal("EPSILON", "regularization for PP interaction",0.0);

    //SAAMG and in case of FFT with dirichlet BC in x and y
    itsAttr[GEOMETRY] = Attributes::makeString("GEOMETRY", "GEOMETRY to be used as domain boundary", "");
    itsAttr[ITSOLVER]  = Attributes::makeString("ITSOLVER", "Type of iterative solver [CG | BiCGSTAB | GMRES]", "CG");
    itsAttr[INTERPL]  = Attributes::makeString("INTERPL", "interpolation used for boundary points [CONSTANT | LINEAR | QUADRATIC]", "LINEAR");
    itsAttr[TOL] = Attributes::makeReal("TOL", "Tolerance for iterative solver", 1e-8);
    itsAttr[MAXITERS] = Attributes::makeReal("MAXITERS", "Maximum number of iterations of iterative solver", 100);
    itsAttr[PRECMODE]  = Attributes::makeString("PRECMODE", "Preconditioner Mode [STD | HIERARCHY | REUSE]", "HIERARCHY");

    mesh_m = 0;
    FL_m = 0;
    PL_m = 0;

    solver_m = 0;
}


FieldSolver::FieldSolver(const std::string &name, FieldSolver *parent):
    Definition(name, parent)
{
    mesh_m = 0;
    FL_m = 0;
    PL_m = 0;
    solver_m = 0;
}


FieldSolver::~FieldSolver() {
    if (mesh_m) {
        delete mesh_m;
        mesh_m = 0;
    }
    if (FL_m) {
        delete FL_m;
        FL_m = 0;
    }
    if (solver_m) {
       delete solver_m;
       solver_m = 0;
    }
}

FieldSolver *FieldSolver::clone(const std::string &name) {
    return new FieldSolver(name, this);
}

void FieldSolver::execute() {
    update();
}

FieldSolver *FieldSolver::find(const std::string &name) {
    FieldSolver *fs = dynamic_cast<FieldSolver *>(OpalData::getInstance()->find(name));

    if(fs == 0) {
        throw OpalException("FieldSolver::find()", "FieldSolver \"" + name + "\" not found.");
    }
    return fs;
}


double FieldSolver::getMX() const {
    return Attributes::getReal(itsAttr[MX]);
}

double FieldSolver::getMY() const {
    return Attributes::getReal(itsAttr[MY]);
}

double FieldSolver::getMT() const {
    return Attributes::getReal(itsAttr[MT]);
}

void FieldSolver::setMX(double value) {
    Attributes::setReal(itsAttr[MX], value);
}

void FieldSolver::setMY(double value) {
    Attributes::setReal(itsAttr[MY], value);
}

void FieldSolver::setMT(double value) {
    Attributes::setReal(itsAttr[MT], value);
}

void FieldSolver::update() {

}

void FieldSolver::initCartesianFields() {

    e_dim_tag decomp[3] = {SERIAL, SERIAL, SERIAL};

    NDIndex<3> domain;
    domain[0] = Index((int)getMX() + 1);
    domain[1] = Index((int)getMY() + 1);
    domain[2] = Index((int)getMT() + 1);

    if(Attributes::getBool(itsAttr[PARFFTX]))
        decomp[0] = PARALLEL;
    if(Attributes::getBool(itsAttr[PARFFTY]))
        decomp[1] = PARALLEL;
    if(Attributes::getBool(itsAttr[PARFFTT]))
        decomp[2] = PARALLEL;

    if(Attributes::getString(itsAttr[FSTYPE]) == "FFTPERIODIC") {
        decomp[0] = decomp[1] = SERIAL;
        decomp[2] = PARALLEL;
    }
    // create prototype mesh and layout objects for this problem domain
    mesh_m   = new Mesh_t(domain);
    FL_m     = new FieldLayout_t(*mesh_m, decomp);
    PL_m     = new Layout_t(*FL_m, *mesh_m);
    // OpalData::getInstance()->setMesh(mesh_m);
    // OpalData::getInstance()->setFieldLayout(FL_m);
    // OpalData::getInstance()->setLayout(PL_m);
}

bool FieldSolver::hasPeriodicZ() {
  return Attributes::getString(itsAttr[BCFFTT])==std::string("PERIODIC");
}

bool FieldSolver::isAMRSolver() {
  return Attributes::getString(itsAttr[FSTYPE])==std::string("AMR");
}


void FieldSolver::initSolver(PartBunch &b) {
  itsBunch_m = &b;
  std::string bcx = Attributes::getString(itsAttr[BCFFTX]);
  std::string bcy = Attributes::getString(itsAttr[BCFFTY]);
  std::string bcz = Attributes::getString(itsAttr[BCFFTT]);
  
  if(Attributes::getString(itsAttr[FSTYPE]) == "FFT") { 
    bool sinTrafo = ((bcx == std::string("DIRICHLET")) && (bcy == std::string("DIRICHLET")) && (bcz == std::string("DIRICHLET")));
    if(sinTrafo) {
      std::cout << "FFTBOX ACTIVE" << std::endl;
      //we go over all geometries and add the Geometry Elements to the geometry list
      std::string geoms = Attributes::getString(itsAttr[GEOMETRY]);
      std::string tmp = "";
      //split and add all to list
      std::vector<BoundaryGeometry *> geometries;
      for(unsigned int i = 0; i <= geoms.length(); i++) {
	if(geoms[i] == ',' || i == geoms.length()) {
	  BoundaryGeometry *geom = BoundaryGeometry::find(tmp);
	  if(geom != 0)
	    geometries.push_back(geom);
	  tmp.clear();
	} else
	  tmp += geoms[i];
      }
      BoundaryGeometry *ttmp = geometries[0];
      solver_m = new FFTBoxPoissonSolver(mesh_m, FL_m, Attributes::getString(itsAttr[GREENSF]), ttmp->getA());
      itsBunch_m->set_meshEnlargement(Attributes::getReal(itsAttr[BBOXINCR]) / 100.0);
      fsType_m = "FFTBOX";
    } else {
      solver_m = new FFTPoissonSolver(mesh_m, FL_m, Attributes::getString(itsAttr[GREENSF]), bcz);
      itsBunch_m->set_meshEnlargement(Attributes::getReal(itsAttr[BBOXINCR]) / 100.0);
      fsType_m = "FFT";
    }
  } else if (Attributes::getString(itsAttr[FSTYPE]) == "P3M") {
    solver_m = new P3MPoissonSolver(mesh_m, FL_m, Attributes::getString(itsAttr[GREENSF]), bcz);
    fsType_m = "P3M";
  } else if(Attributes::getString(itsAttr[FSTYPE]) == "SAAMG") {
#ifdef HAVE_SAAMG_SOLVER
    //we go over all geometries and add the Geometry Elements to the geometry list
    std::string geoms = Attributes::getString(itsAttr[GEOMETRY]);
    std::string tmp = "";
    //split and add all to list
    std::vector<BoundaryGeometry *> geometries;
    for(unsigned int i = 0; i <= geoms.length(); i++) {
      if(geoms[i] == ',' || i == geoms.length()) {
	BoundaryGeometry *geom = OpalData::getInstance()->getGlobalGeometry();
	if(geom != 0) {
	  geometries.push_back(geom);
	}
	tmp.clear();
      } else
	tmp += geoms[i];
    }
    solver_m = new MGPoissonSolver(b, mesh_m, FL_m, geometries, Attributes::getString(itsAttr[ITSOLVER]),
				   Attributes::getString(itsAttr[INTERPL]),
				   Attributes::getReal(itsAttr[TOL]),
				   Attributes::getReal(itsAttr[MAXITERS]),
				   Attributes::getString(itsAttr[PRECMODE]));
    itsBunch_m->set_meshEnlargement(Attributes::getReal(itsAttr[BBOXINCR]) / 100.0);
    fsType_m = "SAAMG";
#else
    INFOMSG("SAAMG Solver not enabled! Please build OPAL with -DENABLE_SAAMG_SOLVER=1" << endl);
    INFOMSG("switching to FFT solver..." << endl);
    solver_m = new FFTPoissonSolver(mesh_m, FL_m, Attributes::getString(itsAttr[GREENSF]),bcz);
    fsType_m = "FFT";
#endif
  } 
  else {
    solver_m = 0;
    INFOMSG("no solver attached" << endl);
  }    
}

bool FieldSolver::hasValidSolver() {
    return (solver_m != 0);
}

Inform &FieldSolver::printInfo(Inform &os) const {
    std::string fsType;
    if (Attributes::getString(itsAttr[BCFFTT])==std::string("PERIODIC"))
	fsType = Attributes::getString(itsAttr[FSTYPE])+"-zPeriodic";
    else
	fsType = Attributes::getString(itsAttr[FSTYPE]);

    os << "* ************* F I E L D S O L V E R ********************************************** " << endl;
    os << "* FIELDSOLVER  " << getOpalName() << '\n'
       << "* TYPE         " << fsType << '\n'
       << "* N-PROCESSORS " << Ippl::getNodes() << '\n'
       << "* MX           " << Attributes::getReal(itsAttr[MX])   << '\n'
       << "* MY           " << Attributes::getReal(itsAttr[MY])   << '\n'
       << "* MT           " << Attributes::getReal(itsAttr[MT])   << '\n'
       << "* BBOXINCR     " << Attributes::getReal(itsAttr[BBOXINCR]) << endl;

    if(fsType == "P3M")
      
      os << "* RC         " << Attributes::getReal(itsAttr[RC]) << '\n'
	 << "* ALPHA      " << Attributes::getReal(itsAttr[ALPHA]) << '\n'
	 << "* EPSILON    " << Attributes::getReal(itsAttr[EPSILON]) << endl;
    

    if(fsType == "FFT") {
        os << "* GRRENSF      " << Attributes::getString(itsAttr[GREENSF]) << endl;
    } else if (fsType != "NONE") {
        os << "* GEOMETRY     " << Attributes::getString(itsAttr[GEOMETRY]) << '\n'
           << "* ITSOLVER     " << Attributes::getString(itsAttr[ITSOLVER])   << '\n'
           << "* INTERPL      " << Attributes::getString(itsAttr[INTERPL])  << '\n'
           << "* TOL          " << Attributes::getReal(itsAttr[TOL])        << '\n'
           << "* MAXITERS     " << Attributes::getReal(itsAttr[MAXITERS]) << '\n'
           << "* PRECMODE     " << Attributes::getString(itsAttr[PRECMODE])   << endl;
    }

    if(Attributes::getBool(itsAttr[PARFFTX]))
        os << "* XDIM         parallel  " << endl;
    else
        os << "* XDIM         serial  " << endl;

    if(Attributes::getBool(itsAttr[PARFFTY]))
        os << "* YDIM         parallel  " << endl;
    else
        os << "* YDIM         serial  " << endl;

    if(Attributes::getBool(itsAttr[PARFFTT]))
        os << "* Z(T)DIM      parallel  " << endl;
    else
        os << "* Z(T)DIM      serial  " << endl;

    INFOMSG(level3 << *mesh_m << endl);
    INFOMSG(level3 << *PL_m << endl);
    if(solver_m)
        os << *solver_m << endl;
    os << "* ********************************************************************************** " << endl;
    return os;
}
