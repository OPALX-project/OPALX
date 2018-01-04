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
#include "Utilities/Util.h"
#include "BoundaryGeometry.h"
#include "AbstractObjects/Element.h"
#include "Algorithms/PartBunchBase.h"

#ifdef ENABLE_AMR
    #include "Amr/AmrBoxLib.h"
    #include "Solvers/BoxLibSolvers/FMGPoissonSolver.h"
    #include "Amr/AmrDefs.h"
    #include "Algorithms/AmrPartBunch.h"
#endif

#ifdef HAVE_AMR_MG_SOLVER
    #include "Solvers/TAMG/AmrMultiGrid.h"
#endif

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
        BCFFTX,     // boundary condition in x [FFT + TAMG only]
        BCFFTY,     // boundary condition in y [FFT + TAMG only]
        BCFFTT,     // boundary condition in z [FFT + TAMG only]
        GREENSF,    // holds greensfunction to be used [FFT only]
        BBOXINCR,   // how much the boundingbox is increased
        GEOMETRY,   // geometry of boundary [SAAMG only]
        ITSOLVER,   // iterative solver [SAAMG + TAMG]
        INTERPL,    // interpolation used for boundary points [SAAMG only]
        TOL,        // tolerance of the SAAMG preconditioned solver [SAAMG only]
        MAXITERS,   // max number of iterations [SAAMG only]
        PRECMODE,   // preconditioner mode [SAAMG only]
        RC,         // cutoff radius for PP interactions
	ALPHA,      // Green’s function splitting parameter
	EPSILON,    // regularization for PP interaction
#ifdef ENABLE_AMR
        AMR_MAXLEVEL,       // AMR, maximum refinement level
        AMR_REFX,           // AMR, refinement ratio in x
        AMR_REFY,           // AMR, refinement ratio in y
        AMR_REFT,           // AMR, refinement ration in z
        AMR_SUBCYCLE,       // AMR, subcycling in time for refined levels (default: false)
        AMR_MAXGRID,        // AMR, maximum grid size (default: 16)
        AMR_TAGGING,
        AMR_DENSITY,
        AMR_MAX_NUM_PART,
        AMR_MIN_NUM_PART,
        AMR_SCALING,
#endif
#ifdef HAVE_AMR_MG_SOLVER
        // TAMG = Trilinos-Adaptive-Multi-Grid
        TAMG_SMOOTHER,      // AMR, smoother for level solution
        TAMG_NSWEEPS,       // AMR, number of smoothing sweeps
        TAMG_PREC,          // AMR, preconditioner for bottom solver
        TAMG_INTERP,        // AMR, interpolater for solution from level l to l+1
        TAMG_NORM,          // AMR, norm convergence criteria
#endif
        // FOR XXX BASED SOLVER
        SIZE
    };
}


FieldSolver::FieldSolver():
    Definition(SIZE, "FIELDSOLVER",
               "The \"FIELDSOLVER\" statement defines data for a the field solver ") {

    itsAttr[FSTYPE] = Attributes::makeString("FSTYPE",
                                             "Name of the attached field solver: FFT, FFTPERIODIC, SAAMG, AMR, and NONE ");

    itsAttr[MX] = Attributes::makeReal("MX", "Meshsize in x");
    itsAttr[MY] = Attributes::makeReal("MY", "Meshsize in y");
    itsAttr[MT] = Attributes::makeReal("MT", "Meshsize in z(t)");

    itsAttr[PARFFTX] = Attributes::makeBool("PARFFTX",
                                            "True, dimension 0 i.e x is parallelized",
                                            false);
    
    itsAttr[PARFFTY] = Attributes::makeBool("PARFFTY",
                                            "True, dimension 1 i.e y is parallelized",
                                            false);
    
    itsAttr[PARFFTT] = Attributes::makeBool("PARFFTT",
                                            "True, dimension 2 i.e z(t) is parallelized",
                                            true);

    //FFT ONLY:
    itsAttr[BCFFTX] = Attributes::makeString("BCFFTX",
                                             "Boundary conditions in x: open, dirichlet (box), periodic");
    
    itsAttr[BCFFTY] = Attributes::makeString("BCFFTY",
                                             "Boundary conditions in y: open, dirichlet (box), periodic");
    
    itsAttr[BCFFTT] = Attributes::makeString("BCFFTT",
                                             "Boundary conditions in z(t): open, periodic");

    itsAttr[GREENSF]  = Attributes::makeString("GREENSF",
                                               "Which Greensfunction to be used [STANDARD | INTEGRATED]",
                                               "INTEGRATED");
    
    itsAttr[BBOXINCR] = Attributes::makeReal("BBOXINCR",
                                             "Increase of bounding box in % ",
                                             2.0);

    // P3M only:
    itsAttr[RC] = Attributes::makeReal("RC",
                                       "cutoff radius for PP interactions",
                                       0.0);
    
    itsAttr[ALPHA] = Attributes::makeReal("ALPHA",
                                          "Green’s function splitting parameter",
                                          0.0);
    
    itsAttr[EPSILON] = Attributes::makeReal("EPSILON",
                                            "regularization for PP interaction",
                                            0.0);

    //SAAMG and in case of FFT with dirichlet BC in x and y
    itsAttr[GEOMETRY] = Attributes::makeString("GEOMETRY",
                                               "GEOMETRY to be used as domain boundary",
                                               "");
    
    itsAttr[ITSOLVER] = Attributes::makeString("ITSOLVER",
                                               "Type of iterative solver [CG | BiCGSTAB | GMRES]",
                                               "CG");
    
    itsAttr[INTERPL] = Attributes::makeString("INTERPL",
                                              "interpolation used for boundary points [CONSTANT | LINEAR | QUADRATIC]",
                                              "LINEAR");
    
    itsAttr[TOL] = Attributes::makeReal("TOL",
                                        "Tolerance for iterative solver",
                                        1e-8);
    
    itsAttr[MAXITERS] = Attributes::makeReal("MAXITERS",
                                             "Maximum number of iterations of iterative solver",
                                             100);
    
    itsAttr[PRECMODE] = Attributes::makeString("PRECMODE",
                                               "Preconditioner Mode [STD | HIERARCHY | REUSE]",
                                               "HIERARCHY");

    // AMR
#ifdef ENABLE_AMR
    itsAttr[AMR_MAXLEVEL] = Attributes::makeReal("AMR_MAXLEVEL",
                                                 "Maximum number of levels in AMR",
                                                 0);
    
    itsAttr[AMR_REFX] = Attributes::makeReal("AMR_REFX",
                                             "Refinement ration in x-direction in AMR",
                                             2);
    
    itsAttr[AMR_REFY] = Attributes::makeReal("AMR_REFY",
                                             "Refinement ration in y-direction in AMR",
                                             2);
    
    itsAttr[AMR_REFT] = Attributes::makeReal("AMR_REFT",
                                             "Refinement ration in z-direction in AMR",
                                             2);
    
    itsAttr[AMR_SUBCYCLE] = Attributes::makeBool("AMR_SUBCYCLE",
                                                 "Subcycling in time for refined levels in AMR",
                                                 false);
    
    itsAttr[AMR_MAXGRID] = Attributes::makeReal("AMR_MAXGRID",
                                                "Maximum grid size in AMR",
                                                16);
    
    itsAttr[AMR_TAGGING] = Attributes::makeString("AMR_TAGGING",
                                                  "Refinement criteria [CHARGE_DENSITY | POTENTIAL | EFIELD]",
                                                  "CHARGE_DENSITY");
    
    itsAttr[AMR_DENSITY] = Attributes::makeReal("AMR_DENSITY",
                                               "Tagging value for charge density refinement [C / cell volume]",
                                               1.0e-14);
    
    itsAttr[AMR_MAX_NUM_PART] = Attributes::makeReal("AMR_MAX_NUM_PART",
                                                     "Tagging value for max. #particles",
                                                     1);
    
    itsAttr[AMR_MIN_NUM_PART] = Attributes::makeReal("AMR_MIN_NUM_PART",
                                                     "Tagging value for min. #particles",
                                                     1);
    
    itsAttr[AMR_SCALING] = Attributes::makeReal("AMR_SCALING",
                                                "Scaling value for maximum value tagging "
                                                "(only POTENTIAL / CHARGE_DENSITY / "
                                                "MOMENTA", 0.75);
#endif
    
#ifdef HAVE_AMR_MG_SOLVER
    itsAttr[TAMG_SMOOTHER] = Attributes::makeString("TAMG_SMOOTHER",
                                                    "Smoothing of level solution", "GS");
    
    itsAttr[TAMG_NSWEEPS] = Attributes::makeReal("TAMG_NSWEEPS",
                                                 "Number of relaxation steps",
                                                 8);
    
    itsAttr[TAMG_PREC] = Attributes::makeString("TAMG_PREC",
                                                "Preconditioner of bottom solver",
                                                "NONE");
    
    itsAttr[TAMG_INTERP] = Attributes::makeString("TAMG_INTERP",
                                                  "Interpolater between levels",
                                                  "PC");
    
    itsAttr[TAMG_NORM] = Attributes::makeString("TAMG_NORM",
                                                "Norm for convergence criteria",
                                                "L2");
#endif

    mesh_m = 0;
    FL_m = 0;
    PL_m.reset(nullptr);

    solver_m = 0;

    registerOwnership(AttributeHandler::STATEMENT);
}


FieldSolver::FieldSolver(const std::string &name, FieldSolver *parent):
    Definition(name, parent)
{
    mesh_m = 0;
    FL_m = 0;
    PL_m.reset(nullptr);
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

std::string FieldSolver::getType() {
    return Attributes::getString(itsAttr[FSTYPE]);
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

    if(Util::toUpper(Attributes::getString(itsAttr[FSTYPE])) == "FFTPERIODIC") {
        decomp[0] = decomp[1] = SERIAL;
        decomp[2] = PARALLEL;
    }
    // create prototype mesh and layout objects for this problem domain
    mesh_m   = new Mesh_t(domain);
    FL_m     = new FieldLayout_t(*mesh_m, decomp);
    PL_m.reset(new Layout_t(*FL_m, *mesh_m));
    // OpalData::getInstance()->setMesh(mesh_m);
    // OpalData::getInstance()->setFieldLayout(FL_m);
    // OpalData::getInstance()->setLayout(PL_m);
}

bool FieldSolver::hasPeriodicZ() {
    return Util::toUpper(Attributes::getString(itsAttr[BCFFTT])) == std::string("PERIODIC");
}

#ifdef ENABLE_AMR
inline bool FieldSolver::isAmrSolverType() const {
    return Options::amr;
}

inline int FieldSolver::getAmrMaxLevel() const {
    return Attributes::getReal(itsAttr[AMR_MAXLEVEL]);
}

inline int FieldSolver::getAmrRefRatioX() const {
    return Attributes::getReal(itsAttr[AMR_REFX]);
}

inline int FieldSolver::getAmrRefRatioY() const {
    return Attributes::getReal(itsAttr[AMR_REFY]);
}

inline int FieldSolver::getAmrRefRatioT() const {
    return Attributes::getReal(itsAttr[AMR_REFT]);
}

inline bool FieldSolver::isAmrSubCycling() const {
    return Attributes::getBool(itsAttr[AMR_SUBCYCLE]);
}

inline int FieldSolver::getAmrMaxGridSize() const {
    return Attributes::getReal(itsAttr[AMR_MAXGRID]);
}
#endif

void FieldSolver::initSolver(PartBunchBase<double, 3> *b) {
    itsBunch_m = b;
    std::string bcx = Attributes::getString(itsAttr[BCFFTX]);
    std::string bcy = Attributes::getString(itsAttr[BCFFTY]);
    std::string bcz = Attributes::getString(itsAttr[BCFFTT]);
    
    #ifdef ENABLE_AMR
    if ( isAmrSolverType() ) {
        Inform m("FieldSolver::initAmrSolver");
        fsType_m = "AMR";
        
        initAmrObject_m();
        
        initAmrSolver_m();
        
    } else if(Attributes::getString(itsAttr[FSTYPE]) == "FFT") {
#else
    if(Attributes::getString(itsAttr[FSTYPE]) == "FFT") {
#endif
        bool sinTrafo = ((bcx == std::string("DIRICHLET")) && (bcy == std::string("DIRICHLET")) && (bcz == std::string("DIRICHLET")));
        if(sinTrafo) {
            std::cout << "FFTBOX ACTIVE" << std::endl;
            //we go over all geometries and add the Geometry Elements to the geometry list
            std::string geoms = Attributes::getString(itsAttr[GEOMETRY]);
            std::string tmp = "";
            //split and add all to list
            std::vector<BoundaryGeometry *> geometries;
            for(unsigned int i = 0; i <= geoms.length(); i++) {
                if(i == geoms.length() || geoms[i] == ',') {
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
    } else if (Util::toUpper(Attributes::getString(itsAttr[FSTYPE])) == "P3M") {
        solver_m = new P3MPoissonSolver(mesh_m,
                                        FL_m,
                                        Attributes::getReal(itsAttr[RC]),
                                        Attributes::getReal(itsAttr[ALPHA]),
                                        Attributes::getReal(itsAttr[EPSILON]));

        PL_m->setAllCacheDimensions(Attributes::getReal(itsAttr[RC]));
        PL_m->enableCaching();

        fsType_m = "P3M";
    } else if(Util::toUpper(Attributes::getString(itsAttr[FSTYPE])) == "SAAMG") {
#ifdef HAVE_SAAMG_SOLVER
        //we go over all geometries and add the Geometry Elements to the geometry list
        std::string geoms = Attributes::getString(itsAttr[GEOMETRY]);
        std::string tmp = "";
        //split and add all to list
        std::vector<BoundaryGeometry *> geometries;
        for(unsigned int i = 0; i <= geoms.length(); i++) {
            if(i == geoms.length() || geoms[i] == ',') {
                BoundaryGeometry *geom = OpalData::getInstance()->getGlobalGeometry();
                if(geom != 0) {
                    geometries.push_back(geom);
                }
                tmp.clear();
            } else
            tmp += geoms[i];
        }
        solver_m = new MGPoissonSolver(dynamic_cast<PartBunch*>(itsBunch_m), mesh_m, FL_m,
                                       geometries, Attributes::getString(itsAttr[ITSOLVER]),
                                       Attributes::getString(itsAttr[INTERPL]),
                                       Attributes::getReal(itsAttr[TOL]),
                                       Attributes::getReal(itsAttr[MAXITERS]),
                                       Attributes::getString(itsAttr[PRECMODE]));
        itsBunch_m->set_meshEnlargement(Attributes::getReal(itsAttr[BBOXINCR]) / 100.0);
        fsType_m = "SAAMG";
#else
        INFOMSG("SAAMG Solver not enabled! Please build OPAL with -DENABLE_SAAMG_SOLVER=1" << endl);
        INFOMSG("switching to FFT solver..." << endl);
        solver_m = new FFTPoissonSolver(mesh_m, FL_m, Util::toUpper(Attributes::getString(itsAttr[GREENSF])),bcz);
        fsType_m = "FFT";
#endif
    } else {
        solver_m = 0;
        INFOMSG("no solver attached" << endl);
    }
}

bool FieldSolver::hasValidSolver() {
    return (solver_m != 0);
}

Inform &FieldSolver::printInfo(Inform &os) const {
    std::string fsType = Util::toUpper(Attributes::getString(itsAttr[FSTYPE]));

    os << "* ************* F I E L D S O L V E R ********************************************** " << endl;
    os << "* FIELDSOLVER  " << getOpalName() << '\n'
       << "* TYPE         " << fsType << '\n'
       << "* N-PROCESSORS " << Ippl::getNodes() << '\n'
       << "* MX           " << Attributes::getReal(itsAttr[MX])   << '\n'
       << "* MY           " << Attributes::getReal(itsAttr[MY])   << '\n'
       << "* MT           " << Attributes::getReal(itsAttr[MT])   << '\n'
       << "* BBOXINCR     " << Attributes::getReal(itsAttr[BBOXINCR]) << endl;

    if(fsType == "P3M")
        os << "* RC           " << Attributes::getReal(itsAttr[RC]) << '\n'
           << "* ALPHA        " << Attributes::getReal(itsAttr[ALPHA]) << '\n'
           << "* EPSILON      " << Attributes::getReal(itsAttr[EPSILON]) << endl;


    if(fsType == "FFT") {
        os << "* GRRENSF      " << Util::toUpper(Attributes::getString(itsAttr[GREENSF])) << endl;
    } else if (fsType == "SAAMG") {
        os << "* GEOMETRY     " << Attributes::getString(itsAttr[GEOMETRY]) << '\n'
           << "* ITSOLVER     " << Util::toUpper(Attributes::getString(itsAttr[ITSOLVER]))   << '\n'
           << "* INTERPL      " << Util::toUpper(Attributes::getString(itsAttr[INTERPL]))  << '\n'
           << "* TOL          " << Attributes::getReal(itsAttr[TOL])        << '\n'
           << "* MAXITERS     " << Attributes::getReal(itsAttr[MAXITERS]) << '\n'
           << "* PRECMODE     " << Attributes::getString(itsAttr[PRECMODE])   << endl;
    }
#ifdef ENABLE_AMR
    else if (fsType == "AMR" || Options::amr) {
        os << "* AMR_MAXLEVEL     " << Attributes::getReal(itsAttr[AMR_MAXLEVEL]) << '\n'
           << "* AMR_REFX         " << Attributes::getReal(itsAttr[AMR_REFX]) << '\n'
           << "* AMR_REFY         " << Attributes::getReal(itsAttr[AMR_REFY]) << '\n'
           << "* AMR_REFT         " << Attributes::getReal(itsAttr[AMR_REFT]) << '\n'
           << "* AMR_SUBCYCLE     " << Attributes::getBool(itsAttr[AMR_SUBCYCLE]) << '\n'
           << "* AMR_MAXGRID      " << Attributes::getReal(itsAttr[AMR_MAXGRID]) << '\n'
           << "* AMR_TAGGING      " << Attributes::getString(itsAttr[AMR_TAGGING]) <<'\n'
           << "* AMR_DENSITY      " << Attributes::getReal(itsAttr[AMR_DENSITY]) << '\n'
           << "* AMR_MAX_NUM_PART " << Attributes::getReal(itsAttr[AMR_MAX_NUM_PART]) << '\n'
           << "* AMR_MIN_NUM_PART " << Attributes::getReal(itsAttr[AMR_MIN_NUM_PART]) << '\n'
           << "* AMR_DENSITY      " << Attributes::getReal(itsAttr[AMR_DENSITY]) << '\n'
           << "* AMR_SCALING      " << Attributes::getReal(itsAttr[AMR_SCALING]) << endl;
    }
#endif

#ifdef HAVE_AMR_MG_SOLVER
    if (fsType == "TAMG") {
        os << "* ITSOLVER (TAMG)  "
           << Util::toUpper(Attributes::getString(itsAttr[ITSOLVER])) << '\n'
           << "* TAMG_PREC        "
           << Util::toUpper(Attributes::getString(itsAttr[TAMG_PREC])) << '\n'
           << "* TAMG_SMOOTHER    "
           << Util::toUpper(Attributes::getString(itsAttr[TAMG_SMOOTHER])) << '\n'
           << "* TAMG_NSWEEPS     "
           << Attributes::getReal(itsAttr[TAMG_NSWEEPS]) << '\n'
           << "* TAMG_INTERP      "
           << Util::toUpper(Attributes::getString(itsAttr[TAMG_INTERP])) << '\n'
           << "* TAMG_NORM        "
           << Util::toUpper(Attributes::getString(itsAttr[TAMG_NORM])) << '\n'
           << "* BCFFTX           "
           << Util::toUpper(Attributes::getString(itsAttr[BCFFTX])) << '\n'
           << "* BCFFTY           "
           << Util::toUpper(Attributes::getString(itsAttr[BCFFTY])) << '\n'
           << "* BCFFTT           "
           << Util::toUpper(Attributes::getString(itsAttr[BCFFTT])) << endl;
    }
#endif

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

#ifdef ENABLE_AMR
void FieldSolver::initAmrObject_m() {
    
    itsBunch_m->set_meshEnlargement(Attributes::getReal(itsAttr[BBOXINCR]) * 0.01);
    
    // setup initial info for creating the object
    AmrBoxLib::AmrInitialInfo info;
    info.gridx      = (int)this->getMX();
    info.gridy      = (int)this->getMY();
    info.gridz      = (int)this->getMT();
    info.maxgrid    = this->getAmrMaxGridSize();
    info.maxlevel   = this->getAmrMaxLevel();
    info.refratx    = this->getAmrRefRatioX();
    info.refraty    = this->getAmrRefRatioY();
    info.refratz    = this->getAmrRefRatioT();
    
    
    itsAmrObject_mp = AmrBoxLib::create(info, dynamic_cast<AmrPartBunch*>(itsBunch_m));
    
    itsAmrObject_mp->setTagging( Attributes::getString(itsAttr[AMR_TAGGING]) );
    
    itsAmrObject_mp->setScalingFactor( Attributes::getReal(itsAttr[AMR_SCALING]) );
    
    itsAmrObject_mp->setChargeDensity( Attributes::getReal(itsAttr[AMR_DENSITY]) );
    
    itsAmrObject_mp->setMaxNumParticles(
        Attributes::getReal(itsAttr[AMR_MAX_NUM_PART])
    );
    
    itsAmrObject_mp->setMinNumParticles(
        Attributes::getReal(itsAttr[AMR_MIN_NUM_PART])
    );
}


void FieldSolver::initAmrSolver_m() {
    std::string fsType = Util::toUpper(Attributes::getString(itsAttr[FSTYPE]));
    if (fsType == "FMG") {
        
        if ( dynamic_cast<AmrBoxLib*>( itsAmrObject_mp.get() ) == 0 )
            throw OpalException("FieldSolver::initAmrSolver_m()",
                                "FMultiGrid solver requires BoxLib.");
        
        solver_m = new FMGPoissonSolver(static_cast<AmrBoxLib*>(itsAmrObject_mp.get()));
        
    } else if (fsType == "HYPRE") {
        throw OpalException("FieldSolver::initAmrSolver_m()",
                            "HYPRE solver not yet implemented.");
    } else if (fsType == "HPGMG") {
        throw OpalException("FieldSolver::initAmrSolver_m()",
                            "HPGMG solver not yet implemented.");
    } else if (fsType == "TAMG") {
        if ( dynamic_cast<AmrBoxLib*>( itsAmrObject_mp.get() ) == 0 )
            throw OpalException("FieldSolver::initAmrSolver_m()",
                                "FMultiGrid solver requires BoxLib.");
        
        solver_m = new AmrMultiGrid(static_cast<AmrBoxLib*>(itsAmrObject_mp.get()),
                                    Attributes::getString(itsAttr[ITSOLVER]),
                                    Attributes::getString(itsAttr[TAMG_PREC]),
                                    Attributes::getString(itsAttr[BCFFTX]),
                                    Attributes::getString(itsAttr[BCFFTY]),
                                    Attributes::getString(itsAttr[BCFFTT]),
                                    Attributes::getString(itsAttr[TAMG_SMOOTHER]),
                                    Attributes::getReal(itsAttr[TAMG_NSWEEPS]),
                                    Attributes::getString(itsAttr[TAMG_INTERP]),
                                    Attributes::getString(itsAttr[TAMG_NORM]));
    } else
        throw OpalException("FieldSolver::initAmrSolver_m()",
                            "Unknown solver " + fsType + ".");
}
#endif
