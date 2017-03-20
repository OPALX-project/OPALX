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
#include "Algorithms/PartBunch.h"

#ifdef HAVE_AMR_SOLVER
    #include <ParmParse.H>
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
#ifdef HAVE_AMR_SOLVER
        AMRMAXLEVEL, // AMR, maximum refinement level
        AMRREFX,     // AMR, refinement ratio in x
        AMRREFY,     // AMR, refinement ratio in y
        AMRREFT,     // AMR, refinement ration in z
        AMRSUBCYCLE, // AMR, subcycling in time for refined levels (default: false)
        AMRMAXGRID,  // AMR, maximum grid size (default: 16)
#endif
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
    itsAttr[BCFFTX] = Attributes::makeString("BCFFTX", "Boundary conditions in x: open, dirichlet (box), periodic");
    itsAttr[BCFFTY] = Attributes::makeString("BCFFTY", "Boundary conditions in y: open, dirichlet (box), periodic");
    itsAttr[BCFFTT] = Attributes::makeString("BCFFTT", "Boundary conditions in z(t): open, periodic");

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

    // AMR
#ifdef HAVE_AMR_SOLVER
    itsAttr[AMRMAXLEVEL] = Attributes::makeReal("AMRMAXLEVEL", "Maximum number of levels in AMR", 0);
    itsAttr[AMRREFX] = Attributes::makeReal("AMRREFX", "Refinement ration in x-direction in AMR", 2);
    itsAttr[AMRREFY] = Attributes::makeReal("AMRREFY", "Refinement ration in y-direction in AMR", 2);
    itsAttr[AMRREFT] = Attributes::makeReal("AMRREFT", "Refinement ration in z-direction in AMR", 2);
    itsAttr[AMRSUBCYCLE] = Attributes::makeBool("AMRSUBCYCLE",
                                                "Subcycling in time for refined levels in AMR", false);
    itsAttr[AMRMAXGRID] = Attributes::makeReal("AMRMAXGRID", "Maximum grid size in AMR", 16);

#endif

    mesh_m = 0;
    FL_m = 0;
    PL_m = 0;

    solver_m = 0;

    registerOwnership(AttributeHandler::STATEMENT);
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

    if(Util::toUpper(Attributes::getString(itsAttr[FSTYPE])) == "FFTPERIODIC") {
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
    return Util::toUpper(Attributes::getString(itsAttr[BCFFTT])) == std::string("PERIODIC");
}

#ifdef HAVE_AMR_SOLVER
bool FieldSolver::isAMRSolver() {
    return Util::toUpper(Attributes::getString(itsAttr[FSTYPE])) == std::string("AMR");
}

int FieldSolver::amrMaxLevel() {
    return Attributes::getReal(itsAttr[AMRMAXLEVEL]);
}

int FieldSolver::amrRefRatioX() {
    return Attributes::getReal(itsAttr[AMRREFX]);
}

int FieldSolver::amrRefRatioY() {
    return Attributes::getReal(itsAttr[AMRREFY]);
}

int FieldSolver::amrRefRatioT() {
    return Attributes::getReal(itsAttr[AMRREFT]);
}

bool FieldSolver::amrSubCycling() {
    return Attributes::getBool(itsAttr[AMRSUBCYCLE]);
}

int FieldSolver::amrMaxGridSize() {
    return Attributes::getReal(itsAttr[AMRMAXGRID]);
}
#endif

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
        solver_m = new FFTPoissonSolver(mesh_m, FL_m, Util::toUpper(Attributes::getString(itsAttr[GREENSF])),bcz);
        fsType_m = "FFT";
#endif
    }
#ifdef HAVE_AMR_SOLVER
    else if (Attributes::getString(itsAttr[FSTYPE]) == "AMR") {
        Inform m("FieldSolver::initSolver-amr ");
        fsType_m = "AMR";

	/*
        // Add the parsed AMR attributes to BoxLib (please check BoxLib/Src/C_AMRLib/Amr.cpp)
        ParmParse pp("amr");

        pp.add("max_level", Attributes::getReal(itsAttr[AMRMAXLEVEL]));

        pp.add("ref_ratio", (int)Attributes::getReal(itsAttr[AMRREFX])); //FIXME

        pp.add("max_grid_size", Attributes::getReal(itsAttr[AMRMAXGRID]));

        FieldLayout<3>::iterator_iv locDomBegin = FL_m->begin_iv();
        FieldLayout<3>::iterator_iv locDomEnd = FL_m->end_iv();
        FieldLayout<3>::iterator_dv globDomBegin = FL_m->begin_rdv();
        FieldLayout<3>::iterator_dv globDomEnd = FL_m->end_rdv();

        BoxArray lev0_grids(Ippl::getNodes());

        Array<int> procMap;
        procMap.resize(lev0_grids.size()+1); // +1 is a historical thing, do not ask

        // first iterate over the local owned domain(s)
        for(FieldLayout<3>::const_iterator_iv v_i = locDomBegin ; v_i != locDomEnd; ++v_i) {
            std::ostringstream stream;
            stream << *((*v_i).second);

            std::pair<Box,unsigned int> res = getBlGrids(stream.str());
            lev0_grids.set(res.second,res.first);
            procMap[res.second] = Ippl::myNode();
        }

        // then iterate over the non-local domain(s)
        for(FieldLayout<3>::iterator_dv v_i = globDomBegin ; v_i != globDomEnd; ++v_i) {
            std::ostringstream stream;
            stream << *((*v_i).second);

            std::pair<Box,unsigned int> res = getBlGrids(stream.str());
            lev0_grids.set(res.second,res.first);
            procMap[res.second] = res.second;
        }
        procMap[lev0_grids.size()] = Ippl::myNode();

        // This init call will cache the distribution map as determined by procMap
        // so that all data will end up on the right processor
        RealBox rb;
        Array<Real> prob_lo(3);
        Array<Real> prob_hi(3);

        prob_lo[0] = -0.02; //-0.08);
        prob_lo[1] = -0.02; //-0.08);
        prob_lo[2] =  0.0; //-0.12);
        prob_hi[0] =  0.02; //0.08);
        prob_hi[1] =  0.02; //0.08);
        prob_hi[2] =  0.04; //0.16);

        rb.setLo(prob_lo);
        rb.setHi(prob_hi);

        int coord_sys = 0;

        NDIndex<3> ipplDom = FL_m->getDomain();

        Array<int> ncell(3);
        ncell[0] = ipplDom[0].length();
        ncell[1] = ipplDom[1].length();
        ncell[2] = ipplDom[2].length();

        std::vector<int   > nr(3);
        std::vector<double> hr(3);
        std::vector<double> prob_lo_in(3);
        for (int i = 0; i < 3; i++) {
            nr[i] = ncell[i];
            hr[i] = (prob_hi[i] - prob_lo[i]) / ncell[i];
            prob_lo_in[i] = prob_lo[i];
        }

        int maxLevel = -1;
        amrptr_m = new Amr(&rb,maxLevel,ncell,coord_sys);

        if(amrptr_m)
            m << fsType_m << " solver: amrptr_m ready " << endl;

        Real strt_time = 0.0;
        Real stop_time = 1.0;

        amrptr_m->InitializeInit(strt_time, stop_time, &lev0_grids, &procMap);

        if(amrptr_m)
            m << fsType_m << " solver: amrptr_m Init done " << endl;
	*/
    }
#endif
    else {
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
#ifdef HAVE_AMR_SOLVER
    else if (fsType == "AMR") {
        os << "* AMRMAXLEVEL  " << Attributes::getReal(itsAttr[AMRMAXLEVEL]) << '\n'
           << "* AMRREFX      " << Attributes::getReal(itsAttr[AMRREFX]) << '\n'
           << "* AMRREFY      " << Attributes::getReal(itsAttr[AMRREFY]) << '\n'
           << "* AMRREFT      " << Attributes::getReal(itsAttr[AMRREFT]) << '\n'
           << "* AMRSUBCYCLE  " << Attributes::getBool(itsAttr[AMRSUBCYCLE]) << '\n'
           << "* AMRMAXGRID   " << Attributes::getReal(itsAttr[AMRMAXGRID]) << endl;
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
