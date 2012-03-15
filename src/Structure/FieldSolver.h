#ifndef OPAL_FieldSolver_HH
#define OPAL_FieldSolver_HH

// ------------------------------------------------------------------------
// $RCSfile: FieldSolver.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FieldSolver
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

class FieldSolver;
#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
//class PartBunch;
#include "Algorithms/PartBunch.h"
//#include "Ippl.h"
//class FFTPosissonSolver;
#include "Solvers/PoissonSolver.h"
#include "Solvers/FFTPoissonSolver.h"
#include "Solvers/FFTBoxPoissonSolver.h"
#ifdef HAVE_ML_SOLVER
#include "Solvers/MGPoissonSolver.h"
#endif


// Class FieldSolver
// ------------------------------------------------------------------------
/// The FieldSolver definition.
//  A FieldSolver definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

// typedef IntCIC  IntrplCIC_t;
// typedef IntNGP  IntrplNGP_t;
// typedef IntSUDS IntrplSUDS_t;

// typedef ParticleSpatialLayout<double, 3>::ParticlePos_t Ppos_t;
// typedef ParticleSpatialLayout<double, 3>::ParticleIndex_t PID_t;

// typedef ParticleAttrib<double> Pscalar_t;

// typedef InterpolatorTraits<double, 3, IntrplCIC_t>::Cache_t Pcache_t;

// typedef UniformCartesian<3, double> Mesh_t;

// typedef ParticleSpatialLayout<double, 3>::SingleParticlePos_t Vector_t;

// typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;

// typedef Cell                                       Center_t;

// typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
// typedef Field<double, 3, Mesh_t, Center_t>       Field_t;
// typedef Field<Vector_t, 3, Mesh_t, Center_t>     VField_t;

class FieldSolver: public Definition {

public:

    /// Exemplar constructor.
    FieldSolver();

    virtual ~FieldSolver();

    /// Make clone.
    virtual FieldSolver *clone(const string &name);

    /// Find named FieldSolver.
    static FieldSolver *find(const string &name);

    /// Return meshsize
    double getMX() const;

    /// Return meshsize
    double getMY() const;

    /// Return meshsize
    double getMT() const;

    /// Store emittance for mode 1.
    void setMX(double);

    /// Store emittance for mode 2.
    void setMY(double);

    /// Store emittance for mode 3.
    void setMT(double);

    /// Update the field solver data.
    virtual void update();

    /// Execute (init) the field solver data.
    virtual void execute();

    void initCartesianFields();

    void initSolver(PartBunch &b);

    bool hasValidSolver();

    std::string getFieldSolverType() {return fsType_m; }

    inline Layout_t &getParticleLayout() { return *PL_m; }

    Inform &printInfo(Inform &os) const;
    unsigned int getInteractionRadius() {return (unsigned int) rpp_m; }

    /// the actual solver, should be a base object
    PoissonSolver *solver_m;

private:

    // Not implemented.
    FieldSolver(const FieldSolver &);
    void operator=(const FieldSolver &);

    // Clone constructor.
    FieldSolver(const string &name, FieldSolver *parent);

    /// The cartesian mesh
    Mesh_t *mesh_m;

    /// The field layout f
    FieldLayout_t *FL_m;

    /// The particle layout
    Layout_t *PL_m;

    /// all the particles are here ...
    PartBunch *itsBunch_m;

    std::string fsType_m;

    double rpp_m;

};

inline Inform &operator<<(Inform &os, const FieldSolver &fs) {
    return fs.printInfo(os);
}

#endif // OPAL_FieldSolver_HH
