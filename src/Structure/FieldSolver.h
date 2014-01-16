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
#include "Solvers/PoissonSolver.h"

class PartBunch;


// Class FieldSolver
// ------------------------------------------------------------------------
/// The FieldSolver definition.
//  A FieldSolver definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

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

    bool hasPeriodicZ();

    bool isAMRSolver();

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
