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

template <class T, unsigned Dim>
class PartBunchBase;


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
    virtual FieldSolver *clone(const std::string &name);

    /// Find named FieldSolver.
    static FieldSolver *find(const std::string &name);

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

    void initSolver(PartBunchBase<double, 3> *b);

    bool hasValidSolver();

    std::string getFieldSolverType() {return fsType_m; }

    inline Layout_t &getParticleLayout() { return *PL_m; }
    
    FieldLayout_t *getFieldLayout() { return FL_m; }
    
    Inform &printInfo(Inform &os) const;
    unsigned int getInteractionRadius() {return (unsigned int) rpp_m; }

    bool hasPeriodicZ();

#ifdef HAVE_AMR_SOLVER
    bool isAmrSolver();

    int getAmrMaxLevel();

    int getAmrRefRatioX();

    int getAmrRefRatioY();

    int getAmrRefRatioT();

    bool isAmrSubCycling();

    int getAmrMaxGridSize();
#endif

    /// the actual solver, should be a base object
    PoissonSolver *solver_m;

private:

    // Not implemented.
    FieldSolver(const FieldSolver &);
    void operator=(const FieldSolver &);

    // Clone constructor.
    FieldSolver(const std::string &name, FieldSolver *parent);

    /// The cartesian mesh
    Mesh_t *mesh_m;

    /// The field layout f
    FieldLayout_t *FL_m;

    /// The particle layout
    Layout_t *PL_m;

    /// all the particles are here ...
    PartBunchBase<double, 3> *itsBunch_m;

    std::string fsType_m;

    double rpp_m;

};

inline Inform &operator<<(Inform &os, const FieldSolver &fs) {
    return fs.printInfo(os);
}

#endif // OPAL_FieldSolver_HH