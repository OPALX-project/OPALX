
#ifndef POISSON_SOLVER_H_
#define POISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "Algorithms/PBunchDefs.h"

#ifdef HAVE_AMR_SOLVER
    #include "Utilities/OpalException.h"
#endif
//////////////////////////////////////////////////////////////
class PartBunch;
//use Barton and Nackman Trick to avoid virtual functions
//template <class T_Leaftype>
class PoissonSolver {
public:

    //T_Leaftype& asLeaf() { return static_cast<T_Leaftype&>(*this); }

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential in open space
    //void computePotential(Field_t &rho, Vector_t hr) {asLeaf().computePotential(rho, hr);}
    //void computePotential(Field_t &rho, Vector_t hr, double zshift) {asLeaf().computePotential(&rho, hr, zshift);}

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential in open space
    virtual void computePotential(Field_t &rho, Vector_t hr) = 0;
    
    /*!
     * @param rho specifies the charge-density field
     * @param baseLevel of adaptive mesh refinement solvers (AMR). Used in case of sub-cycling in time.
     * @param finestLevel of AMR.
     */
#ifdef HAVE_AMR_SOLVER
    virtual void solve(AmrFieldContainer_t &rho,
                       AmrFieldContainer_t &efield,
                       unsigned short baseLevel,
                       unsigned short finestLevel)
    {
        throw OpalException("PoissonSolver::solve()", "Not supported for non-AMR code.");
    };
#endif
                                  
    virtual void computePotential(Field_t &rho, Vector_t hr, double zshift) = 0;

    virtual double getXRangeMin(unsigned short level = 0) = 0;
    virtual double getXRangeMax(unsigned short level = 0) = 0;
    virtual double getYRangeMin(unsigned short level = 0) = 0;
    virtual double getYRangeMax(unsigned short level = 0) = 0;
    virtual double getZRangeMin(unsigned short level = 0) = 0;
    virtual double getZRangeMax(unsigned short level = 0) = 0;
    virtual void test(PartBunch &bunch) = 0 ;
    virtual ~PoissonSolver(){};

};

inline Inform &operator<<(Inform &os, const PoissonSolver &fs) {
    return os << "";
}

#endif
