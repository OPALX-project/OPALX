
#ifndef POISSON_SOLVER_H_
#define POISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "Algorithms/PBunchDefs.h"

#ifdef ENABLE_AMR
    #include "Utilities/OpalException.h"
#endif
//////////////////////////////////////////////////////////////
template <class T, unsigned Dim>
class PartBunchBase;
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
    
#ifdef ENABLE_AMR
    /**
     * AMR solver calls
     * 
     * @param rho right-hand side charge density on grid [C / m]
     * @param phi electrostatic potential (unknown) [V]
     * @param efield electric field [V / m]
     * @param baseLevel for solve
     * @param finestLevel for solve
     * @param prevAsGuess use of previous solution as initial guess
     */
    virtual void solve(AmrFieldContainer_t &rho,
                       AmrFieldContainer_t &phi,
                       AmrFieldContainer_t &efield,
                       unsigned short baseLevel,
                       unsigned short finestLevel,
                       bool prevAsGuess = true)
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
    virtual void test(PartBunchBase<double, 3> *bunch) = 0 ;
    virtual ~PoissonSolver(){};

};

inline Inform &operator<<(Inform &os, const PoissonSolver &fs) {
    return os << "";
}

#endif
