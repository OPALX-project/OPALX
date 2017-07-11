#ifndef RK4_H
#define RK4_H

#include "Stepper.h"
#include "Physics/Physics.h"

/// 4-th order Runnge-Kutta stepper
template <typename FieldFunction, typename ... Arguments>
class RK4 : public Stepper<FieldFunction, Arguments...> {
    
public:
    
    RK4(const FieldFunction& fieldfunc) : Stepper<FieldFunction, Arguments ...>(fieldfunc) { }

    bool advance(PartBunch* bunch,
                 const size_t& i,
                 const double& t,
                 const double dt,
                 Arguments& ... args) const;
    
private:
    /**
     * 
     *
     * @param y
     * @param t
     * @param yp
     * @param Pindex
     *
     * @return
     */
    bool derivate_m(PartBunch* bunch,
                    double *y,
                    const double& t,
                    double* yp,
                    const size_t& i,
                    Arguments& ... args) const;
    
    
    void copyTo(const Vector_t& R, const Vector_t& P, double* x) const;
    
    void copyFrom(Vector_t& R, Vector_t& P, double* x) const;
    
    const double c_mmtns = Physics::c * 1.0e-6; // m/s --> mm/ns
    const double c_mtns = Physics::c * 1.0e-9;  // m/s --> m/ns
    const double mass_coeff = 1.0e18 * Physics::q_e / Physics::c / Physics::c; // from GeV/c^2 to basic unit: GV*C*s^2/m^2
};

#include "RK4.hpp"

#endif
