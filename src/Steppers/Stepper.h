#ifndef STEPPER_H
#define STEPPER_H

#include "Algorithms/PartBunch.h"
#include "Algorithms/Vektor.h"

#include <functional>

/*!
 * @precondition The field function has to return a
 * boolean and take at least the following arguments in that order:
 *  - double    specifying the time
 *  - int       specifying the i-th particle
 *  - Vector_t  specifying the electric field
 *  - Vector_t  specifying the magnetic field
 */

template <typename FieldFunction, typename ... Arguments>
class Stepper {
    
public:
    
    Stepper(const FieldFunction& fieldfunc) : fieldfunc_m(fieldfunc) { }
    
    virtual bool advance(PartBunch* bunch,
                         const size_t& i,
                         const double& t,
                         const double dt,
                         Arguments& ... args) const = 0;

protected:
    const FieldFunction& fieldfunc_m;
};

#endif
