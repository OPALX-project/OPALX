#ifndef LF2_H
#define LF2_H

#include "Stepper.h"
#include "Physics/Physics.h"

/// Leap-Frog 2nd order
template <typename FieldFunction, typename ... Arguments>
class LF2 : public Stepper<FieldFunction, Arguments...> {
    
public:
    
    LF2(const FieldFunction& fieldfunc) : Stepper<FieldFunction, Arguments ...>(fieldfunc) { }
    
    bool advance(PartBunchBase<double, 3>* bunch,
                 const size_t& i,
                 const double& t,
                 const double dt,
                 Arguments& ... args) const;
    
private:
    
    void push_m(Vector_t& R, const Vector_t& P, const double& h) const;
    
    bool kick_m(PartBunchBase<double, 3>* bunch, const size_t& i,
                const double& t, const double& h,
                Arguments& ... args) const;
};

#include "LF2.hpp"

#endif
