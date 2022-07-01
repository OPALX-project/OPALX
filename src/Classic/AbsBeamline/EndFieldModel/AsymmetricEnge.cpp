#include "AbsBeamline/EndFieldModel/AsymmetricEnge.h"

namespace endfieldmodel {

AsymmetricEnge::AsymmetricEnge() : engeStart_m(new Enge()),
    engeEnd_m(new Enge()), x0Start_m(0.0), x0End_m(0.0) {
}

AsymmetricEnge::AsymmetricEnge(const AsymmetricEnge& rhs) 
  : engeStart_m(rhs.engeStart_m->clone()), engeEnd_m(rhs.engeEnd_m->clone()),
    x0Start_m(rhs.x0Start_m), x0End_m(rhs.x0End_m) {
}

AsymmetricEnge::AsymmetricEnge(const std::vector<double> aStart,
                       double x0Start,
                       double lambdaStart, 
                       const std::vector<double> aEnd,
                       double x0End,
                       double lambdaEnd) : engeStart_m(new Enge()),
    engeEnd_m(new Enge()), x0Start_m(x0Start), x0End_m(x0End) {
        engeStart_m->setCoefficients(aStart);
        engeStart_m->setLambda(lambdaStart);
        // x0 is held in this
        engeEnd_m->setCoefficients(aEnd);
        engeEnd_m->setLambda(lambdaEnd);
}

void AsymmetricEnge::rescale(double scaleFactor) {
    engeStart_m->rescale(scaleFactor);
    engeEnd_m->rescale(scaleFactor);
    x0Start_m *= scaleFactor;
    x0End_m *= scaleFactor;
}


}
