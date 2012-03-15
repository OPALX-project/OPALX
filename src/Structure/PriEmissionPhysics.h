#ifndef OPAL_PRIEMISSION_PHYSICS_HH
#define OPAL_PRIEMISSION_PHYSICS_HH

class OpalBeamline;
class ElementBase;

#include "Algorithms/PartBunch.h"
#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "Ippl.h"
class PriEmissionPhysics {

public:
    PriEmissionPhysics();

    ~PriEmissionPhysics();

    void Fieldemission(PartBunch *itsBunch, const double &fa, const double &Enormal, const double &parameterFNB,
                       const double &workFunction, const double &parameterFNVYZe, const double &parameterFNVYSe,
                       const double &parameterFNY, const double &fieldEnhancement, const double &maxFNemission,
                       const double &TriArea, const vector<Vector_t> &vertex, const Vector_t TriNormal, size_t &Nstp);

};




#endif // OPAL_PRIEMISSION_PHYSICS_HH
