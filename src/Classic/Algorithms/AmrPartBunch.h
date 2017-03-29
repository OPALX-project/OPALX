#ifndef AMR_PART_BUNCH_H
#define AMR_PART_BUNCH_H

#include "Algorithms/PartBunchBase.h"
#include "Amr/BoxLibParticle.h"
#include "Amr/BoxLibLayout.h"

class AmrPartBunch : public PartBunchBase,
                     public BoxLibParticle<BoxLibLayout<double, 3> >
{
public:
    void computeSelfFields();
    
    void computeSelfFields(int b);
    
    void computeSelfFields_cycl(double gamma);
    
    void computeSelfFields_cycl(int b);
    
};

#endif
