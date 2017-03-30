#ifndef AMR_PART_BUNCH_H
#define AMR_PART_BUNCH_H

#include "Algorithms/PartBunch.h"
#include "Amr/BoxLibParticle.h"
#include "Amr/BoxLibLayout.h"

class AmrPartBunch : public PartBunch,
                     public BoxLibParticle<BoxLibLayout<double, 3> >
{
public:
    
    AmrPartBunch(const PartData *ref);

    /// Conversion.
    AmrPartBunch(const std::vector<OpalParticle> &, const PartData *ref);

    AmrPartBunch(const AmrPartBunch &);
    
    void computeSelfFields();
    
    void computeSelfFields(int b);
    
    void computeSelfFields_cycl(double gamma);
    
    void computeSelfFields_cycl(int b);
    
    void update();
    
    void update(const ParticleAttrib<char>& canSwap);
};

#endif
