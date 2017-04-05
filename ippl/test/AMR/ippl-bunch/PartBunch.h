#ifndef PART_BUNCH_H
#define PART_BUNCH_H

#include "PartBunchBase.h"

class PartBunch : public PartBunchBase<double, 3>
{
public:
    
//     PartBunch(IpplParticleBase<ParticleSpatialLayout<double, 3> >* pb) : PartBunchBase<double, 3>(pb) {
    PartBunch(AbstractParticle<double, 3>* pb) : PartBunchBase<double, 3>(pb) {
//         dynamic_cast<IpplParticleBase<ParticleSpatialLayout<double, 3> >* >(pb)->initialize( ); // TODO Where is this done?
//         this->addAttribute(Bin);
        
    }
    
};

#endif
