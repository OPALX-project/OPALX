#ifndef AMR_PART_BUNCH_H
#define AMR_PART_BUNCH_H

#include "PartBunchBase.h"

class AmrPartBunch : public PartBunchBase<double, 3>
{
public:
    
    AmrPartBunch(AbstractParticle<double, 3>* pb) : PartBunchBase<double, 3>(pb) {
//         this->initialize(new BoxLibLayout<double, 3>()); // TODO Where is this done?
//         this->addAttribute(Bin);
    }
    
};

#endif
