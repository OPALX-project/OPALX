#ifndef PART_BUNCH_H
#define PART_BUNCH_H

#include "PartBunchBase.h"

class PartBunch : public PartBunchBase<double, 3>
{
public:
    
//     PartBunch(IpplParticleBase<ParticleSpatialLayout<double, 3> >* pb) : PartBunchBase<double, 3>(pb) {
    PartBunch(AbstractParticle<double, 3>* pb) : PartBunchBase<double, 3>(pb),
                                                 mesh_m(), fieldlayout_m(mesh_m) //FIXME
    {
//         dynamic_cast<IpplParticleBase<ParticleSpatialLayout<double, 3> >* >(pb)->initialize( ); // TODO Where is this done?
//         this->addAttribute(Bin);
    }
    
    const Mesh_t &getMesh() const {
        return mesh_m;
    }

    Mesh_t &getMesh() {
        return mesh_m;
    }

    FieldLayout_t &getFieldLayout() {
        return fieldlayout_m;
    }
    
private:
    Mesh_t mesh_m;
    FieldLayout_t fieldlayout_m;
    
};

#endif
