#ifndef AMR_PART_BUNCH_H
#define AMR_PART_BUNCH_H

#include <type_traits>

#include "PartBunchBase.h"

class AmrPartBunch : public PartBunchBase<double, 3>
{
public:
    
    AmrPartBunch(AbstractParticle<double, 3>* pb) : PartBunchBase<double, 3>(pb),
                                                    mesh_m(), fieldlayout_m(mesh_m) //FIXME
    {
//         this->initialize(new BoxLibLayout<double, 3>()); // TODO Where is this done?
//         this->addAttribute(Bin);
    }
    
    
    Mesh_t &getMesh() {
        return mesh_m;
    }

    FieldLayout_t &getFieldLayout() {
        return fieldlayout_m;
    }
    
private:
    static_assert(std::is_same<UniformCartesian<3, double>, Mesh_t>::value,
                  "Mesh_t has to be of type UniformCartesian<3, double>");
    
    
    Mesh_t mesh_m;
    FieldLayout_t fieldlayout_m;
    
};

#endif
