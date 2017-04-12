#ifndef PART_BUNCH_H
#define PART_BUNCH_H

#include "PartBunchBase.h"

class PartBunch : public PartBunchBase<double, 3>
{
public:
    typedef IpplParticleBase<ParticleSpatialLayout<double, 3> > pbase_t;
public:
    
//     PartBunch(IpplParticleBase<ParticleSpatialLayout<double, 3> >* pb) : PartBunchBase<double, 3>(pb) {
    PartBunch(AbstractParticle<double, 3>* pb) : PartBunchBase<double, 3>(pb)/*,
                                                 mesh_m(), fieldlayout_m(mesh_m) //FIXME*/
    {
//         dynamic_cast<IpplParticleBase<ParticleSpatialLayout<double, 3> >* >(pb)->initialize( ); // TODO Where is this done?
//         this->addAttribute(Bin);
    }
    
//     PartBunch() : PartBunchBase() { }
    
//     const Mesh_t &getMesh() const {
//         return mesh_m;
//     }

//     Mesh_t &getMesh() {
//         return mesh_m;
//     }

//     FieldLayout_t &getFieldLayout() {
//         return fieldlayout_m;
//     }

    void setFieldLayout(FieldLayout_t& fl) {
        Layout_t* layout = static_cast<Layout_t*>(&getLayout());
        if ( layout->getLayout().initialized() ) {
            std::cout << "Initialized" << std::endl;
        } else {
//             layout->getLayout().setFieldLayout(&fl);
            
            layout->getLayout().changeDomain(fl);
            std::cout << "Uninitialized" << std::endl;
        }
    }
    
    FieldLayout_t &getFieldLayout() {
        Layout_t* layout = static_cast<Layout_t*>(&getLayout());
        return dynamic_cast<FieldLayout_t &>(layout->getLayout().getFieldLayout());
    }
    
// private:
//     Mesh_t mesh_m;
//     FieldLayout_t fieldlayout_m;
    
};

#endif
