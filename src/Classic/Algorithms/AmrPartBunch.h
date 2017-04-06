#ifndef AMR_PART_BUNCH_H
#define AMR_PART_BUNCH_H

#include "Algorithms/PartBunchBase.h"
#include "Amr/BoxLibParticle.h"
#include "Amr/BoxLibLayout.h"

class AmrPartBunch : public PartBunchBase<double, 3>
{
public:
    
    AmrPartBunch(const PartData *ref);

    /// Conversion.
    AmrPartBunch(const std::vector<OpalParticle> &, const PartData *ref);

    AmrPartBunch(const AmrPartBunch &);
    
    VectorPair_t PartBunch::getEExtrema();
    
    double getRho(int x, int y, int z);
    
//     const Mesh_t &getMesh() const;

//     Mesh_t &getMesh();
    
//     FieldLayout_t &getFieldLayout();
    
    
    void computeSelfFields();
    
    void computeSelfFields(int b);
    
    void computeSelfFields_cycl(double gamma);
    
    void computeSelfFields_cycl(int b);
    
private:
    void updateFieldContainers_m();
    
    void updateDomainLength(Vector_t& grid);
    
private:
    
    /// charge density on the grid for all levels
    AmrFieldContainer_t rho_m;
    
    /// scalar potential on the grid for all levels
    AmrFieldContainer_t phi_m;
    
    /// vector field on the grid for all levels
    AmrFieldContainer_t eg_m;
    
};

#endif
