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
    
    /* virtual member functions of PartBunch that
     * we reimplement
     */
    
    void runTests();
    
    void calcLineDensity(unsigned int nBins,
                         std::vector<double> &lineDensity,
                         std::pair<double, double> &meshInfo);
    
//     const Mesh_t &getMesh() const;

//     Mesh_t &getMesh();
    
//     FieldLayout_t &getFieldLayout();
    
    void boundp();
    
    void boundp_destroy();
    
    size_t boundp_destroyT();
    
    void computeSelfFields();
    
    void computeSelfFields(int b);
    
    void computeSelfFields_cycl(double gamma);
    
    void computeSelfFields_cycl(int b);
    
    /* virtual member functions of IpplParticleBase
     * that we reimplement
     */
    
    void update();
    
    void update(const ParticleAttrib<char>& canSwap);
};

#endif
