#ifndef PARTBUNCH_H
#define PARTBUNCH_H

template<PLayout>
class PartBunch : public ParticleBase<PLayout>
                  public PartBunchBase
{

public:
    
    PartBunch(PLayout* pl, e_dim_tag decomp[Dim], bool gCells) :
        ParticleBase<PLayout>(pl),
        fieldNotInitialized_m(true),
        doRepart_m(true),
        withGuardCells_m(gCells)
    {
        // register the particle attributes
        this->addAttribute(qm);
        this->addAttribute(P);
        this->addAttribute(E);
        this->addAttribute(B);
        setupBCs();
        
        for (unsigned int i=0; i<Dim; i++)
            decomp_m[i]=decomp[i];
    
    
private:
    ParticleAttrib<double>     qm; // charge-to-mass ratio
    typename PL::ParticlePos_t P;  // particle velocity
    typename PL::ParticlePos_t E;  // electric field at particle position
    typename PL::ParticlePos_t B;  // magnetic field at particle position
    
    
    
    
    
};

#endif