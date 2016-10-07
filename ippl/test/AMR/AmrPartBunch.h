#ifndef AMRPARTBUNCH_H
#define AMRPARTBUNCH_H

class AmrPartBunch : public PartBunchBase,
                     /* 0: charge-to-mass ratio
                      * 1:3: particle velocity
                      * 4:6: electric field at particle position
                      * 7:9: magnetic field at particle position
                      */
                     public ParticleContainer<10>{
    
    
    
};

#endif