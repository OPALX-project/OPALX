#ifndef PARTBUNCHAMR_H
#define PARTBUNCHAMR_H

#include "AmrParticleBase.h"

template<class PLayout>
class PartBunchAmr : public AmrParticleBase<PLayout>
{

public:
  
  ParticleAttrib<double> qm;
  typename PLayout::ParticlePos_t E;
  typename PLayout::ParticlePos_t P;

  PartBunchAmr() {
    this->addAttribute(qm);
    this->addAttribute(E);
    this->addAttribute(P);
  }

  ~PartBunchAmr() {}

};


#endif
