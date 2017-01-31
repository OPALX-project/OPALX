#ifndef PARTBUNCHAMR_H
#define PARTBUNCHAMR_H

#include "AmrParticleBase.h"

template<class PLayout>
class PartBunchAmr : public AmrParticleBase<PLayout>
{

public:
  
  ParticleAttrib<double> qm;

  PartBunchAmr() {
    this->addAttribute(qm);
  }

  ~PartBunchAmr() {}

};


#endif
