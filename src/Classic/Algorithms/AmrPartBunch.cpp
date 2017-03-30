#include "AmrPartBunch.h"


AmrPartBunch::AmrPartBunch(const PartData *ref) : PartBunch(ref)
{
    initializeAmr();
}

AmrPartBunch::AmrPartBunch(const std::vector<OpalParticle> &rhs, const PartData *ref)
    : PartBunch(rhs, ref)
{
    initializeAmr();
}

AmrPartBunch::AmrPartBunch(const AmrPartBunch &rhs) : PartBunch(rhs)
{
    initializeAmr();
}

void AmrPartBunch::computeSelfFields() {
    
    
}


void AmrPartBunch::computeSelfFields(int b) {
    
    
    
}


void AmrPartBunch::computeSelfFields_cycl(double gamma) {
    
    
}


void AmrPartBunch::computeSelfFields_cycl(int b) {
    
    
    
}


void AmrPartBunch::update() {
    BoxLibParticle<BoxLibLayout<double, 3> >::update();
}

void AmrPartBunch::update(const ParticleAttrib<char>& canSwap) {
    BoxLibParticle<BoxLibLayout<double, 3> >::update(canSwap);
}
