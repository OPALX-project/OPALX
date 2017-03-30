#include "AmrPartBunch.h"

#include "Utilities/OpalException.h"

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
    //TODO Implement
    throw OpalException("AmrPartBunch::computeSelfFields() ", "Not yet Implemented.");
}


void AmrPartBunch::computeSelfFields(int b) {
    //TODO Implement
    throw OpalException("AmrPartBunch::computeSelfFields(int) ", "Not yet Implemented.");
}


void AmrPartBunch::computeSelfFields_cycl(double gamma) {
    //TODO Implement
    throw OpalException("AmrPartBunch::computeSelfFields_cycl(double) ", "Not yet Implemented.");
}


void AmrPartBunch::computeSelfFields_cycl(int b) {
    //TODO Implement
    throw OpalException("AmrPartBunch::computeSelfFields_cycl(int) ", "Not yet Implemented.");
}


void AmrPartBunch::update() {
    BoxLibParticle<BoxLibLayout<double, 3> >::update();
}


void AmrPartBunch::update(const ParticleAttrib<char>& canSwap) {
    BoxLibParticle<BoxLibLayout<double, 3> >::update(canSwap);
}
