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


void AmrPartBunch::runTests() {
    //TODO Implement
    throw OpalException("AmrPartBunch::runTests() ", "Not yet Implemented.");
}


void AmrPartBunch::calcLineDensity(unsigned int nBins,
                                   std::vector<double> &lineDensity,
                                   std::pair<double, double> &meshInfo)
{
    //TODO Implement
    throw OpalException("&AmrPartBunch::calcLineDensity() ", "Not yet Implemented.");
}


// const Mesh_t &AmrPartBunch::getMesh() const {
//     //TODO Implement
//     throw OpalException("&AmrPartBunch::getMesh() ", "Not yet Implemented.");
// }


// Mesh_t &AmrPartBunch::getMesh() {
//     //TODO Implement
//     throw OpalException("AmrPartBunch::getMesh() ", "Not yet Implemented.");
// }


// FieldLayout_t &AmrPartBunch::getFieldLayout() {
//     //TODO Implement
//     throw OpalException("&AmrPartBunch::getFieldLayout() ", "Not yet Implemented.");
// }


void AmrPartBunch::boundp() {
    //TODO Implement
    throw OpalException("AmrPartBunch::boundp() ", "Not yet Implemented.");
}


void AmrPartBunch::boundp_destroy() {
    //TODO Implement
    throw OpalException("AmrPartBunch::boundp_destroy() ", "Not yet Implemented.");
}


size_t AmrPartBunch::boundp_destroyT() {
    //TODO Implement
    throw OpalException("AmrPartBunch::boundp_destroyT() ", "Not yet Implemented.");
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
