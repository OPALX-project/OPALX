#include "AmrPartBunch.h"

#include "Utilities/OpalException.h"

AmrPartBunch::AmrPartBunch(const PartData *ref)
    : PartBunchBase<double, 3>(ref),
      rho_m(PArrayManage),
      phi_m(PArrayManage),
      eg_m(PArrayManage)
{
    initializeAmr();
}


AmrPartBunch::AmrPartBunch(const std::vector<OpalParticle> &rhs, const PartData *ref)
    : PartBunchBase<double, 3>(rhs, ref),
      rho_m(PArrayManage),
      phi_m(PArrayManage),
      eg_m(PArrayManage)
{
    initializeAmr();
}


AmrPartBunch::AmrPartBunch(const AmrPartBunch &rhs)
    : PartBunchBase<double, 3>(rhs),
      rho_m(PArrayManage),
      phi_m(PArrayManage),
      eg_m(PArrayManage)
{
    initializeAmr();
}

AmrPartBunch::VectorPair_t AmrPartBunch::getEExtrema() {
    return VectorPair_t();
}


double AmrPartBunch::getRho(int x, int y, int z) {
    throw OpalException("AmrPartBunch::getRho(x, y, z) ", "Not yet Implemented.");
    return 0.0;
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


void AmrPartBunch::computeSelfFields() {
    
//     if ( !fs_m->hasValidSolver() )
//         throw OpalException("AmrPartBunch::computeSelfFields() ", "No field solver.");
//     
//     IpplTimings::startTimer(selfFieldTimer_m);
//     
//     //scatter charges onto grid
//     this->Q *= this->dt;
//     this->scatter(this->Q, this->rho_m, /*BoxLibParticle<BoxLibLayout<double, 3u> >::*/this->R, 0, -1);
//     this->Q /= this->dt;
//     
//     int baseLevel = 0;
//     int finestLevel = (&this->getAmrLayout())->finestLevel();
//     
//     int nLevel = finestLevel + 1;
//     rho_m.resize(nLevel);
//     phi_m.resize(nLevel);
//     eg_m.resize(nLevel);
//     
//     double invDt = 1.0 / getdT() * getCouplingConstant();
//     for (int i = 0; i <= finestLevel; ++i) {
//         this->rho_m[i].mult(invDt, 0, 1);
//     }
//     
//     // charge density is in rho_m
//     IpplTimings::startTimer(compPotenTimer_m);
//     
// //     fs_m->solver_m->solve(rho_m, eg_m, baseLevel, finestLevel);
// 
//     IpplTimings::stopTimer(compPotenTimer_m);
//     
//     
//     
//     
//     IpplTimings::stopTimer(selfFieldTimer_m);
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


void AmrPartBunch::updateFieldContainers_m() {
    
    
    
}

void AmrPartBunch::updateDomainLength(Vector_t& grid) {
    
}
