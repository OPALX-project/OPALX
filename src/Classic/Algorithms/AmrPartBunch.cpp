#include "AmrPartBunch.h"

#include "Utilities/OpalException.h"

AmrPartBunch::AmrPartBunch(const PartData *ref)
    : PartBunchBase<double, 3>(new AmrPartBunch::pbase_t(new AmrLayout_t()), ref),
      mesh_m(), fieldlayout_m(mesh_m), //FIXME
      amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase))
{
    amrpbase_mp->initializeAmr();
}


AmrPartBunch::AmrPartBunch(const std::vector<OpalParticle> &rhs,
                           const PartData *ref)
    : PartBunchBase<double, 3>(new AmrPartBunch::pbase_t(new AmrLayout_t()), rhs, ref),
      mesh_m(), fieldlayout_m(mesh_m), //FIXME
      amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase))
{
    amrpbase_mp->initializeAmr();
}


AmrPartBunch::AmrPartBunch(const AmrPartBunch &rhs)
    : PartBunchBase<double, 3>(rhs),
      mesh_m(), fieldlayout_m(mesh_m), //FIXME
      amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase))
{
    amrpbase_mp->initializeAmr();
}

AmrPartBunch::~AmrPartBunch() {

}


AmrPartBunch::pbase_t *AmrPartBunch::getAmrParticleBase() {
    return amrpbase_mp;
}


const AmrPartBunch::pbase_t *AmrPartBunch::getAmrParticleBase() const {
    return amrpbase_mp;
}

// AmrPartBunch::pbase_t* AmrPartBunch::clone() {
//     return new pbase_t(new BoxLibLayout<double, 3>());
// }


void AmrPartBunch::initialize(FieldLayout_t *fLayout) {
    Layout_t* layout = static_cast<Layout_t*>(&getLayout());
    layout->getLayout().changeDomain(*fLayout);
}


AmrPartBunch::VectorPair_t AmrPartBunch::getEExtrema() {
    return amrobj_mp->getEExtrema();
}

double AmrPartBunch::getRho(int x, int y, int z) {
    /* This function is called in
     * H5PartWrapperForPC::writeStepData(PartBunchBase<double, 3>* bunch)
     * and
     * H5PartWrapperForPT::writeStepData(PartBunchBase<double, 3>* bunch)
     * in case of Options::rhoDump = true.
     * 
     * Currently, we do not support writing multilevel grid data that's why
     * we average the values down to the coarsest level.
     */
    return amrobj_mp->getRho(x, y, z);
}

// const Mesh_t &AmrPartBunch::getMesh() const {
//     //TODO Implement
//     throw OpalException("&AmrPartBunch::getMesh() ", "Not yet Implemented.");
// }


// Mesh_t &AmrPartBunch::getMesh() {
//     //TODO Implement
//     throw OpalException("AmrPartBunch::getMesh() ", "Not yet Implemented.");
// }


// void AmrPartBunch::setMesh(Mesh_t* mesh) {
//     Layout_t* layout = static_cast<Layout_t*>(&getLayout());
//     layout->getLayout().setMesh(mesh);
// }
// 
// 
// void AmrPartBunch::setFieldLayout(FieldLayout_t* fLayout) {
//     Layout_t* layout = static_cast<Layout_t*>(&getLayout());
//     layout->getLayout().setFieldLayout(fLayout);
//     layout->getLayout().changeDomain(*fLayout);
// }


FieldLayout_t &AmrPartBunch::getFieldLayout() {
    //TODO Implement
    throw OpalException("&AmrPartBunch::getFieldLayout() ", "Not yet Implemented.");
    return fieldlayout_m;
}


void AmrPartBunch::computeSelfFields() {
    if ( !fs_m->hasValidSolver() )
        throw OpalException("AmrPartBunch::computeSelfFields() ", "No field solver.");
    
    amrobj_mp->computeSelfFields();
}


void AmrPartBunch::computeSelfFields(int b) {
    amrobj_mp->computeSelfFields(b);
}


void AmrPartBunch::computeSelfFields_cycl(double gamma) {
    amrobj_mp->computeSelfFields_cycl(gamma);
}


void AmrPartBunch::computeSelfFields_cycl(int b) {
    amrobj_mp->computeSelfFields_cycl(b);
}


void AmrPartBunch::updateFieldContainers_m() {
    
    
    
}

void AmrPartBunch::updateDomainLength(Vektor<int, 3>& grid) {
    
}


void AmrPartBunch::updateFields(const Vector_t& hr, const Vector_t& origin) {
    //TODO regrid; called in boundp()
    
}