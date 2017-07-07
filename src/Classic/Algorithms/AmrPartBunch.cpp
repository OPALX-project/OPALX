#include "AmrPartBunch.h"

#include "Utilities/OpalException.h"

AmrPartBunch::AmrPartBunch(const PartData *ref)
    : PartBunchBase<double, 3>(new AmrPartBunch::pbase_t(new AmrLayout_t()), ref),
      fieldlayout_m(nullptr),
      amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase))
{
    amrpbase_mp->initializeAmr();
}


AmrPartBunch::AmrPartBunch(const std::vector<OpalParticle> &rhs,
                           const PartData *ref)
    : PartBunchBase<double, 3>(new AmrPartBunch::pbase_t(new AmrLayout_t()), rhs, ref),
      fieldlayout_m(nullptr),
      amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase))
{
    amrpbase_mp->initializeAmr();
}


AmrPartBunch::AmrPartBunch(const AmrPartBunch &rhs)
    : PartBunchBase<double, 3>(rhs),
      fieldlayout_m(nullptr),
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
//     layout->getLayout().changeDomain(*fLayout);
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
    return *fieldlayout_m;
}


void AmrPartBunch::boundp() {
    IpplTimings::startTimer(boundpTimer_m);
    //if(!R.isDirty() && stateOfLastBoundP_ == unit_state_) return;
    if ( !(R.isDirty() || ID.isDirty() ) && stateOfLastBoundP_ == unit_state_) return; //-DW

    stateOfLastBoundP_ = unit_state_;
    
    /* we do an explicit domain mapping of the particles and then
     * forbid it during the regrid process, this way it's only
     * executed ones --> saves computation
     */
    AmrLayout_t* layout_p = &amrpbase_mp->getAmrLayout();
    
    layout_p->domainMapping(*amrpbase_mp);
    
    layout_p->setForbidTransform(true);
    
    // update all level
    update();
    
    if ( amrobj_mp ) {
        int maxLevel = amrobj_mp->maxLevel();
        
        for (int i = 0; i <= amrobj_mp->finestLevel() && i < maxLevel; ++i) {
            amrobj_mp->regrid(i, maxLevel, t_m * 1.0e9 /*time [ns] */);
            /* update to multilevel --> update GDB
             * Only update current and new finest level
             */
            amrpbase_mp->update(i, amrobj_mp->finestLevel());
        }
    } else {
        // At this point an amrobj_mp needs already be set
        throw GeneralClassicException("AmrPartBunch::boundp() ",
                                      "AmrObject pointer is not set.");
    }
    
    layout_p->setForbidTransform(false);
    
    // map particles back
    layout_p->domainMapping(*amrpbase_mp, true);
    
    R.resetDirtyFlag();
    
    IpplTimings::stopTimer(boundpTimer_m);
}


void AmrPartBunch::computeSelfFields() {
    IpplTimings::startTimer(selfFieldTimer_m);
    
    if ( !fs_m->hasValidSolver() )
        throw OpalException("AmrPartBunch::computeSelfFields() ",
                            "No field solver.");
    
    amrobj_mp->computeSelfFields();
    
    IpplTimings::stopTimer(selfFieldTimer_m);
}


void AmrPartBunch::computeSelfFields(int bin) {
    IpplTimings::startTimer(selfFieldTimer_m);
    amrobj_mp->computeSelfFields(bin);
    IpplTimings::stopTimer(selfFieldTimer_m);
}


void AmrPartBunch::computeSelfFields_cycl(double gamma) {
    IpplTimings::startTimer(selfFieldTimer_m);
    amrobj_mp->computeSelfFields_cycl(gamma);
    IpplTimings::stopTimer(selfFieldTimer_m);
}


void AmrPartBunch::computeSelfFields_cycl(int bin) {
    IpplTimings::startTimer(selfFieldTimer_m);
    amrobj_mp->computeSelfFields_cycl(bin);
    IpplTimings::stopTimer(selfFieldTimer_m);
}


void AmrPartBunch::gatherLevelStatistics() {
    int nLevel = (&amrpbase_mp->getAmrLayout())->maxLevel() + 1;
    
    std::unique_ptr<size_t[]> partPerLevel( new size_t[nLevel] );
    globalPartPerLevel_m.reset( new size_t[nLevel] );
    
    for (int i = 0; i < nLevel; ++i)
        partPerLevel[i] = globalPartPerLevel_m[i] = 0.0;
    
    // do not modify LocalNumPerLevel in here!!!
    auto& LocalNumPerLevel = amrpbase_mp->getLocalNumPerLevel();
        
    for (size_t i = 0; i < LocalNumPerLevel.size(); ++i)
        partPerLevel[i] = LocalNumPerLevel[i];
        
    reduce(partPerLevel.get(),
           partPerLevel.get() + nLevel,
           globalPartPerLevel_m.get(),
           OpAddAssign());
}


const size_t& AmrPartBunch::getLevelStatistics(int l) const {
    return globalPartPerLevel_m[l];
}


void AmrPartBunch::updateFieldContainers_m() {
    
}

void AmrPartBunch::updateDomainLength(Vektor<int, 3>& grid) {
    grid = amrobj_mp->getBaseLevelGridPoints();
}


void AmrPartBunch::updateFields(const Vector_t& hr, const Vector_t& origin) {
    //TODO regrid; called in boundp()
//     amrobj_mp->updateMesh();
}