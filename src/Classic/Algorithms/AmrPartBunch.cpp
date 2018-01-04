#include "AmrPartBunch.h"

#include "Utilities/OpalException.h"

AmrPartBunch::AmrPartBunch(const PartData *ref)
    : PartBunchBase<double, 3>(new AmrPartBunch::pbase_t(new AmrLayout_t()), ref),
      amrobj_mp(nullptr),
      amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase)),
      fieldlayout_m(nullptr)
{
    amrpbase_mp->initializeAmr();
}


AmrPartBunch::AmrPartBunch(const std::vector<OpalParticle> &rhs,
                           const PartData *ref)
    : PartBunchBase<double, 3>(new AmrPartBunch::pbase_t(new AmrLayout_t()), rhs, ref),
      amrobj_mp(nullptr),
      amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase)),
      fieldlayout_m(nullptr)
{
    amrpbase_mp->initializeAmr();
}


AmrPartBunch::AmrPartBunch(const AmrPartBunch &rhs)
    : PartBunchBase<double, 3>(rhs),
      amrobj_mp(nullptr),
      amrpbase_mp(dynamic_cast<AmrPartBunch::pbase_t*>(pbase)),
      fieldlayout_m(nullptr)
{
    amrpbase_mp->initializeAmr();
}

AmrPartBunch::~AmrPartBunch() {
    delete amrpbase_mp;
}


AmrPartBunch::pbase_t *AmrPartBunch::getAmrParticleBase() {
    return amrpbase_mp;
}


const AmrPartBunch::pbase_t *AmrPartBunch::getAmrParticleBase() const {
    return amrpbase_mp;
}


void AmrPartBunch::initialize(FieldLayout_t *fLayout) {
    Layout_t* layout = static_cast<Layout_t*>(&getLayout());
}


void AmrPartBunch::do_binaryRepart() {
    
    if ( amrobj_mp ) {
        
        const int& maxLevel = amrobj_mp->maxLevel();
        
        if ( !amrobj_mp->isRefined() ) {
            /* In the first call to this function we
             * intialize all fine levels
             */
            amrobj_mp->initFineLevels();
            
        } else if ( maxLevel > 0 && this->numBunch_m > 1 )  {
            /* we do an explicit domain mapping of the particles and then
             * forbid it during the regrid process, this way it's only
             * executed ones --> saves computation
             * 
             * allow only multi-level in case of multi-bunch
             * simulation
             */
            bool isForbidTransform = amrpbase_mp->isForbidTransform();
            
            if ( !isForbidTransform ) {
                amrpbase_mp->domainMapping();
                amrpbase_mp->setForbidTransform(true);
            }
            
            /* Update first in order to make
             * sure that the particles belong to the right
             * level and grid
             */
            this->update();
            
            int lev_top = std::min(amrobj_mp->finestLevel(), maxLevel - 1);
            
            *gmsg << "* Start regriding:" << endl
                  << "*     Old finest level: "
                  << amrobj_mp->finestLevel() << endl;
            
            /* ATTENTION: The bunch has to be updated during
             * the regrid process!
             * We regrid from base level 0 up to the finest level.
             */
            amrobj_mp->regrid(0, lev_top, t_m * 1.0e9 /*time [ns] */);
                
            *gmsg << "*     New finest level: "
                  << amrobj_mp->finestLevel() << endl
                  << "* Finished regriding" << endl;
            
            if ( !isForbidTransform ) {
                amrpbase_mp->setForbidTransform(false);
                // map particles back
                amrpbase_mp->domainMapping(true);
            }
        }
    }
//     amrobj_mp->redistributeGrids(-1 /*KnapSack*/);
//     update();
}


Vector_t AmrPartBunch::get_hr() const {
    const double& scalefactor = amrpbase_mp->getScalingFactor();
    return hr_m * scalefactor;
}


void AmrPartBunch::set_meshEnlargement(double dh) {
    // set dh_m = dh
    PartBunchBase<double, 3>::set_meshEnlargement(dh);
    
    // update base geometry with new mesh enlargement
    AmrLayout_t* layout_p = &amrpbase_mp->getAmrLayout();
    layout_p->setBoundingBox(dh);
    
    // if amrobj_mp != nullptr --> we need to regrid
    this->do_binaryRepart();
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
    
    if ( amrobj_mp ) {
        /* we do an explicit domain mapping of the particles and then
         * forbid it during the regrid process, this way it's only
         * executed ones --> saves computation
         */
        bool isForbidTransform = amrpbase_mp->isForbidTransform();
            
        if ( !isForbidTransform ) {
            amrpbase_mp->domainMapping();
            amrpbase_mp->setForbidTransform(true);
        }
        
        this->update();
        
        if ( !isForbidTransform ) {
            amrpbase_mp->setForbidTransform(false);
            // map particles back
            amrpbase_mp->domainMapping(true);
        }
        
    } else {
        // At this point an amrobj_mp needs already be set
        throw GeneralClassicException("AmrPartBunch::boundp() ",
                                      "AmrObject pointer is not set.");
    }
    
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
    
    reduce(*partPerLevel.get(),
           *globalPartPerLevel_m.get(),
           nLevel, std::plus<size_t>());
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
