#ifndef AMR_PARTICLE_BASE_HPP
#define AMR_PARTICLE_BASE_HPP

#include <numeric>
#include <algorithm>

template<class PLayout>
AmrParticleBase<PLayout>::AmrParticleBase() : forbidTransform_m(false),
                                              scale_m(1.0),
                                              LocalNumPerLevel_m()
{
    UpdateParticlesTimer_m = IpplTimings::getTimer("AMR update particles");
    SortParticlesTimer_m = IpplTimings::getTimer("AMR sort particles");
}


template<class PLayout>
AmrParticleBase<PLayout>::AmrParticleBase(PLayout* layout)
    : IpplParticleBase<PLayout>(layout),
      forbidTransform_m(false),
      scale_m(1.0),
      LocalNumPerLevel_m()
{
    UpdateParticlesTimer_m = IpplTimings::getTimer("AMR update particles");
    SortParticlesTimer_m = IpplTimings::getTimer("AMR sort particles");
}


template<class PLayout>
const typename AmrParticleBase<PLayout>::ParticleLevelCounter_t&
    AmrParticleBase<PLayout>::getLocalNumPerLevel() const
{
    return LocalNumPerLevel_m;
}


template<class PLayout>
typename AmrParticleBase<PLayout>::ParticleLevelCounter_t&
    AmrParticleBase<PLayout>::getLocalNumPerLevel()
{
    return LocalNumPerLevel_m;
}


template<class PLayout>
void AmrParticleBase<PLayout>::setLocalNumPerLevel(
    const ParticleLevelCounter_t& LocalNumPerLevel)
{
    LocalNumPerLevel_m = LocalNumPerLevel;
}


template<class PLayout>
void AmrParticleBase<PLayout>::destroy(size_t M, size_t I, bool doNow) {
    /* if the particles are deleted directly
     * we need to update the particle level count
     */
    if (M > 0) {
        if ( doNow ) {
            for (size_t ip = I; ip < M + I; ++ip)
                --LocalNumPerLevel_m[ Level[ip] ];
        }
        IpplParticleBase<PLayout>::destroy(M, I, doNow);
    }
}


template<class PLayout>
void AmrParticleBase<PLayout>::performDestroy(bool updateLocalNum) {
    // nothing to do if destroy list is empty
    if ( this->DestroyList.empty() )
        return;
    
    if ( updateLocalNum ) {
        typedef std::vector< std::pair<size_t,size_t> > dlist_t;
        dlist_t::const_iterator curr = this->DestroyList.begin();
        const dlist_t::const_iterator last = this->DestroyList.end();
        
        while ( curr != last ) {
            for (size_t ip = curr->first;
                 ip < curr->first + curr->second;
                 ++ip)
            {
                --LocalNumPerLevel_m[ Level[ip] ];
            }
            ++curr;
        }
    }
    IpplParticleBase<PLayout>::performDestroy(updateLocalNum);
}


template<class PLayout>
void AmrParticleBase<PLayout>::create(size_t M) {
    // particles are created at the coarsest level
    LocalNumPerLevel_m[0] += M;
    
    IpplParticleBase<PLayout>::create(M);
}


template<class PLayout>
void AmrParticleBase<PLayout>::createWithID(unsigned id) {
    ++LocalNumPerLevel_m[0];
    
    IpplParticleBase<PLayout>::createWithID(id);
}


template<class PLayout>
void AmrParticleBase<PLayout>::update() {
    // update all level
    this->update(0, -1);
}


template<class PLayout>
void AmrParticleBase<PLayout>::update(int lev_min, int lev_max) {
    
    IpplTimings::startTimer(UpdateParticlesTimer_m);

    // make sure we've been initialized
    PLayout *Layout = &this->getLayout();

    PAssert(Layout != 0);
    
    // ask the layout manager to update our atoms, etc.
    Layout->update(*this, lev_min, lev_max);
    
    //sort the particles by grid and level
    sort();
    
    INCIPPLSTAT(incParticleUpdates);
    
    IpplTimings::stopTimer(UpdateParticlesTimer_m);
}

template<class PLayout>
void AmrParticleBase<PLayout>::update(const ParticleAttrib<char>& canSwap) {
    
    IpplTimings::startTimer(UpdateParticlesTimer_m);

    // make sure we've been initialized
    PLayout *Layout = &this->getLayout();
    PAssert(Layout != 0);
    
    // ask the layout manager to update our atoms, etc.
    Layout->update(*this, &canSwap);
    
    //sort the particles by grid and level
    sort();
    
    INCIPPLSTAT(incParticleUpdates);
    
    IpplTimings::stopTimer(UpdateParticlesTimer_m);
}


template<class PLayout>
void AmrParticleBase<PLayout>::sort() {
    
    IpplTimings::startTimer(SortParticlesTimer_m);
    size_t LocalNum = this->getLocalNum();
    SortList_t slist1(LocalNum); //slist1 holds the index of where each element should go
    SortList_t slist2(LocalNum); //slist2 holds the index of which element should go in this position

    //sort the lists by grid and level
    //slist1 hold the index of where each element should go in the list
    std::iota(slist1.begin(), slist1.end(), 0);
    std::sort(slist1.begin(), slist1.end(), [this](const SortListIndex_t &i, 
                                                   const SortListIndex_t &j)
    {
        return (this->Level[i] < this->Level[j] ||
               (this->Level[i] == this->Level[j] && this->Grid[i] < this->Grid[j]));
    });

    //slist2 holds the index of which element should go in this position
    for (unsigned int i = 0; i < LocalNum; ++i)
        slist2[slist1[i]] = i;

    //sort the array according to slist2
    this->sort(slist2);

    IpplTimings::stopTimer(SortParticlesTimer_m);
}


template<class PLayout>
void AmrParticleBase<PLayout>::sort(SortList_t &sortlist) {
    attrib_container_t::iterator abeg = this->begin();
        attrib_container_t::iterator aend = this->end();
        for ( ; abeg != aend; ++abeg )
            (*abeg)->sort(sortlist);
}


template<class PLayout>
void AmrParticleBase<PLayout>::setForbidTransform(bool forbidTransform) {
    forbidTransform_m = forbidTransform;
}


template<class PLayout>
bool AmrParticleBase<PLayout>::isForbidTransform() const {
    return forbidTransform_m;
}


template<class PLayout>
const double& AmrParticleBase<PLayout>::domainMapping(bool inverse) {
    Vector_t rmin, rmax;
    bounds(this->R, rmin, rmax);
    
    double scale = scale_m;
    
    if ( !inverse ) {
        Vector_t tmp = Vector_t(std::max( std::abs(rmin[0]), std::abs(rmax[0]) ),
                                std::max( std::abs(rmin[1]), std::abs(rmax[1]) ),
                                std::max( std::abs(rmin[2]), std::abs(rmax[2]) )
                               );
        
        scale = std::max( tmp[0], tmp[1] );
        scale = std::max( scale, tmp[2] );
    }
    
    Vector_t vscale = Vector_t(scale, scale, scale);
    
    for (unsigned int i = 0; i < this->getLocalNum(); ++i)
        this->R[i] /= vscale;
    
    
    scale_m = 1.0 / scale;
    
    return scale_m;
}


template<class PLayout>
const double& AmrParticleBase<PLayout>::getScalingFactor() const {
    return scale_m;
}

#endif
