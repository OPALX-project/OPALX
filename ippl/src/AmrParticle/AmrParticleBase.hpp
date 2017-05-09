#ifndef AMR_PARTICLE_BASE_HPP
#define AMR_PARTICLE_BASE_HPP

#include <numeric>
#include <algorithm>

template<class PLayout>
AmrParticleBase<PLayout>::AmrParticleBase() {
    UpdateParticlesTimer_m = IpplTimings::getTimer("AMR update particles");
    SortParticlesTimer_m = IpplTimings::getTimer("AMR sort particles");
}


template<class PLayout>
void AmrParticleBase<PLayout>::update() {
    
    IpplTimings::startTimer(UpdateParticlesTimer_m);

    // make sure we've been initialized
    PLayout *Layout = &this->getLayout();

    PAssert(Layout != 0);
    
    // ask the layout manager to update our atoms, etc.
    Layout->update(*this);
    
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

#endif
