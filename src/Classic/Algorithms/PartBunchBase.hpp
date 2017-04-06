#ifndef PART_BUNCH_BASE_HPP
#define PART_BUNCH_BASE_HPP

#include "PartBunchBase.h"

#include "Utilities/OpalException.h"

template <class T, unsigned Dim>
PartBunchBase<T, Dim>::PartBunchBase(AbstractParticle<T, Dim>* pb)
    : pbase(pb),
      R(*(pb->R_p)),
      ID(*(pb->ID_p))
{
    pb->addAttribute(P);
    pb->addAttribute(Q);
    pb->addAttribute(M);
    pb->addAttribute(Phi);
    pb->addAttribute(Ef);
    pb->addAttribute(Eftmp);
    
    pb->addAttribute(Bf);
    pb->addAttribute(Bin);
    pb->addAttribute(dt);
    pb->addAttribute(PType);
    pb->addAttribute(TriID);
}


/*
 * Bunch common member functions
 */

template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::getIfBeamEmitting() {
    if (dist_m != NULL) {
        size_t isBeamEmitted = dist_m->getIfDistEmitting();
        reduce(isBeamEmitted, isBeamEmitted, OpAddAssign());
        if (isBeamEmitted > 0)
            return true;
        else
            return false;
    } else
        return false;
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getLastEmittedEnergyBin() {
    /*
     * Get maximum last emitted energy bin.
     */
    int lastEmittedBin = dist_m->getLastEmittedEnergyBin();
    reduce(lastEmittedBin, lastEmittedBin, OpMaxAssign());
    return lastEmittedBin;
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::getNumberOfEmissionSteps() {
    return dist_m->getNumberOfEmissionSteps();
}


template <class T, unsigned Dim>
int PartBunchBase<T, Dim>::getNumberOfEnergyBins() {
    return dist_m->getNumberOfEnergyBins();
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::Rebin() {

    size_t isBeamEmitting = dist_m->getIfDistEmitting();
    reduce(isBeamEmitting, isBeamEmitting, OpAddAssign());
    if (isBeamEmitting > 0) {
        *gmsg << "*****************************************************" << endl
              << "Warning: attempted to rebin, but not all distribution" << endl
              << "particles have been emitted. Rebin failed." << endl
              << "*****************************************************" << endl;
    } else {
        if (dist_m->Rebin())
            this->Bin = 0;
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setEnergyBins(int numberOfEnergyBins) {
    bingamma_m = std::unique_ptr<double[]>(new double[numberOfEnergyBins]);
    binemitted_m = std::unique_ptr<size_t[]>(new size_t[numberOfEnergyBins]);
    for(int i = 0; i < numberOfEnergyBins; i++)
        binemitted_m[i] = 0;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::weHaveEnergyBins() {

    if (dist_m != NULL)
        return dist_m->getNumberOfEnergyBins() > 0;
    else
        return false;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::switchToUnitlessPositions(bool use_dt_per_particle) {

    if(unit_state_ == unitless)
        throw SwitcherError("PartBunch::switchToUnitlessPositions",
                            "Cannot make a unitless PartBunch unitless");

    bool hasToReset = false;
    if(!R.isDirty()) hasToReset = true;

    for(size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if(use_dt_per_particle)
            dt = this->dt[i];

        R[i] /= Vector_t(Physics::c * dt);
    }

    unit_state_ = unitless;

    if(hasToReset) R.resetDirtyFlag();
}


//FIXME: unify methods, use convention that all particles have own dt
template <class T, unsigned Dim>
// inline
void PartBunch::switchOffUnitlessPositions(bool use_dt_per_particle) {

    if(unit_state_ == units)
        throw SwitcherError("PartBunch::switchOffUnitlessPositions",
                            "Cannot apply units twice to PartBunch");

    bool hasToReset = false;
    if(!R.isDirty()) hasToReset = true;

    for(size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if(use_dt_per_particle)
            dt = this->dt[i];

        R[i] *= Vector_t(Physics::c * dt);
    }

    unit_state_ = units;

    if(hasToReset) R.resetDirtyFlag();
}


/** \brief After each Schottky scan we delete all the particles.

 */
template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::cleanUpParticles() {

    size_t np = getTotalNum();
    bool scan = false;

    destroy(getLocalNum(), 0, true);

    dist_m->createOpalT(*this, np, scan);

    update();
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::resetIfScan()
/*
  In case of a scan we have
  to reset some variables
 */
{
    dt = 0.0;
    localTrackStep_m = 0;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::do_binaryRepart() {
    get_bounds(rmin_m, rmax_m);
    BinaryRepartition(*this);
    update();
    get_bounds(rmin_m, rmax_m);
    boundp();
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setDistribution(Distribution *d,
                                            std::vector<Distribution *> addedDistributions,
                                            size_t &np,
                                            bool scan)
{
    Inform m("setDistribution " );
    dist_m = d;

    dist_m->createOpalT(*this, addedDistributions, np, scan);

//    if (Options::cZero)
//        dist_m->create(*this, addedDistributions, np / 2, scan);
//    else
//        dist_m->create(*this, addedDistributions, np, scan);
}


template <class T, unsigned Dim>
// inline
bool PartBunchBase<T, Dim>::isGridFixed() {
    return fixed_grid;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setTEmission(double t) {
    tEmission_m = t;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getTEmission() {
    return tEmission_m;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::doEmission() {
    if (dist_m != NULL)
        return dist_m->getIfDistEmitting();
    else
        return false;
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::weHaveBins() const {

    if(pbin_m != NULL)
        return pbin_m->weHaveBins();
    else
        return false;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setPBins(PartBins *pbin) {
    pbin_m = pbin;
    *gmsg << *pbin_m << endl;
    setEnergyBins(pbin_m->getNBins());
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setPBins(PartBinsCyc *pbin) {
    pbin_m = pbin;
    setEnergyBins(pbin_m->getNBins());
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::emitParticles(double eZ) {
    return dist_m->emitParticles(*this, eZ);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::updateNumTotal() {
    size_t numParticles = getLocalNum();
    reduce(numParticles, numParticles, OpAddAssign());
    setTotalNum(numParticles);
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::rebin() {
    this->Bin = 0;
    pbin_m->resetBins();
    // delete pbin_m; we did not allocate it!
    pbin_m = NULL;
}


template <class T, unsigned Dim>
// inline
int PartBunchBase<T, Dim>::getNumBins() {
    if(pbin_m != NULL)
        return pbin_m->getNBins();
    else
        return 0;
}


template <class T, unsigned Dim>
// inline
int PartBunchBase<T, Dim>::getLastemittedBin() {
    if(pbin_m != NULL)
        return pbin_m->getLastemittedBin();
    else
        return 0;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcGammas() {

    const int emittedBins = dist_m->getNumberOfEnergyBins();

    size_t s = 0;

    for(int i = 0; i < emittedBins; i++)
        bingamma_m[i] = 0.0;

    for(unsigned int n = 0; n < getLocalNum(); n++)
        bingamma_m[this->Bin[n]] += sqrt(1.0 + dot(this->P[n], this->P[n]));

    std::unique_ptr<size_t[]> particlesInBin(new size_t[emittedBins]);
    reduce(bingamma_m.get(), bingamma_m.get() + emittedBins, bingamma_m.get(), OpAddAssign());
    reduce(binemitted_m.get(), binemitted_m.get() + emittedBins, particlesInBin.get(), OpAddAssign());
    for(int i = 0; i < emittedBins; i++) {
        size_t &pInBin = particlesInBin[i];
        if(pInBin != 0) {
            bingamma_m[i] /= pInBin;
            INFOMSG(level2 << "Bin " << std::setw(3) << i
                           << " gamma = " << std::setw(8) << std::scientific
                           << std::setprecision(5) << bingamma_m[i]
                           << "; NpInBin= " << std::setw(8)
                           << std::setfill(' ') << pInBin << endl);
        } else {
            bingamma_m[i] = 1.0;
            INFOMSG(level2 << "Bin " << std::setw(3) << i << " has no particles " << endl);
        }
        s += pInBin;
    }
    particlesInBin.reset();


    if(s != getTotalNum() && !OpalData::getInstance()->hasGlobalGeometry())
        ERRORMSG("sum(Bins)= " << s << " != sum(R)= " << getTotalNum() << endl;);

    if(emittedBins >= 2) {
        for(int i = 1; i < emittedBins; i++) {
            if(binemitted_m[i - 1] != 0 && binemitted_m[i] != 0)
                INFOMSG(level2 << "d(gamma)= " << 100.0 * std::abs(bingamma_m[i - 1] - bingamma_m[i]) / bingamma_m[i] << " [%] "
                        << "between bin " << i - 1 << " and " << i << endl);
        }
    }
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcGammas_cycl() {

    const int emittedBins = pbin_m->getLastemittedBin();

    for(int i = 0; i < emittedBins; i++)
        bingamma_m[i] = 0.0;
    for(unsigned int n = 0; n < getLocalNum(); n++)
        bingamma_m[this->Bin[n]] += sqrt(1.0 + dot(this->P[n], this->P[n]));
    for(int i = 0; i < emittedBins; i++) {
        reduce(bingamma_m[i], bingamma_m[i], OpAddAssign());
        if(pbin_m->getTotalNumPerBin(i) > 0)
            bingamma_m[i] /= pbin_m->getTotalNumPerBin(i);
        else
            bingamma_m[i] = 0.0;
        INFOMSG("Bin " << i << " : particle number = " << pbin_m->getTotalNumPerBin(i) << " gamma = " << bingamma_m[i] << endl);
    }

}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getBinGamma(int bin) {
    return bingamma_m[bin];
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setBinCharge(int bin, double q) {
    this->Q = where(eq(this->Bin, bin), q, 0.0);
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setBinCharge(int bin) {
    this->Q = where(eq(this->Bin, bin), this->Q, 0.0);
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::calcNumPartsOutside(Vector_t x) {

    partPerNode_m[Ippl::myNode()] = 0;
    const Vector_t meanR = get_rmean();

    for(unsigned long k = 0; k < getLocalNum(); ++ k)
        if (abs(R[k](0) - meanR(0)) > x(0) ||
            abs(R[k](1) - meanR(1)) > x(1) ||
            abs(R[k](2) - meanR(2)) > x(2)) {

            ++ partPerNode_m[Ippl::myNode()];
        }

    reduce(partPerNode_m.get(), partPerNode_m.get() + Ippl::getNodes(), globalPartPerNode_m.get(), OpAddAssign());

    return *globalPartPerNode_m.get();
}


/**
 * \method calcLineDensity()
 * \brief calculates the 1d line density (not normalized) and append it to a file.
 * \see ParallelTTracker
 * \warning none yet
 *
 * DETAILED TODO
 *
 */
template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcLineDensity(unsigned int nBins,
                                            std::vector<double> &lineDensity,
                                            std::pair<double, double> &meshInfo)
{
    Vector_t rmin, rmax;
    get_bounds(rmin, rmax);

    if (nBins < 2) {
        NDIndex<3> grid;
        this->updateDomainLength(grid);
        nBins = grid[2];
    }

    double length = rmax(2) - rmin(2);
    double zmin = rmin(2) - dh_m * length, zmax = rmax(2) + dh_m * length;
    double hz = (zmax - zmin) / (nBins - 2);
    double perMeter = 1.0 / hz;//(zmax - zmin);
    zmin -= hz;

    lineDensity.resize(nBins, 0.0);
    std::fill(lineDensity.begin(), lineDensity.end(), 0.0);

    const unsigned int lN = getLocalNum();
    for (unsigned int i = 0; i < lN; ++ i) {
        const double z = R[i](2) - 0.5 * hz;
        unsigned int idx = (z - zmin) / hz;
        double tau = z - (zmin + idx * hz);

        lineDensity[idx] += Q[i] * (1.0 - tau) * perMeter;
        lineDensity[idx + 1] += Q[i] * tau * perMeter;
    }

    reduce(&(lineDensity[0]), &(lineDensity[0]) + nBins, &(lineDensity[0]), OpAddAssign());

    meshInfo.first = zmin;
    meshInfo.second = hz;
}



template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::boundp() {
    /*
      Assume rmin_m < 0.0
     */

    IpplTimings::startTimer(boundpTimer_m);
    //if(!R.isDirty() && stateOfLastBoundP_ == unit_state_) return;
    if ( !(R.isDirty() || ID.isDirty() ) && stateOfLastBoundP_ == unit_state_) return; //-DW

    stateOfLastBoundP_ = unit_state_;

    if(!isGridFixed()) {
        const int dimIdx = 3;

        /**
           In case of dcBeam_m && hr_m < 0
           this is the first call to boundp and we
           have to set hr completely i.e. x,y and z.
           */

        const bool fullUpdate = (dcBeam_m && (hr_m[2] < 0.0)) || !dcBeam_m;
        
        this->updateDomainLength(nr_m);
        
        get_bounds(rmin_m, rmax_m);
        Vector_t len = rmax_m - rmin_m;
        
        double volume = 1.0;
        if (fullUpdate) {
            // double volume = 1.0;
            for(int i = 0; i < dimIdx; i++) {
                double length = std::abs(rmax_m[i] - rmin_m[i]);
                if (length < 1e-10) {
                    rmax_m[i] += 1e-10;
                    rmin_m[i] -= 1e-10;
                } else {
                    rmax_m[i] += dh_m * length;
                    rmin_m[i] -= dh_m * length;
                }
                hr_m[i]    = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
                volume *= std::abs(rmax_m[i] - rmin_m[i]);
            }

            if (getIfBeamEmitting() && dist_m != NULL) {
                // keep particles per cell ratio high, don't spread a hand full particles across the whole grid
                double percent = std::max((1.0 + (3 - nr_m[2]) * dh_m) / (nr_m[2] - 1), dist_m->getPercentageEmitted());
                double length   = std::abs(rmax_m[2] - rmin_m[2]);
                if (percent < 1.0 && percent > 0.0) {
                    length /= (1.0 + 2 * dh_m);
                    rmax_m[2] -= dh_m * length;
                    rmin_m[2] = rmax_m[2] * (1.0 - 1.0 / percent);

                    length = std::abs(rmax_m[2] - rmin_m[2]);
                    rmax_m[2] += dh_m * length;
                    rmin_m[2] -= dh_m * length;
                    hr_m[2] = length * (1.0 + 2 * dh_m) / (nr_m[2] - 1);
                }
            }

            if (volume < 1e-21 && getTotalNum() > 1 && std::abs(sum(Q)) > 0.0) {
                WARNMSG(level1 << "!!! Extremly high particle density detected !!!" << endl);
            }
            //INFOMSG("It is a full boundp hz= " << hr_m << " rmax= " << rmax_m << " rmin= " << rmin_m << endl);
        }
        
        if(hr_m[0] * hr_m[1] * hr_m[2] > 0) {
            Vector_t origin = rmin_m - Vector_t(hr_m[0] / 2.0, hr_m[1] / 2.0, hr_m[2] / 2.0);
            this->updateFields(hr, origin);
        } else {
            *gmsg << __DBGMSG__ << std::scientific << volume << "\t" << dh_m << endl;
            
            throw GeneralClassicException("boundp() ", "h<0, can not build a mesh");
        }
    }
    update();
    R.resetDirtyFlag();

    IpplTimings::stopTimer(boundpTimer_m);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::boundp_destroy() {

    Inform gmsgAll("boundp_destroy ", INFORM_ALL_NODES);

    Vector_t len;
    const int dimIdx = 3;
    IpplTimings::startTimer(boundpTimer_m);

    std::unique_ptr<size_t[]> countLost;
    if(weHaveBins()) {
        const int tempN = pbin_m->getLastemittedBin();
        countLost = std::unique_ptr<size_t[]>(new size_t[tempN]);
        for(int ii = 0; ii < tempN; ii++) countLost[ii] = 0;
    }
    
    this->updateDomainLength(nr_m);

    get_bounds(rmin_m, rmax_m);
    len = rmax_m - rmin_m;

    calcBeamParameters();

    int checkfactor = Options::remotePartDel;
    if (checkfactor != 0) {
        //INFOMSG("checkfactor= " << checkfactor << endl);
        // check the bunch if its full size is larger than checkfactor times of its rms size
        if(checkfactor < 0) {
            checkfactor *= -1;
            if (len[0] > checkfactor * rrms_m[0] ||
                len[1] > checkfactor * rrms_m[1] ||
                len[2] > checkfactor * rrms_m[2])
            {
                for(unsigned int ii = 0; ii < this->getLocalNum(); ii++) {
                    /* delete the particle if the ditance to the beam center
                     * is larger than 8 times of beam's rms size
                     */
                    if (abs(R[ii](0) - rmean_m(0)) > checkfactor * rrms_m[0] ||
                        abs(R[ii](1) - rmean_m(1)) > checkfactor * rrms_m[1] ||
                        abs(R[ii](2) - rmean_m(2)) > checkfactor * rrms_m[2])
                    {
                        // put particle onto deletion list
                        destroy(1, ii);
                        //update bin parameter
                        if (weHaveBins())
                            countLost[Bin[ii]] += 1 ;
                        /* INFOMSG("REMOTE PARTICLE DELETION: ID = " << ID[ii] << ", R = " << R[ii]
                         * << ", beam rms = " << rrms_m << endl;);
                         */
                    }
                }
            }
        }
        else {
            if (len[0] > checkfactor * rrms_m[0] ||
                len[2] > checkfactor * rrms_m[2])
            {
                for(unsigned int ii = 0; ii < this->getLocalNum(); ii++) {
                    /* delete the particle if the ditance to the beam center
                     * is larger than 8 times of beam's rms size
                     */
                    if (abs(R[ii](0) - rmean_m(0)) > checkfactor * rrms_m[0] ||
                        abs(R[ii](2) - rmean_m(2)) > checkfactor * rrms_m[2])
                    {
                        // put particle onto deletion list
                        destroy(1, ii);
                        //update bin parameter
                        if (weHaveBins())
                            countLost[Bin[ii]] += 1 ;
                        /* INFOMSG("REMOTE PARTICLE DELETION: ID = " << ID[ii] << ", R = " << R[ii]
                         * << ", beam rms = " << rrms_m << endl;);
                         */
                    }
                }
            }
        }
    }

    for(int i = 0; i < dimIdx; i++) {
        double length = std::abs(rmax_m[i] - rmin_m[i]);
        rmax_m[i] += dh_m * length;
        rmin_m[i] -= dh_m * length;
        hr_m[i]    = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
    }

    // rescale mesh
    this->updateFields(hr, rmin_m);
    
    if(weHaveBins()) {
        pbin_m->updatePartInBin_cyc(countLost.get());
    }
    
    update();
    
    IpplTimings::stopTimer(boundpTimer_m);
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::boundp_destroyT() {

    const unsigned int minNumParticlesPerCore = getMinimumNumberOfParticlesPerCore();

    this->updateDomainLength(nr_m);

    std::unique_ptr<size_t[]> tmpbinemitted;

    boundp();

    size_t ne = 0;
    const size_t localNum = getLocalNum();

    if(weHaveEnergyBins()) {
        tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[getNumberOfEnergyBins()]);
        for(int i = 0; i < getNumberOfEnergyBins(); i++) {
            tmpbinemitted[i] = 0;
        }
        for(unsigned int i = 0; i < localNum; i++) {
            if (Bin[i] < 0) {
                ne++;
                destroy(1, i);
            } else
                tmpbinemitted[Bin[i]]++;
        }
    } else {
        for(unsigned int i = 0; i < localNum; i++) {
            if((Bin[i] < 0) && ((localNum - ne) > minNumParticlesPerCore)) {   // need in minimum x particles per node
                ne++;
                destroy(1, i);
            }
        }
        lowParticleCount_m = ((localNum - ne) <= minNumParticlesPerCore);
        reduce(lowParticleCount_m, lowParticleCount_m, OpOr());
    }

    if (lowParticleCount_m) {
        Inform m ("boundp_destroyT a) ", INFORM_ALL_NODES);
        m << level3 << "Warning low number of particles on some cores localNum= "
          << localNum << " ne= " << ne << " NLocal= " << getLocalNum() << endl;
    } else {
        boundp();
    }
    calcBeamParameters();
    gatherLoadBalanceStatistics();

    if(weHaveEnergyBins()) {
        const int lastBin = dist_m->getLastEmittedEnergyBin();
        for(int i = 0; i < lastBin; i++) {
            binemitted_m[i] = tmpbinemitted[i];
        }
    }
    reduce(ne, ne, OpAddAssign());
    return ne;
}


template <class T, unsigned Dim>
size_t PartBunchBase<T, Dim>::destroyT() {

    const unsigned int minNumParticlesPerCore = getMinimumNumberOfParticlesPerCore();
    std::unique_ptr<size_t[]> tmpbinemitted;

    const size_t localNum = getLocalNum();
    const size_t totalNum = getTotalNum();

    if(weHaveEnergyBins()) {
        tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[getNumberOfEnergyBins()]);
        for(int i = 0; i < getNumberOfEnergyBins(); i++) {
            tmpbinemitted[i] = 0;
        }
        for(unsigned int i = 0; i < localNum; i++) {
            if (Bin[i] < 0) {
                destroy(1, i);
            } else
                tmpbinemitted[Bin[i]]++;
        }
    } else {
        Inform dmsg("destroy: ", INFORM_ALL_NODES);
        size_t ne = 0;
        for(size_t i = 0; i < localNum; i++) {
            if((Bin[i] < 0) && ((localNum - ne) > minNumParticlesPerCore)) {   // need in minimum x particles per node
                ne++;
                destroy(1, i);
            }
        }
        lowParticleCount_m = ((localNum - ne) <= minNumParticlesPerCore);
        reduce(lowParticleCount_m, lowParticleCount_m, OpOr());
        // unsigned int i = 0, j = localNum;
        // while (j > 0 && Bin[j - 1] < 0) -- j;

        // while (i + 1 < j) {
        //     if (Bin[i] < 0) {
        //         this->swap(i,j - 1);
        //         -- j;

        //         while (i + 1 < j && Bin[j - 1] < 0) -- j;
        //     }
        //     ++ i;
        // }

        // j = std::max(j, minNumParticlesPerCore);
        // for(unsigned int i = localNum; i > j; -- i) {
        //     destroy(1, i-1, true);
        // }
        // lowParticleCount_m = (j == minNumParticlesPerCore);
        // reduce(lowParticleCount_m, lowParticleCount_m, OpOr());

        if (ne > 0) {
            performDestroy(true);
        }
    }

    calcBeamParameters();
    gatherLoadBalanceStatistics();

    if (weHaveEnergyBins()) {
        const int lastBin = dist_m->getLastEmittedEnergyBin();
        for(int i = 0; i < lastBin; i++) {
            binemitted_m[i] = tmpbinemitted[i];
        }
    }
    size_t newTotalNum = getLocalNum();
    reduce(newTotalNum, newTotalNum, OpAddAssign());

    setTotalNum(newTotalNum);

    return totalNum - newTotalNum;
}

template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPx(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPy(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getPz(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getPx0(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getPy0(int i) {
    return 0;
}


//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getX(int i) {
    return this->R[i](0);
}


//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getY(int i) {
    return this->R[i](1);
}


//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getZ(int i) {
    return this->R[i](2);
}


//ff
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getX0(int i) {
    return 0.0;
}


//ff
template <class T, unsigned Dim>
double PartBunch<T, Dim>::getY0(int i) {
    return 0.0;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setZ(int i, double zcoo)
{
    // nothing done here
};


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::get_bounds(Vector_t &rmin, Vector_t &rmax) {
    bounds(this->R, rmin, rmax);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::getLocalBounds(Vector_t &rmin, Vector_t &rmax) {
    const size_t localNum = getLocalNum();
    if (localNum == 0) return;

    rmin = R[0];
    rmax = R[0];
    for (size_t i = 1; i < localNum; ++ i) {
        for (unsigned short d = 0; d < 3u; ++ d) {
            if (rmin(d) > R[i](d)) rmin(d) = R[i](d);
            else if (rmax(d) < R[i](2)) rmax(d) = R[i](d);
        }
    }
}


template <class T, unsigned Dim>
// inline
std::pair<Vector_t, double> PartBunchBase<T, Dim>::getBoundingSphere() {
    Vector_t rmin, rmax;
    get_bounds(rmin, rmax);

    std::pair<Vector_t, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


template <class T, unsigned Dim>
// inline
std::pair<Vector_t, double> PartBunchBase<T, Dim>::getLocalBoundingSphere() {
    Vector_t rmin, rmax;
    getLocalBounds(rmin, rmax);

    std::pair<Vector_t, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::push_back(Particle p) {
    Inform msg("PartBunch ");

    create(1);
    size_t i = getTotalNum();

    R[i](0) = p[0];
    R[i](1) = p[1];
    R[i](2) = p[2];

    P[i](0) = p[3];
    P[i](1) = p[4];
    P[i](2) = p[5];

    update();
    msg << "Created one particle i= " << i << endl;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::set_part(FVector<double, 6> z, int ii) {
    R[ii](0) = z[0];
    P[ii](0) = z[1];
    R[ii](1) = z[2];
    P[ii](1) = z[3];
    R[ii](2) = z[4];
    P[ii](2) = z[5];
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::set_part(Particle p, int ii) {
    R[ii](0) = p[0];
    P[ii](0) = p[1];
    R[ii](1) = p[2];
    P[ii](1) = p[3];
    R[ii](2) = p[4];
    P[ii](2) = p[5];
}


template <class T, unsigned Dim>
// inline
Particle PartBunchBase<T, Dim>::get_part(int ii) {
    Particle part;
    part[0] = R[ii](0);
    part[1] = P[ii](0);
    part[2] = R[ii](1);
    part[3] = P[ii](1);
    part[4] = R[ii](2);
    part[5] = P[ii](2);
    return part;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::maximumAmplitudes(const FMatrix<double, 6, 6> &D,
                                              double &axmax, double &aymax) {
    axmax = aymax = 0.0;
    Particle part;

    for(unsigned int ii = 0; ii < getLocalNum(); ii++) {

        part = get_part(ii);

        double xnor =
            D(0, 0) * part.x()  + D(0, 1) * part.px() + D(0, 2) * part.y() +
            D(0, 3) * part.py() + D(0, 4) * part.t()  + D(0, 5) * part.pt();
        double pxnor =
            D(1, 0) * part.x()  + D(1, 1) * part.px() + D(1, 2) * part.y() +
            D(1, 3) * part.py() + D(1, 4) * part.t()  + D(1, 5) * part.pt();
        double ynor =
            D(2, 0) * part.x()  + D(2, 1) * part.px() + D(2, 2) * part.y() +
            D(2, 3) * part.py() + D(2, 4) * part.t()  + D(2, 5) * part.pt();
        double pynor =
            D(3, 0) * part.x()  + D(3, 1) * part.px() + D(3, 2) * part.y() +
            D(3, 3) * part.py() + D(3, 4) * part.t()  + D(3, 5) * part.pt();

        axmax = std::max(axmax, (xnor * xnor + pxnor * pxnor));
        aymax = std::max(aymax, (ynor * ynor + pynor * pynor));
    }
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setdT(double dt) {
    dt_m = dt;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getdT() const {
    return dt_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setT(double t) {
    t_m = t;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getT() const {
    return t_m;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::get_sPos() {
    if(sum(PType != ParticleType::REGULAR)) {
        const size_t n = getLocalNum();
        size_t numPrimBeamParts = 0;
        double z = 0.0;
        if(n != 0) {
            for(size_t i = 0; i < n; i++) {
                if(PType[i] == ParticleType::REGULAR) {
                    z += R[i](2);
                    numPrimBeamParts++;
                }
            }
        }
        reduce(z, z, OpAddAssign());
	reduce(numPrimBeamParts, numPrimBeamParts, OpAddAssign());
	if(numPrimBeamParts != 0)
            z = z / numPrimBeamParts;
        return z;
    } else {
        return spos_m;
    }
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::set_sPos(double s) {
    spos_m = s;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::get_gamma() const {
    return eKin_m / (getM()*1e-6) + 1.0;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::get_meanKineticEnergy() const {
    return eKin_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_origin() const {
    return rmin_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_maxExtent() const {
    return rmax_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_centroid() const {
    return rmean_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_rrms() const {
    return rrms_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_rprms() const {
    return rprms_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_rmean() const {
    return rmean_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_prms() const {
    return prms_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_pmean() const {
    return pmean_m;
}


template <class T, unsigned Dim>
Vector_t PartBunchBase<T, Dim>::get_pmean_Distribution() const {
    if (dist_m && dist_m->getType() != DistrTypeT::FROMFILE)
        return dist_m->get_pmean();

    double gamma = 0.1 / getM() + 1; // set default 0.1 eV
    return Vector_t(0, 0, sqrt(std::pow(gamma, 2) - 1));
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_emit() const {
    return eps_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_norm_emit() const {
    return eps_norm_m;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::get_hr() const {
    return hr_m;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::get_Dx() const {
    return Dx_m;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::get_Dy() const {
    return Dy_m;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::get_DDx() const {
    return DDx_m;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::get_DDy() const {
    return DDy_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::set_meshEnlargement(double dh) {
    dh_m = dh;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::gatherLoadBalanceStatistics() {
    minLocNum_m =  std::numeric_limits<size_t>::max();

    for(int i = 0; i < Ippl::getNodes(); i++)
        partPerNode_m[i] = globalPartPerNode_m[i] = 0;

    partPerNode_m[Ippl::myNode()] = getLocalNum();

    reduce(partPerNode_m.get(), partPerNode_m.get() + Ippl::getNodes(), globalPartPerNode_m.get(), OpAddAssign());

    for(int i = 0; i < Ippl::getNodes(); i++)
        if (globalPartPerNode_m[i] <  minLocNum_m)
            minLocNum_m = globalPartPerNode_m[i];
}


template <class T, unsigned Dim>
// inline
size_t PartBunchBase<T, Dim>::getLoadBalance(int p) const {
    return globalPartPerNode_m[p];
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::get_PBounds(Vector_t &min, Vector_t &max) const {
    bounds(this->P, min, max);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcBeamParameters() {
    using Physics::c;

    Vector_t eps2, fac, rsqsum, psqsum, rpsum;

    const double m0 = getM() * 1.E-6;

    IpplTimings::startTimer(statParamTimer_m);

    const size_t locNp = getLocalNum();
    const double N =  static_cast<double>(getTotalNum());
    const double zero = 0.0;

    get_bounds(rmin_m, rmax_m);

    if(N == 0) {
        for(unsigned int i = 0 ; i < Dim; i++) {
            rmean_m(i) = 0.0;
            pmean_m(i) = 0.0;
            rrms_m(i) = 0.0;
            prms_m(i) = 0.0;
            eps_norm_m(i)  = 0.0;
        }
        rprms_m = 0.0;
        eKin_m = 0.0;
        eps_m = 0.0;
        IpplTimings::stopTimer(statParamTimer_m);
        return;
    }

    calcMoments();

    for(unsigned int i = 0 ; i < Dim; i++) {
        rmean_m(i) = centroid_m[2 * i] / N;
        pmean_m(i) = centroid_m[(2 * i) + 1] / N;
        rsqsum(i) = moments_m(2 * i, 2 * i) - N * rmean_m(i) * rmean_m(i);
        psqsum(i) = moments_m((2 * i) + 1, (2 * i) + 1) - N * pmean_m(i) * pmean_m(i);
        if(psqsum(i) < 0)
            psqsum(i) = 0;
        rpsum(i) = moments_m((2 * i), (2 * i) + 1) - N * rmean_m(i) * pmean_m(i);
    }
    eps2 = (rsqsum * psqsum - rpsum * rpsum) / (N * N);
    rpsum /= N;

    for(unsigned int i = 0 ; i < Dim; i++) {
        rrms_m(i) = sqrt(rsqsum(i) / N);
        prms_m(i) = sqrt(psqsum(i) / N);
        eps_norm_m(i)  =  std::sqrt(std::max(eps2(i), zero));
        double tmp = rrms_m(i) * prms_m(i);
        fac(i) = (tmp == 0) ? zero : 1.0 / tmp;
    }

    rprms_m = rpsum * fac;

    Dx_m = moments_m(0, 5) / N;
    DDx_m = moments_m(1, 5) / N;

    Dy_m = moments_m(2, 5) / N;
    DDy_m = moments_m(3, 5) / N;

    // Find unnormalized emittance.
    double gamma = 0.0;
    for(size_t i = 0; i < locNp; i++)
        gamma += sqrt(1.0 + dot(P[i], P[i]));

    reduce(gamma, gamma, OpAddAssign());
    gamma /= N;

    calcEMean();

    // calculate energy spread
    dE_m = prms_m(2) * sqrt(eKin_m * (eKin_m + 2.*m0) / (1. + eKin_m * (eKin_m + 2.*m0) / (m0 * m0)));

    eps_m = eps_norm_m / Vector_t(gamma * sqrt(1.0 - 1.0 / (gamma * gamma)));
    IpplTimings::stopTimer(statParamTimer_m);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcBeamParametersInitial() {
    using Physics::c;

    const double N =  static_cast<double>(getTotalNum());

    if(N == 0) {
        rmean_m = Vector_t(0.0);
        pmean_m = Vector_t(0.0);
        rrms_m  = Vector_t(0.0);
        prms_m  = Vector_t(0.0);
        eps_m   = Vector_t(0.0);
    } else {
        if(Ippl::myNode() == 0) {
            // fixme:  the following code is crap!
            // Only use one node as this function will get called only once before
            // particles have been emitted (at least in principle).
            Vector_t eps2, fac, rsqsum, psqsum, rpsum;

            const double zero = 0.0;
            const double  N =  static_cast<double>(pbin_m->getNp());
            calcMomentsInitial();

            for(unsigned int i = 0 ; i < Dim; i++) {
                rmean_m(i) = centroid_m[2 * i] / N;
                pmean_m(i) = centroid_m[(2 * i) + 1] / N;
                rsqsum(i) = moments_m(2 * i, 2 * i) - N * rmean_m(i) * rmean_m(i);
                psqsum(i) = moments_m((2 * i) + 1, (2 * i) + 1) - N * pmean_m(i) * pmean_m(i);
                if(psqsum(i) < 0)
                    psqsum(i) = 0;
                rpsum(i) =  moments_m((2 * i), (2 * i) + 1) - N * rmean_m(i) * pmean_m(i);
            }
            eps2 = (rsqsum * psqsum - rpsum * rpsum) / (N * N);
            rpsum /= N;

            for(unsigned int i = 0 ; i < Dim; i++) {

                rrms_m(i) = sqrt(rsqsum(i) / N);
                prms_m(i) = sqrt(psqsum(i) / N);
                eps_m(i)  = sqrt(std::max(eps2(i), zero));
                double tmp = rrms_m(i) * prms_m(i);
                fac(i) = (tmp == 0) ? zero : 1.0 / tmp;
            }
            rprms_m = rpsum * fac;
        }
    }
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getCouplingConstant() const {
    return couplingConstant_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setCouplingConstant(double c) {
    couplingConstant_m = c;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setCharge(double q) {
    qi_m = q;
    if(getTotalNum() != 0)
        Q = qi_m;
    else
        WARNMSG("Could not set total charge in PartBunch::setCharge based on getTotalNum" << endl);
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setChargeZeroPart(double q) {
    qi_m = q;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setMass(double mass) {
    M = mass;
}


/// \brief Need Ek for the Schottky effect calculation (eV)
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getEkin() const {
    if(dist_m)
        return dist_m->getEkin();
    else
        return 0.0;
}

/// \brief Need the work function for the Schottky effect calculation (eV)
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getWorkFunctionRf() const {
    if(dist_m)
        return dist_m->getWorkFunctionRf();
    else
        return 0.0;
}
/// \brief Need the laser energy for the Schottky effect calculation (eV)
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getLaserEnergy() const {
    if(dist_m)
        return dist_m->getLaserEnergy();
    else
        return 0.0;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getCharge() const {
    return sum(Q);
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getChargePerParticle() const {
    return qi_m;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setSolver(FieldSolver *fs) {
    fs_m = fs;
    fs_m->initSolver(*this);
    /**
       CAN not re-inizialize ParticleLayout
       this is an IPPL issue
     */
    if(!OpalData::getInstance()->hasBunchAllocated())
        initialize(&(fs_m->getParticleLayout()));
}


template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::hasFieldSolver() {
    if(fs_m)
        return fs_m->hasValidSolver();
    else
        return false;
}


/// \brief Return the fieldsolver type if we have a fieldsolver
template <class T, unsigned Dim>
std::string PartBunchBase<T, Dim>::getFieldSolverType() const {
    if(fs_m)
        return fs_m->getFieldSolverType();
    else
        return "";
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setLPath(double s) {
    lPath_m = s;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getLPath() const {
    return lPath_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setStepsPerTurn(int n) {
    stepsPerTurn_m = n;
}


template <class T, unsigned Dim>
// inline
int PartBunchBase<T, Dim>::getStepsPerTurn() const {
    return stepsPerTurn_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setGlobalTrackStep(long long n) {
    globalTrackStep_m = n;
}


template <class T, unsigned Dim>
// inline
long long PartBunchBase<T, Dim>::getGlobalTrackStep() const {
    return globalTrackStep_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setLocalTrackStep(long long n) {
    localTrackStep_m = n;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::incTrackSteps() {
    globalTrackStep_m++; localTrackStep_m++;
}


template <class T, unsigned Dim>
// inline
long long PartBunchBase<T, Dim>::getLocalTrackStep() const {
    return localTrackStep_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setNumBunch(int n) {
    numBunch_m = n;
}


template <class T, unsigned Dim>
// inline
int PartBunchBase<T, Dim>::getNumBunch() const {
    return numBunch_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setGlobalMeanR(Vector_t globalMeanR) {
    globalMeanR_m = globalMeanR;
}


template <class T, unsigned Dim>
// inline
Vector_t PartBunchBase<T, Dim>::getGlobalMeanR() {
    return globalMeanR_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion) {

    globalToLocalQuaternion_m = globalToLocalQuaternion;
}


template <class T, unsigned Dim>
// inline
Quaternion_t PartBunchBase<T, Dim>::getGlobalToLocalQuaternion() {
    return globalToLocalQuaternion_m;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::setSteptoLastInj(int n) {
    SteptoLastInj_m = n;
}


template <class T, unsigned Dim>
// inline
int PartBunchBase<T, Dim>::getSteptoLastInj() {
    return SteptoLastInj_m;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::calcMeanPhi() {

    const int emittedBins = pbin_m->getLastemittedBin();
    double phi[emittedBins];
    double px[emittedBins];
    double py[emittedBins];
    double meanPhi = 0.0;

    for(int ii = 0; ii < emittedBins; ii++) {
        phi[ii] = 0.0;
        px[ii] = 0.0;
        py[ii] = 0.0;
    }

    for(unsigned int ii = 0; ii < getLocalNum(); ii++) {

        px[Bin[ii]] += P[ii](0);
        py[Bin[ii]] += P[ii](1);
    }

    reduce(px, px + emittedBins, px, OpAddAssign());
    reduce(py, py + emittedBins, py, OpAddAssign());
    for(int ii = 0; ii < emittedBins; ii++) {
        phi[ii] = calculateAngle(px[ii], py[ii]);
        meanPhi += phi[ii];
        INFOMSG("Bin " << ii  << "  mean phi = " << phi[ii] * 180.0 / pi - 90.0 << "[degree]" << endl);
    }

    meanPhi /= emittedBins;

    INFOMSG("mean phi of all particles " <<  meanPhi * 180.0 / pi - 90.0 << "[degree]" << endl);

    return meanPhi;
}


// this function reset the BinID for each particles according to its current beta*gamma
// it is for multi-turn extraction cyclotron with small energy gain
// the bin number can be different with the bunch number

template <class T, unsigned Dim>
bool PartBunchBase<T, Dim>::resetPartBinID2(const double eta) {


    INFOMSG("Before reset Bin: " << endl);
    calcGammas_cycl();
    int maxbin = pbin_m->getNBins();
    size_t partInBin[maxbin];
    for(int ii = 0; ii < maxbin; ii++) partInBin[ii] = 0;

    double pMin0 = 1.0e9;
    double pMin = 0.0;
    double maxbinIndex = 0;

    for(unsigned long int n = 0; n < getLocalNum(); n++) {
        double temp_betagamma = sqrt(pow(P[n](0), 2) + pow(P[n](1), 2));
        if(pMin0 > temp_betagamma)
            pMin0 = temp_betagamma;
    }
    reduce(pMin0, pMin, OpMinAssign());
    INFOMSG("minimal beta*gamma = " << pMin << endl);

    double asinh0 = asinh(pMin);
    for(unsigned long int n = 0; n < getLocalNum(); n++) {

        double temp_betagamma = sqrt(pow(P[n](0), 2) + pow(P[n](1), 2));

        int itsBinID = floor((asinh(temp_betagamma) - asinh0) / eta + 1.0E-6);
        Bin[n] = itsBinID;
        if(maxbinIndex < itsBinID) {
            maxbinIndex = itsBinID;
        }

        if(itsBinID >= maxbin) {
            ERRORMSG("The bin number limit is " << maxbin << ", please increase the energy interval and try again" << endl);
            return false;
        } else
            partInBin[itsBinID]++;

    }

    // partInBin only count particle on the local node.
    pbin_m->resetPartInBin_cyc(partInBin, maxbinIndex);

    // after reset Particle Bin ID, update mass gamma of each bin again
    INFOMSG("After reset Bin: " << endl);
    calcGammas_cycl();

    return true;

}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getQ() const {
    return reference->getQ();
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getM() const {
    return reference->getM();
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getP() const {
    return reference->getP();
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getE() const {
    return reference->getE();
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::resetQ(double q)  {
    const_cast<PartData *>(reference)->setQ(q);
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::resetM(double m)  {
    const_cast<PartData *>(reference)->setM(m);
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getdE() {
    return dE_m;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getInitialBeta() const {
    return reference->getBeta();
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getInitialGamma() const {
    return reference->getGamma();
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getGamma(int i) {
    return 0;
}


template <class T, unsigned Dim>
// inline
double PartBunchBase<T, Dim>::getBeta(int i) {
    return 0;
}


template <class T, unsigned Dim>
// inline
void PartBunchBase<T, Dim>::actT()
{
    // do nothing here
};


template <class T, unsigned Dim>
// inline
const PartData *PartBunchBase<T, Dim>::getReference() const {
    return reference;
}


template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::getEmissionDeltaT() {
    return dist_m->getEmissionDeltaT();
}


template <class T, unsigned Dim>
/*inline */Quaternion_t PartBunchBase<T, Dim>::getQKs3D() {
    return QKs3D_m;
}


template <class T, unsigned Dim>
/*inline*/ void PartBunchBase<T, Dim>::setQKs3D(Quaternion_t q) {
    QKs3D_m=q;
}


template <class T, unsigned Dim>
/*inline*/ Vector_t PartBunchBase<T, Dim>::getKs3DRefr() {
    return Ks3DRefr_m;
}


template <class T, unsigned Dim>
/*inline*/ void PartBunchBase<T, Dim>::setKs3DRefr(Vector_t r) {
    Ks3DRefr_m=r;
}


template <class T, unsigned Dim>
/*inline*/ Vector_t PartBunchBase<T, Dim>::getKs3DRefp() {
    return Ks3DRefp_m;
}


template <class T, unsigned Dim>
/*inline*/ void PartBunchBase<T, Dim>::setKs3DRefp(Vector_t p) {
    Ks3DRefp_m=p;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::iterateEmittedBin(int binNumber) {
    binemitted_m[binNumber]++;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::calcEMean() {

    const double totalNp = static_cast<double>(getTotalNum());
    const double locNp = static_cast<double>(getLocalNum());

    //Vector_t meanP_temp = Vector_t(0.0);

    eKin_m = 0.0;

    for(unsigned int k = 0; k < locNp; k++)
        //meanP_temp += P[k];
        eKin_m += sqrt(dot(P[k], P[k]) + 1.0);

    eKin_m -= locNp;
    eKin_m *= getM() * 1.0e-6;

    //reduce(meanP_temp, meanP_temp, OpAddAssign());
    reduce(eKin_m, eKin_m, OpAddAssign());

    //meanP_temp /= totalNp;
    eKin_m /= totalNp;

    //eKin_m = (sqrt(dot(meanP_temp, meanP_temp) + 1.0) - 1.0) * getM() * 1.0e-6;
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::correctEnergy(double avrgp_m) {

    const double totalNp = static_cast<double>(getTotalNum());
    const double locNp = static_cast<double>(getLocalNum());

    double avrgp = 0.0;
    for(unsigned int k = 0; k < locNp; k++)
        avrgp += sqrt(dot(P[k], P[k]));

    reduce(avrgp, avrgp, OpAddAssign());
    avrgp /= totalNp;

    for(unsigned int k = 0; k < locNp; k++)
        P[k](2) =  P[k](2) - avrgp + avrgp_m;
}


// angle range [0~2PI) degree
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::calculateAngle(double x, double y) {
    double thetaXY = atan2(y, x);

    // if(x < 0)                   thetaXY = pi + atan(y / x);
    // else if((x > 0) && (y >= 0))  thetaXY = atan(y / x);
    // else if((x > 0) && (y < 0))   thetaXY = 2.0 * pi + atan(y / x);
    // else if((x == 0) && (y > 0)) thetaXY = pi / 2.0;
    // else if((x == 0) && (y < 0)) thetaXY = 3.0 / 2.0 * pi;

    return thetaXY >= 0 ? thetaXY : thetaXY + Physics::two_pi;

}


// angle range [-PI~PI) degree
template <class T, unsigned Dim>
double PartBunchBase<T, Dim>::calculateAngle2(double x, double y) {

    // double thetaXY = atan2(y, x);
    // if(x > 0)              thetaXY = atan(y / x);
    // else if((x < 0)  && (y > 0)) thetaXY = pi + atan(y / x);
    // else if((x < 0)  && (y <= 0)) thetaXY = -pi + atan(y / x);
    // else if((x == 0) && (y > 0)) thetaXY = pi / 2.0;
    // else if((x == 0) && (y < 0)) thetaXY = -pi / 2.0;

    return atan2(y, x);

}



template <class T, unsigned Dim>
// inline
Inform &operator<<(Inform &os, PartBunchBase<T, Dim> &p) {
    return p.print(os);
}


/*
 * Virtual member functions 
 */

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::runTests() {
    throw OpalException("PartBunchBase<T, Dim>::runTests() ", "No test supported.");
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::resetInterpolationCache(bool clearCache) {
    
}

template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::swap(unsigned int i, unsigned int j) {
    if (i >= getLocalNum() || j >= getLocalNum() || i == j) return;

    std::swap(R[i], R[j]);
    std::swap(P[i], P[j]);
    std::swap(Q[i], Q[j]);
    std::swap(M[i], M[j]);
    std::swap(Phi[i], Phi[j]);
    std::swap(Ef[i], Ef[j]);
    std::swap(Eftmp[i], Eftmp[j]);
    std::swap(Bf[i], Bf[j]);
    std::swap(Bin[i], Bin[j]);
    std::swap(dt[i], dt[j]);
    std::swap(PType[i], PType[j]);
    std::swap(TriID[i], TriID[j]);
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setBCAllPeriodic() {
    throw OpalException("PartBunchBase<T, Dim>::setBCAllPeriodic() ", "Not supported BC.");
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setBCAllOpen() {
    throw OpalException("PartBunchBase<T, Dim>::setBCAllOpen() ", "Not supported BC.");
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::setBCForDCBeam() {
    throw OpalException("PartBunchBase<T, Dim>::setBCForDCBeam() ", "Not supported BC.");
}


template <class T, unsigned Dim>
void PartBunchBase<T, Dim>::updateFields()
{
    
}

// --------------------------------------------------------------------





int PartBunchBase::getNumberOfEnergyBins() {
    return dist_m->getNumberOfEnergyBins();
}


void PartBunchBase::Rebin() {

    size_t isBeamEmitting = dist_m->getIfDistEmitting();
    reduce(isBeamEmitting, isBeamEmitting, OpAddAssign());
    if (isBeamEmitting > 0) {
        *gmsg << "*****************************************************" << endl
              << "Warning: attempted to rebin, but not all distribution" << endl
              << "particles have been emitted. Rebin failed." << endl
              << "*****************************************************" << endl;
    } else {
        if (dist_m->Rebin())
            this->Bin = 0;
    }
}


void PartBunchBase::setEnergyBins(int numberOfEnergyBins) {
    bingamma_m = std::unique_ptr<double[]>(new double[numberOfEnergyBins]);
    binemitted_m = std::unique_ptr<size_t[]>(new size_t[numberOfEnergyBins]);
    for(int i = 0; i < numberOfEnergyBins; i++)
        binemitted_m[i] = 0;
}


bool PartBunchBase::weHaveEnergyBins() {

    if (dist_m != NULL)
        return dist_m->getNumberOfEnergyBins() > 0;
    else
        return false;
}


inline
void PartBunchBase::switchOffUnitlessPositions(bool use_dt_per_particle) {

    if(unit_state_ == units)
        throw SwitcherError("PartBunch::switchOffUnitlessPositions",
                            "Cannot apply units twice to PartBunch");

    bool hasToReset = false;
    if(!bunch_mp->R.isDirty()) hasToReset = true;

    for(size_t i = 0; i < bunch_mp->getLocalNum(); i++) {
        double dt = getdT();
        if(use_dt_per_particle)
            dt = bunch_mp->dt[i];

        bunch_mp->R[i] *= Vector_t(Physics::c * dt);
    }

    unit_state_ = units;

    if(hasToReset) bunch_mp->R.resetDirtyFlag();
}

inline
void PartBunchBase::switchToUnitlessPositions(bool use_dt_per_particle) {

    if(unit_state_ == unitless)
        throw SwitcherError("PartBunch::switchToUnitlessPositions",
                            "Cannot make a unitless PartBunch unitless");

    bool hasToReset = false;
    if(!bunch_mp->R.isDirty()) hasToReset = true;

    for(size_t i = 0; i < bunch_mp->getLocalNum(); i++) {
        double dt = getdT();
        if(use_dt_per_particle)
            dt = bunch_mp->dt[i];

        bunch_mp->R[i] /= Vector_t(Physics::c * dt);
    }

    unit_state_ = unitless;

    if(hasToReset) bunch_mp->R.resetDirtyFlag();
}


void PartBunchBase::resetIfScan()
/*
  In case of a scan we have
  to reset some variables
 */
{
    bunch_mp->dt = 0.0;
    localTrackStep_m = 0;
}


inline
void PartBunch::do_binaryRepart() {
    get_bounds(rmin_m, rmax_m);
    BinaryRepartition(*this);
    bunch_mp->update();
    get_bounds(rmin_m, rmax_m);
    bunch_mp->boundp();
}


void PartBunchBase::setDistribution(Distribution *d,
                                    std::vector<Distribution *> addedDistributions,
                                    size_t &np,
                                    bool scan)
{
    Inform m("setDistribution " );
    dist_m = d;

    dist_m->createOpalT(*this, addedDistributions, np, scan);

//    if (Options::cZero)
//        dist_m->create(*this, addedDistributions, np / 2, scan);
//    else
//        dist_m->create(*this, addedDistributions, np, scan);
}


inline
void PartBunchBase::setTEmission(double t) {
    tEmission_m = t;
}


inline
double PartBunch::getTEmission() {
    return tEmission_m;
}


bool PartBunchBase::doEmission() {
    if (dist_m != NULL)
        return dist_m->getIfDistEmitting();
    else
        return false;
}


bool PartBunchBase::weHaveBins() const {

    if(pbin_m != NULL)
        return pbin_m->weHaveBins();
    else
        return false;
}


void PartBunchBase::setPBins(PartBins *pbin) {
    pbin_m = pbin;
    *gmsg << *pbin_m << endl;
    setEnergyBins(pbin_m->getNBins());
}


void PartBunchBase::setPBins(PartBinsCyc *pbin) {

    pbin_m = pbin;
    setEnergyBins(pbin_m->getNBins());
}


size_t PartBunchBase::emitParticles(double eZ) {

    return dist_m->emitParticles(*this, eZ);

}


void PartBunchBase::updateNumTotal() {
    size_t numParticles = bunch_mp->getLocalNum();
    reduce(numParticles, numParticles, OpAddAssign());
    bunch_mp->setTotalNum(numParticles);
}


inline
void PartBunchBase::rebin() {
    this->Bin = 0;
    pbin_m->resetBins();
    // delete pbin_m; we did not allocate it!
    pbin_m = NULL;
}


inline
int PartBunchBase::getNumBins() {
    if(pbin_m != NULL)
        return pbin_m->getNBins();
    else
        return 0;
}


inline
int PartBunchBase::getLastemittedBin() {
    if(pbin_m != NULL)
        return pbin_m->getLastemittedBin();
    else
        return 0;
}


void PartBunchBase::calcGammas() {

    const int emittedBins = dist_m->getNumberOfEnergyBins();

    size_t s = 0;

    for(int i = 0; i < emittedBins; i++)
        bingamma_m[i] = 0.0;

    for(unsigned int n = 0; n < getLocalNum(); n++)
        bingamma_m[this->Bin[n]] += sqrt(1.0 + dot(this->P[n], this->P[n]));

    std::unique_ptr<size_t[]> particlesInBin(new size_t[emittedBins]);
    reduce(bingamma_m.get(), bingamma_m.get() + emittedBins, bingamma_m.get(), OpAddAssign());
    reduce(binemitted_m.get(), binemitted_m.get() + emittedBins, particlesInBin.get(), OpAddAssign());
    for(int i = 0; i < emittedBins; i++) {
        size_t &pInBin = particlesInBin[i];
        if(pInBin != 0) {
            bingamma_m[i] /= pInBin;
            INFOMSG(level2 << "Bin " << std::setw(3) << i << " gamma = " << std::setw(8) << std::scientific << std::setprecision(5) << bingamma_m[i] << "; NpInBin= " << std::setw(8) << std::setfill(' ') << pInBin << endl);
        } else {
            bingamma_m[i] = 1.0;
            INFOMSG(level2 << "Bin " << std::setw(3) << i << " has no particles " << endl);
        }
        s += pInBin;
    }
    particlesInBin.reset();


    if(s != getTotalNum() && !OpalData::getInstance()->hasGlobalGeometry())
        ERRORMSG("sum(Bins)= " << s << " != sum(R)= " << getTotalNum() << endl;);

    if(emittedBins >= 2) {
        for(int i = 1; i < emittedBins; i++) {
            if(binemitted_m[i - 1] != 0 && binemitted_m[i] != 0)
                INFOMSG(level2 << "d(gamma)= " << 100.0 * std::abs(bingamma_m[i - 1] - bingamma_m[i]) / bingamma_m[i] << " [%] "
                        << "between bin " << i - 1 << " and " << i << endl);
        }
    }
}


void PartBunchBase::calcGammas_cycl() {

    const int emittedBins = pbin_m->getLastemittedBin();

    for(int i = 0; i < emittedBins; i++)
        bingamma_m[i] = 0.0;
    for(unsigned int n = 0; n < getLocalNum(); n++)
        bingamma_m[this->Bin[n]] += sqrt(1.0 + dot(this->P[n], this->P[n]));
    for(int i = 0; i < emittedBins; i++) {
        reduce(bingamma_m[i], bingamma_m[i], OpAddAssign());
        if(pbin_m->getTotalNumPerBin(i) > 0)
            bingamma_m[i] /= pbin_m->getTotalNumPerBin(i);
        else
            bingamma_m[i] = 0.0;
        INFOMSG("Bin " << i << " : particle number = " << pbin_m->getTotalNumPerBin(i) << " gamma = " << bingamma_m[i] << endl);
    }

}


inline
void PartBunchBase::setBinCharge(int bin, double q) {
    this->Q = where(eq(this->Bin, bin), q, 0.0);
}

inline
void PartBunchBase::setBinCharge(int bin) {
    this->Q = where(eq(this->Bin, bin), this->Q, 0.0);
}

inline
std::pair<Vector_t, Vector_t> PartBunchBase::getEExtrema() {
    const Vector_t maxE = max(eg_m);
    //      const double maxL = max(dot(eg_m,eg_m));
    const Vector_t minE = min(eg_m);
    // INFOMSG("MaxE= " << maxE << " MinE= " << minE << endl);
    return std::pair<Vector_t, Vector_t>(maxE, minE);
}


size_t PartBunchBase::calcNumPartsOutside(Vector_t x) {

    partPerNode_m[Ippl::myNode()] = 0;
    const Vector_t meanR = get_rmean();

    for(unsigned long k = 0; k < getLocalNum(); ++ k)
        if (abs(R[k](0) - meanR(0)) > x(0) ||
            abs(R[k](1) - meanR(1)) > x(1) ||
            abs(R[k](2) - meanR(2)) > x(2)) {

            ++ partPerNode_m[Ippl::myNode()];
        }

    reduce(partPerNode_m.get(), partPerNode_m.get() + Ippl::getNodes(), globalPartPerNode_m.get(), OpAddAssign());

    return *globalPartPerNode_m.get();
}


size_t PartBunchBase::destroyT() {

    const unsigned int minNumParticlesPerCore = bunch_mp->getMinimumNumberOfParticlesPerCore();
    std::unique_ptr<size_t[]> tmpbinemitted;

    const size_t localNum = bunch_mp->getLocalNum();
    const size_t totalNum = bunch_mp->getTotalNum();

    if(weHaveEnergyBins()) {
        tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[getNumberOfEnergyBins()]);
        for(int i = 0; i < getNumberOfEnergyBins(); i++) {
            tmpbinemitted[i] = 0;
        }
        for(unsigned int i = 0; i < localNum; i++) {
            if (Bin[i] < 0) {
                bunch_mp->destroy(1, i);
            } else
                tmpbinemitted[Bin[i]]++;
        }
    } else {
        Inform dmsg("destroy: ", INFORM_ALL_NODES);
        size_t ne = 0;
        for(size_t i = 0; i < localNum; i++) {
            if((Bin[i] < 0) && ((localNum - ne) > minNumParticlesPerCore)) {   // need in minimum x particles per node
                ne++;
                bunch_mp->destroy(1, i);
            }
        }
        lowParticleCount_m = ((localNum - ne) <= minNumParticlesPerCore);
        reduce(lowParticleCount_m, lowParticleCount_m, OpOr());
        // unsigned int i = 0, j = localNum;
        // while (j > 0 && Bin[j - 1] < 0) -- j;

        // while (i + 1 < j) {
        //     if (Bin[i] < 0) {
        //         this->swap(i,j - 1);
        //         -- j;

        //         while (i + 1 < j && Bin[j - 1] < 0) -- j;
        //     }
        //     ++ i;
        // }

        // j = std::max(j, minNumParticlesPerCore);
        // for(unsigned int i = localNum; i > j; -- i) {
        //     destroy(1, i-1, true);
        // }
        // lowParticleCount_m = (j == minNumParticlesPerCore);
        // reduce(lowParticleCount_m, lowParticleCount_m, OpOr());

        if (ne > 0) {
            bunch_mp->performDestroy(true);
        }
    }

    calcBeamParameters();
    gatherLoadBalanceStatistics();

    if(weHaveEnergyBins()) {
        const int lastBin = dist_m->getLastEmittedEnergyBin();
        for(int i = 0; i < lastBin; i++) {
            binemitted_m[i] = tmpbinemitted[i];
        }
    }
    size_t newTotalNum = getLocalNum();
    reduce(newTotalNum, newTotalNum, OpAddAssign());

    setTotalNum(newTotalNum);

    return totalNum - newTotalNum;
}

double PartBunchBase::getPx(int i) {
    return 0.0;
}

double PartBunchBase::getPy(int i) {
    return 0.0;
}

double PartBunchBase::getPz(int i) {
    return 0.0;
}

//ff
double PartBunchBase::getX(int i) {
    return bunch_mp->R[i](0);
}

//ff
double PartBunchBase::getY(int i) {
    return bunch_mp->R[i](1);
}

//ff
double PartBunchBase::getX0(int i) {
    return 0.0;
}

//ff
double PartBunchBase::getY0(int i) {
    return 0.0;
}

//ff
double PartBunchBase::getZ(int i) {
    return bunch_mp->R[i](2);
}


inline
void PartBunch::setZ(int i, double zcoo) {};


inline
void PartBunch::get_bounds(Vector_t &rmin, Vector_t &rmax) {
    bounds(this->R, rmin, rmax);
}


void PartBunchBase::getLocalBounds(Vector_t &rmin, Vector_t &rmax) {
    const size_t localNum = bunch_mp->getLocalNum();
    if (localNum == 0) return;

    rmin = bunch_mp->R[0];
    rmax = bunch_mp->R[0];
    for (size_t i = 1; i < localNum; ++ i) {
        for (unsigned short d = 0; d < 3u; ++ d) {
            if (rmin(d) > bunch_mp->R[i](d)) rmin(d) = bunch_mp->R[i](d);
            else if (rmax(d) < bunch_mp->R[i](2)) rmax(d) = bunch_mp->R[i](d);
        }
    }
}


inline
std::pair<Vector_t, double> PartBunchBase::getBoundingSphere() {
    Vector_t rmin, rmax;
    get_bounds(rmin, rmax);

    std::pair<Vector_t, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


inline
std::pair<Vector_t, double> PartBunchBase::getLocalBoundingSphere() {
    Vector_t rmin, rmax;
    getLocalBounds(rmin, rmax);

    std::pair<Vector_t, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


inline
void PartBunchBase::push_back(Particle p) {
    Inform msg("PartBunch ");

    create(1);
    size_t i = bunch_mp->getTotalNum();

    R[i](0) = p[0];
    R[i](1) = p[1];
    R[i](2) = p[2];

    P[i](0) = p[3];
    P[i](1) = p[4];
    P[i](2) = p[5];

    bunch_mp->update();
    msg << "Created one particle i= " << i << endl;
}


inline
void PartBunchBase::set_part(FVector<double, 6> z, int ii) {
    R[ii](0) = z[0];
    P[ii](0) = z[1];
    R[ii](1) = z[2];
    P[ii](1) = z[3];
    R[ii](2) = z[4];
    P[ii](2) = z[5];
}

inline
void PartBunchBase::set_part(Particle p, int ii) {
    R[ii](0) = p[0];
    P[ii](0) = p[1];
    R[ii](1) = p[2];
    P[ii](1) = p[3];
    R[ii](2) = p[4];
    P[ii](2) = p[5];
}


inline
Particle PartBunchBase::get_part(int ii) {
    Particle part;
    part[0] = R[ii](0);
    part[1] = P[ii](0);
    part[2] = R[ii](1);
    part[3] = P[ii](1);
    part[4] = R[ii](2);
    part[5] = P[ii](2);
    return part;
}


void PartBunchBase::maximumAmplitudes(const FMatrix<double, 6, 6> &D,
                                      double &axmax, double &aymax) {
    axmax = aymax = 0.0;
    Particle part;

    for(unsigned int ii = 0; ii < getLocalNum(); ii++) {

        part = get_part(ii);

        double xnor =
            D(0, 0) * part.x()  + D(0, 1) * part.px() + D(0, 2) * part.y() +
            D(0, 3) * part.py() + D(0, 4) * part.t()  + D(0, 5) * part.pt();
        double pxnor =
            D(1, 0) * part.x()  + D(1, 1) * part.px() + D(1, 2) * part.y() +
            D(1, 3) * part.py() + D(1, 4) * part.t()  + D(1, 5) * part.pt();
        double ynor =
            D(2, 0) * part.x()  + D(2, 1) * part.px() + D(2, 2) * part.y() +
            D(2, 3) * part.py() + D(2, 4) * part.t()  + D(2, 5) * part.pt();
        double pynor =
            D(3, 0) * part.x()  + D(3, 1) * part.px() + D(3, 2) * part.y() +
            D(3, 3) * part.py() + D(3, 4) * part.t()  + D(3, 5) * part.pt();

        axmax = std::max(axmax, (xnor * xnor + pxnor * pxnor));
        aymax = std::max(aymax, (ynor * ynor + pynor * pynor));
    }
}


inline
void PartBunchBase::setdT(double dt) {
    dt_m = dt;
}

inline
double PartBunchBase::getdT() const {
    return dt_m;
}

inline
void PartBunchBase::setT(double t) {
    t_m = t;
}

inline
double PartBunchBase::getT() const {
    return t_m;
}

#endif
