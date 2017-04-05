#ifndef PART_BUNCH_BASE_HPP
#define PART_BUNCH_BASE_HPP

#include "PartBunchBase.h"

bool PartBunchBase::getIfBeamEmitting() {
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


int PartBunchBase::getLastEmittedEnergyBin() {
    /*
     * Get maximum last emitted energy bin.
     */
    int lastEmittedBin = dist_m->getLastEmittedEnergyBin();
    reduce(lastEmittedBin, lastEmittedBin, OpMaxAssign());
    return lastEmittedBin;
}


size_t PartBunchBase::getNumberOfEmissionSteps() {
    return dist_m->getNumberOfEmissionSteps();
}


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
