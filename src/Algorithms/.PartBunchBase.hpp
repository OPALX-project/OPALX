//
// Class PartBunch
//   Base class for representing particle bunches.
//
// Copyright (c) 2008 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef PART_BUNCH_BASE_HPP
#define PART_BUNCH_BASE_HPP

#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/PartData.h"
#include "Distribution/Distribution.h"
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"
#include "Structure/FieldSolver.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <cmath>
#include <numeric>

extern Inform* gmsg;

template <class T, unsigned Dim>
PartBunch<T, Dim>::PartBunch(AbstractParticle<T, Dim>* pb, const PartData* ref)
    : R(*(pb->R_p)),
      ID(*(pb->ID_p)),
      pbin_m(nullptr),
      pmsg_m(nullptr),
      f_stream(nullptr),
      fixed_grid(false),
      reference(ref),
      unit_state_(units),
      stateOfLastBoundP_(unitless),
      moments_m(),
      dt_m(0.0),
      t_m(0.0),
      spos_m(0.0),
      globalMeanR_m(Vector_t<double, 3>(0.0, 0.0, 0.0)),
      globalToLocalQuaternion_m(Quaternion_t(1.0, 0.0, 0.0, 0.0)),
      rmax_m(0.0),
      rmin_m(0.0),
      hr_m(-1.0),
      nr_m(0),
      fs_m(nullptr),
      couplingConstant_m(0.0),
      qi_m(0.0),
      massPerParticle_m(0.0),
      distDump_m(0),
      dh_m(1e-12),
      bingamma_m(nullptr),
      binemitted_m(nullptr),
      stepsPerTurn_m(0),
      localTrackStep_m(0),
      globalTrackStep_m(0),
      numBunch_m(1),
      bunchTotalNum_m(1),
      bunchLocalNum_m(1),
      SteptoLastInj_m(0),
      globalPartPerNode_m(nullptr),
      dist_m(nullptr),
      dcBeam_m(false),
      periodLength_m(Physics::c / 1e9),
      pbase_m(pb)
{
    setup(pb);
}

/*
 * Bunch common member functions
 */

template <class T, unsigned Dim>
void PartBunch<T, Dim>::switchToUnitlessPositions(bool use_dt_per_particle) {

    bool hasToReset = false;
    if (!R.isDirty()) hasToReset = true;

    for (size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if (use_dt_per_particle)
            dt = this->dt[i];

        R[i] /= Vector_t<double, 3>(Physics::c * dt);
    }

    unit_state_ = unitless;

    if (hasToReset) R.resetDirtyFlag();
}

//FIXME: unify methods, use convention that all particles have own dt
template <class T, unsigned Dim>
void PartBunch<T, Dim>::switchOffUnitlessPositions(bool use_dt_per_particle) {

    bool hasToReset = false;
    if (!R.isDirty()) hasToReset = true;

    for (size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if (use_dt_per_particle)
            dt = this->dt[i];

        R[i] *= Vector_t<double, 3>(Physics::c * dt);
    }

    unit_state_ = units;

    if (hasToReset) R.resetDirtyFlag();
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::do_binaryRepart() {
 
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setDistribution(Distribution* d,
                                            std::vector<Distribution*> addedDistributions,
                                            size_t& np) {
    Inform m("setDistribution " );
    dist_m = d;
    dist_m->createOpalT(this, addedDistributions, np);
}

template <class T, unsigned Dim>
bool PartBunch<T, Dim>::isGridFixed() const {
    return fixed_grid;
}


template <class T, unsigned Dim>
bool PartBunch<T, Dim>::hasBinning() const {
    return (pbin_m != nullptr);
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::updateNumTotal() {
    size_t numParticles = getLocalNum();
    reduce(numParticles, numParticles, OpAddAssign());
    setTotalNum(numParticles);
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::rebin() {
    this->Bin = 0;
    pbin_m->resetBins();
    // delete pbin_m; we did not allocate it!
    pbin_m = nullptr;
}


template <class T, unsigned Dim>
int PartBunch<T, Dim>::getLastemittedBin() {
    if (pbin_m != nullptr)
        return pbin_m->getLastemittedBin();
    else
        return 0;
}

template <class T, unsigned Dim>
int PartBunch<T, Dim>::getNumberOfEnergyBins() {
    return 0;
}

template <class T, unsigned Dim>
size_t PartBunch<T, Dim>::emitParticles(double eZ) {
    return 0;
}

template <class T, unsigned Dim>
double PartBunch<T, Dim>::getEmissionDeltaT() {
    return 0.;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getBinGamma(int bin) {
    return 1.0;
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::Rebin() {
   
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::calcGammas() {
   
}

template <class T, unsigned Dim>
bool PartBunch<T, Dim>::weHaveEnergyBins() {
    return false;
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setLocalBinCount(size_t num, int bin) {
    if (pbin_m != nullptr) {
        pbin_m->setPartNum(bin, num);
    }
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setBinCharge(int bin, double q) {
  this->Q = where(eq(this->Bin, bin), q, 0.0);
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setBinCharge(int bin) {
  this->Q = where(eq(this->Bin, bin), this->qi_m, 0.0);
}


template <class T, unsigned Dim>
size_t PartBunch<T, Dim>::calcNumPartsOutside(Vector_t<double, 3> x) {

    std::size_t localnum = 0;
    const Vector_t<double, 3> meanR = get_rmean();

    for (unsigned long k = 0; k < getLocalNum(); ++ k)
        if (std::abs(R[k](0) - meanR(0)) > x(0) ||
            std::abs(R[k](1) - meanR(1)) > x(1) ||
            std::abs(R[k](2) - meanR(2)) > x(2)) {

            ++localnum;
        }

    gather(&localnum, &globalPartPerNode_m[0], 1);

    size_t npOutside = std::accumulate(globalPartPerNode_m.get(),
                                       globalPartPerNode_m.get() + Ippl::getNodes(), 0,
                                       std::plus<size_t>());

    return npOutside;
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
void PartBunch<T, Dim>::calcLineDensity(unsigned int nBins,
                                            std::vector<double>& lineDensity,
                                            std::pair<double, double>& meshInfo) {
    Vector_t<double, 3> rmin, rmax;
    get_bounds(rmin, rmax);

    if (nBins < 2) {
        Vektor<int, 3>/*NDIndex<3>*/ grid;
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
        double tau = (z - zmin) / hz - idx;

        lineDensity[idx] += Q[i] * (1.0 - tau) * perMeter;
        lineDensity[idx + 1] += Q[i] * tau * perMeter;
    }

    reduce(&(lineDensity[0]), &(lineDensity[0]) + nBins, &(lineDensity[0]), OpAddAssign());

    meshInfo.first = zmin;
    meshInfo.second = hz;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::boundp() {
    /*
      Assume rmin_m < 0.0
     */

    IpplTimings::startTimer(boundpTimer_m);
    //if (!R.isDirty() && stateOfLastBoundP_ == unit_state_) return;
    if ( !(R.isDirty() || ID.isDirty() ) && stateOfLastBoundP_ == unit_state_) return; //-DW

    stateOfLastBoundP_ = unit_state_;

    if (!isGridFixed()) {
        const int dimIdx = (dcBeam_m? 2: 3);

        /**
            In case of dcBeam_m && hr_m < 0
            this is the first call to boundp and we
            have to set hr completely i.e. x,y and z.
         */

        this->updateDomainLength(nr_m);
        IpplTimings::startTimer(boundpBoundsTimer_m);
        get_bounds(rmin_m, rmax_m);
        IpplTimings::stopTimer(boundpBoundsTimer_m);
        Vector_t<double, 3> len = rmax_m - rmin_m;

        double volume = 1.0;
        for (int i = 0; i < dimIdx; i++) {
            double length = std::abs(rmax_m[i] - rmin_m[i]);
            if (length < 1e-10) {
                rmax_m[i] += 1e-10;
                rmin_m[i] -= 1e-10;
            } else {
                rmax_m[i] += dh_m * length;
                rmin_m[i] -= dh_m * length;
            }
            hr_m[i]    = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
        }
        if (dcBeam_m) {
            rmax_m[2] = rmin_m[2] + periodLength_m;
            hr_m[2] = periodLength_m / (nr_m[2] - 1);
        }
        for (int i = 0; i < dimIdx; ++ i) {
            volume *= std::abs(rmax_m[i] - rmin_m[i]);
        }

        if (volume < 1e-21 && getTotalNum() > 1 && std::abs(sum(Q)) > 0.0) {
            WARNMSG(level1 << "!!! Extremely high particle density detected !!!" << endl);
        }
        //INFOMSG("It is a full boundp hz= " << hr_m << " rmax= " << rmax_m << " rmin= " << rmin_m << endl);

        if (hr_m[0] * hr_m[1] * hr_m[2] <= 0) {
            throw GeneralClassicException("boundp() ", "h<0, can not build a mesh");
        }

        Vector_t<double, 3> origin = rmin_m - Vector_t<double, 3>(hr_m[0] / 2.0, hr_m[1] / 2.0, hr_m[2] / 2.0);
        this->updateFields(hr_m, origin);
    }
    IpplTimings::startTimer(boundpUpdateTimer_m);
    update();
    IpplTimings::stopTimer(boundpUpdateTimer_m);
    R.resetDirtyFlag();

    IpplTimings::stopTimer(boundpTimer_m);
}


template <class T, unsigned Dim>
size_t PartBunch<T, Dim>::boundp_destroyT() {

    this->updateDomainLength(nr_m);

    std::vector<size_t> tmpbinemitted;

    boundp();

    size_t ne = 0;
    const size_t localNum = getLocalNum();

    double rzmean = 0.0; //momentsComputer_m.getMeanPosition()(2);
    double rzrms = 0.0; //momentsComputer_m.getStandardDeviationPosition()(2);
    const bool haveEnergyBins = weHaveEnergyBins();
    if (haveEnergyBins) {
        tmpbinemitted.resize(getNumberOfEnergyBins(), 0.0);
    }
    for (unsigned int i = 0; i < localNum; i++) {
        if (Bin[i] < 0 || (Options::remotePartDel > 0 && std::abs(R[i](2) - rzmean) < Options::remotePartDel * rzrms)) {
            ne++;
            destroy(1, i);
        } else if (haveEnergyBins) {
            tmpbinemitted[Bin[i]]++;
        }
    }

    boundp();

    calcBeamParameters();
    gatherLoadBalanceStatistics();

    reduce(ne, ne, OpAddAssign());
    return ne;
}


template <class T, unsigned Dim>
size_t PartBunch<T, Dim>::destroyT() {

    std::unique_ptr<size_t[]> tmpbinemitted;

    const size_t localNum = getLocalNum();
    const size_t totalNum = getTotalNum();
    size_t ne = 0;

    if (weHaveEnergyBins()) {
        tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[getNumberOfEnergyBins()]);
        for (int i = 0; i < getNumberOfEnergyBins(); i++) {
            tmpbinemitted[i] = 0;
        }
        for (unsigned int i = 0; i < localNum; i++) {
            if (Bin[i] < 0) {
                destroy(1, i);
                ++ ne;
            } else
                tmpbinemitted[Bin[i]]++;
        }
    } else {
        Inform dmsg("destroy: ", INFORM_ALL_NODES);
        for (size_t i = 0; i < localNum; i++) {
            if ((Bin[i] < 0)) {
                ne++;
                destroy(1, i);
            }
        }
    }

    if (ne > 0) {
        performDestroy(true);
    }

    calcBeamParameters();
    gatherLoadBalanceStatistics();

    size_t newTotalNum = getLocalNum();
    reduce(newTotalNum, newTotalNum, OpAddAssign());

    setTotalNum(newTotalNum);

    return totalNum - newTotalNum;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getPx(int /*i*/) {
    return 0.0;
}

template <class T, unsigned Dim>
double PartBunch<T, Dim>::getPy(int) {
    return 0.0;
}

template <class T, unsigned Dim>
double PartBunch<T, Dim>::getPz(int) {
    return 0.0;
}

template <class T, unsigned Dim>
double PartBunch<T, Dim>::getPx0(int) {
    return 0.0;
}

template <class T, unsigned Dim>
double PartBunch<T, Dim>::getPy0(int) {
    return 0;
}

//ff
template <class T, unsigned Dim>
double PartBunch<T, Dim>::getX(int i) {
    return this->R[i](0);
}

//ff
template <class T, unsigned Dim>
double PartBunch<T, Dim>::getY(int i) {
    return this->R[i](1);
}

//ff
template <class T, unsigned Dim>
double PartBunch<T, Dim>::getZ(int i) {
    return this->R[i](2);
}

//ff
template <class T, unsigned Dim>
double PartBunch<T, Dim>::getX0(int /*i*/) {
    return 0.0;
}

//ff
template <class T, unsigned Dim>
double PartBunch<T, Dim>::getY0(int /*i*/) {
    return 0.0;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setZ(int /*i*/, double /*zcoo*/) {
};


template <class T, unsigned Dim>
void PartBunch<T, Dim>::get_bounds(Vector_t<double, 3>& rmin, Vector_t<double, 3>& rmax) const {

    this->getLocalBounds(rmin, rmax);

    double min[2*Dim];

    for (unsigned int i = 0; i < Dim; ++i) {
        min[2*i] = rmin[i];
        min[2*i + 1] = -rmax[i];
    }

    allreduce(min, 2*Dim, std::less<double>());

    for (unsigned int i = 0; i < Dim; ++i) {
        rmin[i] = min[2*i];
        rmax[i] = -min[2*i + 1];
    }
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::getLocalBounds(Vector_t<double, 3>& rmin, Vector_t<double, 3>& rmax) const {
    const size_t localNum = getLocalNum();
    if (localNum == 0) {
        double maxValue = 1e8;
        rmin = Vector_t<double, 3>(maxValue, maxValue, maxValue);
        rmax = Vector_t<double, 3>(-maxValue, -maxValue, -maxValue);
        return;
    }

    rmin = R[0];
    rmax = R[0];
    for (size_t i = 1; i < localNum; ++ i) {
        for (unsigned short d = 0; d < 3u; ++ d) {
            if (rmin(d) > R[i](d)) rmin(d) = R[i](d);
            else if (rmax(d) < R[i](d)) rmax(d) = R[i](d);
        }
    }
}


template <class T, unsigned Dim>
std::pair<Vector_t<double, 3>, double> PartBunch<T, Dim>::getBoundingSphere() {
    Vector_t<double, 3> rmin, rmax;
    get_bounds(rmin, rmax);

    std::pair<Vector_t<double, 3>, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = std::sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


template <class T, unsigned Dim>
std::pair<Vector_t<double, 3>, double> PartBunch<T, Dim>::getLocalBoundingSphere() {
    Vector_t<double, 3> rmin, rmax;
    getLocalBounds(rmin, rmax);

    std::pair<Vector_t<double, 3>, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = std::sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::push_back(OpalParticle const& particle) {
    Inform msg("PartBunch ");

    size_t i = getLocalNum();
    create(1);

    R[i] = particle.getR();
    P[i] = particle.getP();

    update();
    msg << "Created one particle i= " << i << endl;
}



template <class T, unsigned Dim>
void PartBunch<T, Dim>::setdT(double dt) {
    dt_m = dt;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getdT() const {
    return dt_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setT(double t) {
    t_m = t;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::incrementT() {
    t_m += dt_m;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getT() const {
    return t_m;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::get_sPos() const {
    return spos_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::set_sPos(double s) {
    spos_m = s;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::get_gamma() const {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::get_meanKineticEnergy() const {
    return 0.0;
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_origin() const {
    return rmin_m;
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_maxExtent() const {
    return rmax_m;
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_centroid() const {
    return 0.0; 
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_rrms() const {
    return 0.0; 
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_rprms() const {
    return 0.0; 
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_rmean() const {
    return 0.0; 
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_prms() const {
    return 0.0; 
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_pmean() const {
    return 0.0;
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_pmean_Distribution() const {
    if (dist_m)// && dist_m->getType() != DistrTypeT::FROMFILE)
        return dist_m->get_pmean();

    double gamma = 0.1 / getM() + 1; // set default 0.1 eV
    return Vector_t<double, 3>(0, 0, std::sqrt(std::pow(gamma, 2) - 1));
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_emit() const {
    return 0.0;
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_norm_emit() const {
    return 0.0;
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_halo() const {
    return 0.0; 
}

template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_68Percentile() const {
    return 0.0;
}

template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_95Percentile() const {
    return 0.0; 
}

template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_99Percentile() const {
    return 0.0; 
}

template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_99_99Percentile() const {
    return 0.0;
}

template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_normalizedEps_68Percentile() const {
    return 0.0;
}

template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_normalizedEps_95Percentile() const {
    return 0.0;
}

template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_normalizedEps_99Percentile() const {
    return 0.0;
}

template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_normalizedEps_99_99Percentile() const {
    return 0.0;
}

template <class T, unsigned Dim>
double PartBunch<T, Dim>::get_Dx() const {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::get_Dy() const {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::get_DDx() const {
    return 0.0;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::get_DDy() const {
    return 0.0;
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::get_hr() const {
    return hr_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::set_meshEnlargement(double dh) {
    dh_m = dh;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::gatherLoadBalanceStatistics() {

    for (int i = 0; i < Ippl::getNodes(); i++)
        globalPartPerNode_m[i] = 0;

    std::size_t localnum = getLocalNum();
    gather(&localnum, &globalPartPerNode_m[0], 1);
}


template <class T, unsigned Dim>
size_t PartBunch<T, Dim>::getLoadBalance(int p) const {
    return globalPartPerNode_m[p];
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::get_PBounds(Vector_t<double, 3> &min, Vector_t<double, 3> &max) const {
    bounds(this->P, min, max);
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::calcBeamParameters() {

    IpplTimings::startTimer(statParamTimer_m);
    get_bounds(rmin_m, rmax_m);
    // momentsComputer_m.compute(*this);
    IpplTimings::stopTimer(statParamTimer_m);
}

template <class T, unsigned Dim>
double PartBunch<T, Dim>::getCouplingConstant() const {
    return couplingConstant_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setCouplingConstant(double c) {
    couplingConstant_m = c;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setCharge(double q) {
    qi_m = q;
    if (getTotalNum() != 0)
        Q = qi_m;
    else
        WARNMSG("Could not set total charge in PartBunch::setCharge based on getTotalNum" << endl);
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setChargeZeroPart(double q) {
    qi_m = q;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setMass(double mass) {
    massPerParticle_m = mass;
    if (getTotalNum() != 0)
        M = mass;
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setMassZeroPart(double mass) {
    massPerParticle_m = mass;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getCharge() const {
    return sum(Q);
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getChargePerParticle() const {
    return qi_m;
}

template <class T, unsigned Dim>
double PartBunch<T, Dim>::getMassPerParticle() const {
    return massPerParticle_m;
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setSolver(FieldSolver *fs) {
    fs_m = fs;
    fs_m->initSolver(this);

    /**
       CAN not re-inizialize ParticleLayout
       this is an IPPL issue
     */
    if (!OpalData::getInstance()->hasBunchAllocated()) {
        this->initialize(fs_m->getFieldLayout());
//         this->setMesh(fs_m->getMesh());
//         this->setFieldLayout(fs_m->getFieldLayout());
    }
}


template <class T, unsigned Dim>
bool PartBunch<T, Dim>::hasFieldSolver() {
    if (fs_m)
        return fs_m->hasValidSolver();
    else
        return false;
}


/// \brief Return the fieldsolver type if we have a fieldsolver
template <class T, unsigned Dim>
FieldSolverType PartBunch<T, Dim>::getFieldSolverType() const {
    if (fs_m) {
        return fs_m->getFieldSolverType();
    } else {
        return FieldSolverType::NONE;
    }
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setStepsPerTurn(int n) {
    stepsPerTurn_m = n;
}


template <class T, unsigned Dim>
int PartBunch<T, Dim>::getStepsPerTurn() const {
    return stepsPerTurn_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setGlobalTrackStep(long long n) {
    globalTrackStep_m = n;
}


template <class T, unsigned Dim>
long long PartBunch<T, Dim>::getGlobalTrackStep() const {
    return globalTrackStep_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setLocalTrackStep(long long n) {
    localTrackStep_m = n;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::incTrackSteps() {
    globalTrackStep_m++; localTrackStep_m++;
}


template <class T, unsigned Dim>
long long PartBunch<T, Dim>::getLocalTrackStep() const {
    return localTrackStep_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setNumBunch(short n) {
    numBunch_m = n;
    bunchTotalNum_m.resize(n);
    bunchLocalNum_m.resize(n);
}


template <class T, unsigned Dim>
short PartBunch<T, Dim>::getNumBunch() const {
    return numBunch_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setGlobalMeanR(Vector_t<double, 3> globalMeanR) {
    globalMeanR_m = globalMeanR;
}


template <class T, unsigned Dim>
Vector_t<double, 3> PartBunch<T, Dim>::getGlobalMeanR() {
    return globalMeanR_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion) {

    globalToLocalQuaternion_m = globalToLocalQuaternion;
}


template <class T, unsigned Dim>
Quaternion_t PartBunch<T, Dim>::getGlobalToLocalQuaternion() {
    return globalToLocalQuaternion_m;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setSteptoLastInj(int n) {
    SteptoLastInj_m = n;
}


template <class T, unsigned Dim>
int PartBunch<T, Dim>::getSteptoLastInj() const {
    return SteptoLastInj_m;
}

template <class T, unsigned Dim>
double PartBunch<T, Dim>::getQ() const {
    return reference->getQ();
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getM() const {
    return reference->getM();
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getP() const {
    return reference->getP();
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getE() const {
    return reference->getE();
}


template <class T, unsigned Dim>
ParticleOrigin PartBunch<T, Dim>::getPOrigin() const {
    return refPOrigin_m;
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setPOrigin(ParticleOrigin origin) {
    refPOrigin_m = origin;
}


template <class T, unsigned Dim>
ParticleType PartBunch<T, Dim>::getPType() const {
    return refPType_m;
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setPType(const std::string& type) {
    refPType_m =  ParticleProperties::getParticleType(type);
}


template <class T, unsigned Dim>
DistributionType PartBunch<T, Dim>::getDistType() const {
    return dist_m->getType();
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::resetQ(double q)  {
    const_cast<PartData *>(reference)->setQ(q);
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::resetM(double m)  {
    const_cast<PartData *>(reference)->setM(m);
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getdE() const {
    return 0.0; 
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getInitialBeta() const {
    return reference->getBeta();
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getInitialGamma() const {
    return reference->getGamma();
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getGamma(int /*i*/) {
    return 0;
}


template <class T, unsigned Dim>
double PartBunch<T, Dim>::getBeta(int /*i*/) {
    return 0;
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::actT() {
    // do nothing here
};


template <class T, unsigned Dim>
const PartData* PartBunch<T, Dim>::getReference() const {
    return reference;
}



template <class T, unsigned Dim>
void PartBunch<T, Dim>::calcEMean() {
    //momentsComputer_m.computeMeanKineticEnergy(*this);
}

template <class T, unsigned Dim>
Inform& PartBunch<T, Dim>::print(Inform& os) {

    if (getTotalNum() != 0) {  // to suppress Nans
        Inform::FmtFlags_t ff = os.flags();

        double pathLength = get_sPos();

        os << std::scientific;
        os << level1 << "\n";
        os << "* ************** B U N C H ********************************************************* \n";
        os << "* NP              = " << getTotalNum() << "\n";
        os << "* Qtot            = " << std::setw(17) << Util::getChargeString(std::abs(sum(Q))) << "         "
        << "Qi    = "             << std::setw(17) << Util::getChargeString(std::abs(qi_m)) << "\n";
        os << "* Ekin            = " << std::setw(17) << Util::getEnergyString(get_meanKineticEnergy()) << "         "
           << "dEkin = "             << std::setw(17) << Util::getEnergyString(getdE()) << "\n";
        os << "* rmax            = " << Util::getLengthString(rmax_m, 5) << "\n";
        os << "* rmin            = " << Util::getLengthString(rmin_m, 5) << "\n";
        if (getTotalNum() >= 2) { // to suppress Nans
            os << "* rms beam size   = " << Util::getLengthString(get_rrms(), 5) << "\n";
            os << "* rms momenta     = " << std::setw(12) << std::setprecision(5) << get_prms() << " [beta gamma]\n";
            os << "* mean position   = " << Util::getLengthString(get_rmean(), 5) << "\n";
            os << "* mean momenta    = " << std::setw(12) << std::setprecision(5) << get_pmean() << " [beta gamma]\n";
            os << "* rms emittance   = " << std::setw(12) << std::setprecision(5) << get_emit() << " (not normalized)\n";
            os << "* rms correlation = " << std::setw(12) << std::setprecision(5) << get_rprms() << "\n";
        }
        os << "* hr              = " << Util::getLengthString(get_hr(), 5) << "\n";
        os << "* dh              = " << std::setw(13) << std::setprecision(5) << dh_m * 100 << " [%]\n";
        os << "* t               = " << std::setw(17) << Util::getTimeString(getT()) << "         "
           << "dT    = "             << std::setw(17) << Util::getTimeString(getdT()) << "\n";
        os << "* spos            = " << std::setw(17) << Util::getLengthString(pathLength) << "\n";
        os << "* ********************************************************************************** " << endl;
        os.flags(ff);
    }
    return os;
}

// angle range [0~2PI) degree
template <class T, unsigned Dim>
double PartBunch<T, Dim>::calculateAngle(double x, double y) {
    double thetaXY = atan2(y, x);

    return thetaXY >= 0 ? thetaXY : thetaXY + Physics::two_pi;
}


template <class T, unsigned Dim>
Inform& operator<<(Inform &os, PartBunch<T, Dim>& p) {
    return p.print(os);
}


/*
 * Virtual member functions
 */

template <class T, unsigned Dim>
void PartBunch<T, Dim>::runTests() {
    throw OpalException("PartBunch<T, Dim>::runTests() ", "No test supported.");
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::resetInterpolationCache(bool /*clearCache*/) {

}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::swap(unsigned int i, unsigned int j) {
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
    std::swap(POrigin[i], POrigin[j]);
    std::swap(TriID[i], TriID[j]);
    std::swap(cavityGapCrossed[i], cavityGapCrossed[j]);
    std::swap(bunchNum[i], bunchNum[j]);
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setBCAllPeriodic() {
    throw OpalException("PartBunch<T, Dim>::setBCAllPeriodic() ", "Not supported BC.");
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setBCAllOpen() {
    throw OpalException("PartBunch<T, Dim>::setBCAllOpen() ", "Not supported BC.");
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setBCForDCBeam() {
    throw OpalException("PartBunch<T, Dim>::setBCForDCBeam() ", "Not supported BC.");
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::updateFields(const Vector_t<double, 3>& /*hr*/, const Vector_t<double, 3>& /*origin*/) {
}


template <class T, unsigned Dim>
void PartBunch<T, Dim>::setup(AbstractParticle<T, Dim>* pb) {

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
    pb->addAttribute(POrigin);
    pb->addAttribute(TriID);
    pb->addAttribute(cavityGapCrossed);
    pb->addAttribute(bunchNum);

    boundpTimer_m       = IpplTimings::getTimer("Boundingbox");
    boundpBoundsTimer_m = IpplTimings::getTimer("Boundingbox-bounds");
    boundpUpdateTimer_m = IpplTimings::getTimer("Boundingbox-update");
    statParamTimer_m    = IpplTimings::getTimer("Compute Statistics");
    selfFieldTimer_m    = IpplTimings::getTimer("SelfField total");

    histoTimer_m        = IpplTimings::getTimer("Histogram");

    distrCreate_m       = IpplTimings::getTimer("Create Distr");
    distrReload_m       = IpplTimings::getTimer("Load Distr");

    globalPartPerNode_m = std::unique_ptr<size_t[]>(new size_t[Ippl::getNodes()]);

    pmsg_m.release();
}


template <class T, unsigned Dim>
size_t PartBunch<T, Dim>::getTotalNum() const {
    return pbase_m->getTotalNum();
}

template <class T, unsigned Dim>
size_t PartBunch<T, Dim>::getLocalNum() const {
    return pbase_m->getLocalNum();
}


template <class T, unsigned Dim>
size_t PartBunch<T, Dim>::getDestroyNum() const {
    return pbase_m->getDestroyNum();
}

template <class T, unsigned Dim>
size_t PartBunch<T, Dim>::getGhostNum() const {
    return pbase_m->getGhostNum();
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setTotalNum(size_t n) {
    pbase_m->setTotalNum(n);
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setLocalNum(size_t n) {
    pbase_m->setLocalNum(n);
}

template <class T, unsigned Dim>
ParticleLayout<T, Dim> & PartBunch<T, Dim>::getLayout() {
    return pbase_m->getLayout();
}

template <class T, unsigned Dim>
const ParticleLayout<T, Dim>& PartBunch<T, Dim>::getLayout() const {
    return pbase_m->getLayout();
}

template <class T, unsigned Dim>
bool PartBunch<T, Dim>::getUpdateFlag(UpdateFlags_t f) const {
    return pbase_m->getUpdateFlag(f);
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setUpdateFlag(UpdateFlags_t f, bool val) {
    pbase_m->setUpdateFlag(f, val);
}

template <class T, unsigned Dim>
bool PartBunch<T, Dim>::singleInitNode() const {
    return pbase_m->singleInitNode();
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::resetID() {
    pbase_m->resetID();
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::update() {
    try {
        pbase_m->update();
    } catch (const IpplException& ex) {
        throw OpalException(ex.where(), ex.what());
    }
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::update(const ParticleAttrib<char>& canSwap) {
    try {
        pbase_m->update(canSwap);
    } catch (const IpplException& ex) {
        throw OpalException(ex.where(), ex.what());
    }
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::createWithID(unsigned id) {
    pbase_m->createWithID(id);
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::create(size_t M) {
    pbase_m->create(M);
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::globalCreate(size_t np) {
    pbase_m->globalCreate(np);
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::destroy(size_t M, size_t I, bool doNow) {
    pbase_m->destroy(M, I, doNow);
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::performDestroy(bool updateLocalNum) {
    pbase_m->performDestroy(updateLocalNum);
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::ghostDestroy(size_t M, size_t I) {
    pbase_m->ghostDestroy(M, I);
}

template <class T, unsigned Dim>
void PartBunch<T, Dim>::setBeamFrequency(double f) {
    periodLength_m = Physics::c / f;
}

template <class T, unsigned Dim>
matrix_t PartBunch<T, Dim>::getSigmaMatrix() const {

}

#endif
