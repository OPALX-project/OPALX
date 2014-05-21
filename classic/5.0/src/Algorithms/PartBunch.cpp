// ------------------------------------------------------------------------
// $RCSfile: PartBunch.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class PartBunch
//   Interface to a particle bunch.
//   Can be used to avoid use of a template in user code.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/PartBunch.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include <iostream>
#include <cfloat>
#include <fstream>
#include <iomanip>

#include "AbstractObjects/OpalData.h"
#include "Distribution/Distribution.h"
#include "Structure/LossDataSink.h"
#include "Structure/FieldSolver.h"
#include "Utilities/Options.h"

#include "ListElem.h"
#include "BasicActions/Option.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_qrng.h>

#ifdef OPAL_NOCPLUSPLUS11_NULLPTR
#define nullptr NULL
#endif

using Physics::pi;

using namespace std;

extern Inform *gmsg;

// Class PartBunch
// ------------------------------------------------------------------------

PartBunch::PartBunch(const PartData *ref):
    myNode_m(Ippl::myNode()),
    nodes_m(Ippl::getNodes()),
    fixed_grid(false),
    pbin_m(nullptr),
    lossDs_m(nullptr),
    pmsg_m(nullptr),
    f_stream(nullptr),
    reference(ref),
    unit_state_(units),
    stateOfLastBoundP_(unitless),
    lineDensity_m(nullptr),
    nBinsLineDensity_m(0),
    moments_m(),
    dt_m(0.0),
    t_m(0.0),
    eKin_m(0.0),
    energy_m(nullptr),
    dE_m(0.0),
    rmax_m(0.0),
    rmin_m(0.0),
    rrms_m(0.0),
    prms_m(0.0),
    rmean_m(0.0),
    pmean_m(0.0),
    eps_m(0.0),
    eps_norm_m(0.0),
    rprms_m(0.0),
    Dx_m(0.0),
    Dy_m(0.0),
    DDx_m(0.0),
    DDy_m(0.0),
    hr_m(-1.0),
    nr_m(0),
    fs_m(nullptr),
    couplingConstant_m(0.0),
    qi_m(0.0),
    interpolationCacheSet_m(false),
    distDump_m(0),
    stash_Nloc_m(0),
    stash_iniR_m(0.0),
    stash_iniP_m(0.0),
    bunchStashed_m(false),
    fieldDBGStep_m(0),
    dh_m(0.0),
    tEmission_m(0.0),
    bingamma_m(nullptr),
    binemitted_m(nullptr),
    lPath_m(0.0),
    stepsPerTurn_m(0),
    localTrackStep_m(0),
    globalTrackStep_m(0),
    numBunch_m(1),
    SteptoLastInj_m(0),
    partPerNode_m(nullptr),
    globalPartPerNode_m(nullptr),
    dist_m(nullptr),
    dcBeam_m(false) {
    addAttribute(X);
    addAttribute(P);
    addAttribute(Q);
    addAttribute(M);
    addAttribute(Ef);
    addAttribute(Eftmp);

    addAttribute(Bf);
    addAttribute(Bin);
    addAttribute(dt);
    addAttribute(LastSection);
    addAttribute(PType);
    addAttribute(TriID);

    selfFieldTimer_m = IpplTimings::getTimer("SelfField");
    boundpTimer_m = IpplTimings::getTimer("Boundingbox");
    statParamTimer_m = IpplTimings::getTimer("Statistics");
    compPotenTimer_m  = IpplTimings::getTimer("Potential");

    histoTimer_m = IpplTimings::getTimer("Histogram");

    distrCreate_m = IpplTimings::getTimer("CreatDistr");
    distrReload_m = IpplTimings::getTimer("LoadDistr");


    partPerNode_m = std::unique_ptr<double[]>(new double[Ippl::getNodes()]);
    globalPartPerNode_m = std::unique_ptr<double[]>(new double[Ippl::getNodes()]);

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(std::string("GlobalLosses"), !Options::asciidump));

    pmsg_m.release();
    //    f_stream.release();
    /*
    if(Ippl::getNodes() == 1) {
        f_stream = std::unique_ptr<ofstream>(new ofstream);
        f_stream->open("data/dist.dat", ios::out);
        pmsg_m = std::unique_ptr<Inform>(new Inform(0, *f_stream, 0));
    }
    */
}

PartBunch::PartBunch(const PartBunch &rhs):
    myNode_m(Ippl::myNode()),
    nodes_m(Ippl::getNodes()),
    fixed_grid(rhs.fixed_grid),
    pbin_m(nullptr),
    lossDs_m(nullptr),
    pmsg_m(nullptr),
    f_stream(nullptr),
    reference(rhs.reference),
    unit_state_(rhs.unit_state_),
    stateOfLastBoundP_(rhs.stateOfLastBoundP_),
    lineDensity_m(nullptr),
    nBinsLineDensity_m(rhs.nBinsLineDensity_m),
    moments_m(rhs.moments_m),
    dt_m(rhs.dt_m),
    t_m(rhs.t_m),
    eKin_m(rhs.eKin_m),
    energy_m(nullptr),
    dE_m(rhs.dE_m),
    rmax_m(rhs.rmax_m),
    rmin_m(rhs.rmin_m),
    rrms_m(rhs.rrms_m),
    prms_m(rhs.prms_m),
    rmean_m(rhs.rmean_m),
    pmean_m(rhs.pmean_m),
    eps_m(rhs.eps_m),
    eps_norm_m(rhs.eps_norm_m),
    rprms_m(rhs.rprms_m),
    Dx_m(rhs.Dx_m),
    Dy_m(rhs.Dy_m),
    DDx_m(rhs.DDx_m),
    DDy_m(rhs.DDy_m),
    hr_m(rhs.hr_m),
    nr_m(rhs.nr_m),
    fs_m(nullptr),
    couplingConstant_m(rhs.couplingConstant_m),
    qi_m(rhs.qi_m),
    interpolationCacheSet_m(rhs.interpolationCacheSet_m),
    distDump_m(rhs.distDump_m),
    stash_Nloc_m(rhs.stash_Nloc_m),
    stash_iniR_m(rhs.stash_iniR_m),
    stash_iniP_m(rhs.stash_iniP_m),
    bunchStashed_m(rhs.bunchStashed_m),
    fieldDBGStep_m(rhs.fieldDBGStep_m),
    dh_m(rhs.dh_m),
    tEmission_m(rhs.tEmission_m),
    bingamma_m(nullptr),
    binemitted_m(nullptr),
    lPath_m(rhs.lPath_m),
    stepsPerTurn_m(rhs.stepsPerTurn_m),
    localTrackStep_m(rhs.localTrackStep_m),
    globalTrackStep_m(rhs.globalTrackStep_m),
    numBunch_m(rhs.numBunch_m),
    SteptoLastInj_m(rhs.SteptoLastInj_m),
    partPerNode_m(nullptr),
    globalPartPerNode_m(nullptr),
    dist_m(nullptr),
    dcBeam_m(rhs.dcBeam_m) {
    ERRORMSG("should not be here: PartBunch::PartBunch(const PartBunch &rhs):" << endl);
    std::exit(0);
}


PartBunch::PartBunch(const std::vector<Particle> &rhs, const PartData *ref):
    myNode_m(Ippl::myNode()),
    nodes_m(Ippl::getNodes()),
    fixed_grid(false),
    pbin_m(nullptr),
    lossDs_m(nullptr),
    pmsg_m(nullptr),
    f_stream(nullptr),
    reference(ref),
    unit_state_(units),
    stateOfLastBoundP_(unitless),
    lineDensity_m(nullptr),
    nBinsLineDensity_m(0),
    moments_m(),
    dt_m(0.0),
    t_m(0.0),
    eKin_m(0.0),
    energy_m(nullptr),
    dE_m(0.0),
    rmax_m(0.0),
    rmin_m(0.0),
    rrms_m(0.0),
    prms_m(0.0),
    rmean_m(0.0),
    pmean_m(0.0),
    eps_m(0.0),
    eps_norm_m(0.0),
    rprms_m(0.0),
    Dx_m(0.0),
    Dy_m(0.0),
    DDx_m(0.0),
    DDy_m(0.0),
    hr_m(-1.0),
    nr_m(0),
    fs_m(nullptr),
    couplingConstant_m(0.0),
    qi_m(0.0),
    interpolationCacheSet_m(false),
    distDump_m(0),
    stash_Nloc_m(0),
    stash_iniR_m(0.0),
    stash_iniP_m(0.0),
    bunchStashed_m(false),
    fieldDBGStep_m(0),
    dh_m(0.0),
    tEmission_m(0.0),
    bingamma_m(nullptr),
    binemitted_m(nullptr),
    lPath_m(0.0),
    stepsPerTurn_m(0),
    localTrackStep_m(0),
    globalTrackStep_m(0),
    numBunch_m(1),
    SteptoLastInj_m(0),
    partPerNode_m(nullptr),
    globalPartPerNode_m(nullptr),
    dist_m(nullptr),
    dcBeam_m(false) {
    ERRORMSG("should not be here: PartBunch::PartBunch(const std::vector<Particle> &rhs, const PartData *ref):" << endl);
}

PartBunch::~PartBunch() {

}

/// \brief make density histograms
void PartBunch::makHistograms()  {
    IpplTimings::startTimer(histoTimer_m);
    const unsigned int bins = 1000;
    if(getTotalNum() > bins) {
        int tag = Ippl::Comm->next_tag(IPPL_APP_TAG1, IPPL_APP_CYCLE);
        gsl_histogram *h = gsl_histogram_alloc(bins);
        const double l = rmax_m[2] - rmin_m[2]; // max => min
        gsl_histogram_set_ranges_uniform(h, 0.0, l);
        const double minz = abs(rmin_m[2]);

        // 1d hitogram z positions
        for(size_t n = 0; n < getLocalNum(); n++)
            gsl_histogram_increment(h, R[n](2) - minz);

        // now we need to reduce

        if(Ippl::myNode() == 0) {
            // wait for msg from all processors (EXEPT NODE 0)
            int notReceived = Ippl::getNodes() - 1;
            double recVal = 0;
            while(notReceived > 0) {
                int node = COMM_ANY_NODE;
                std::unique_ptr<Message> rmsg(Ippl::Comm->receive_block(node, tag));
                if(!bool(rmsg))
                    ERRORMSG("Could not receive from client nodes in makHistograms." << endl);
                for(unsigned int i = 0; i < bins; i++) {
                    rmsg->get(&recVal);
                    gsl_histogram_increment(h, recVal);
                }
                notReceived--;
            }
            stringstream filename_str;
            static unsigned int file_number = 0;
            ++ file_number;
            filename_str << "data/zhist-" << file_number << ".dat";
            FILE *fp;
            fp = fopen(filename_str.str().c_str(), "w");
            gsl_histogram_fprintf(fp, h, "%g", "%g");
            fclose(fp);
        } else {
            Message *smsg = new Message();
            for(unsigned int i = 0; i < bins; i++)
                smsg->put(gsl_histogram_get(h, i));
            bool res = Ippl::Comm->send(smsg, 0, tag);
            if(! res)
                ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
        }
        gsl_histogram_free(h);
    }
    IpplTimings::stopTimer(histoTimer_m);
}


/// \brief Need Ek for the Schottky effect calculation (eV)
double PartBunch::getEkin() const {
    if(dist_m)
        return dist_m->GetEkin();
    else
        return 0.0;
}

/// \brief Need the work function for the Schottky effect calculation (eV)
double PartBunch::getWorkFunctionRf() const {
    if(dist_m)
        return dist_m->GetWorkFunctionRf();
    else
        return 0.0;
}
/// \brief Need the laser energy for the Schottky effect calculation (eV)
double PartBunch::getLaserEnergy() const {
    if(dist_m)
        return dist_m->GetLaserEnergy();
    else
        return 0.0;
}



/** \brief After each Schottky scan we delete all the particles.

 */
void PartBunch::cleanUpParticles() {

    size_t np = getTotalNum();
    bool scan = false;

    destroy(getLocalNum(), 0, true);

    dist_m->CreateOpalT(*this, np, scan);

    update();
}

void PartBunch::setDistribution(Distribution *d,
                                std::vector<Distribution *> addedDistributions,
                                size_t &np,
                                bool scan) {
    Inform m("setDistribution " );
    dist_m = d;
    dist_m->CreateOpalT(*this, addedDistributions, np, scan);
//    if (Options::cZero)
//        dist_m->Create(*this, addedDistributions, np / 2, scan);
//    else
//        dist_m->Create(*this, addedDistributions, np, scan);
}

void PartBunch::resetIfScan()
/*
  In case of a scan we have
  to reset some variables
 */
{
    dt = 0.0;
    localTrackStep_m = 0;
}



bool PartBunch::hasFieldSolver() {
    if(fs_m)
        return fs_m->hasValidSolver();
    else
        return false;
}


bool PartBunch::hasZeroNLP() {
    /**
       Check if a node has no particles
     */
    Inform m("hasZeroNLP() ", INFORM_ALL_NODES);
    int minnlp = 0;
    int nlp = getLocalNum();
    minnlp = 100000;
    reduce(nlp, minnlp, OpMinAssign());
    return (minnlp == 0);
}

double PartBunch::getPx(int i) {
    return 0.0;
}

double PartBunch::getPy(int i) {
    return 0.0;
}

double PartBunch::getPz(int i) {
    return 0.0;
}

//ff
double PartBunch::getX(int i) {
    return this->R[i](0);
}

//ff
double PartBunch::getY(int i) {
    return this->R[i](1);
}

//ff
double PartBunch::getX0(int i) {
    return 0.0;
}

//ff
double PartBunch::getY0(int i) {
    return 0.0;
}

//ff
double PartBunch::getZ(int i) {
    return this->R[i](2);
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
void PartBunch::calcLineDensity() {
    //   e_dim_tag decomp[3];
    list<ListElem> listz;

    //   for (int d=0; d < 3; ++d) {                                    // this does not seem to work properly
    //     decomp[d] = getFieldLayout().getRequestedDistribution(d);
    //   }

    FieldLayout_t &FL  = getFieldLayout();
    double hz = getMesh().get_meshSpacing(2); // * Physics::c * getdT();
    //   FieldLayout_t *FL  = new FieldLayout_t(getMesh(), decomp);

    if(!bool(lineDensity_m)) {
        if(nBinsLineDensity_m == 0)
            nBinsLineDensity_m = nr_m[2];
        lineDensity_m = std::unique_ptr<double[]>(new double[nBinsLineDensity_m]);
    }

    for(unsigned int i = 0; i < nBinsLineDensity_m; ++i)
        lineDensity_m[i] = 0.0;

    rho_m = 0.0;
    this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

    //   NDIndex<Dim> idx = FL->getLocalNDIndex(); // gives the wrong indices!!
    //   NDIndex<Dim> idxdom = FL->getDomain();
    NDIndex<Dim> idx = FL.getLocalNDIndex();
    NDIndex<Dim> idxdom = FL.getDomain();
    NDIndex<Dim> elem;
    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG1, IPPL_APP_CYCLE);
    double spos = get_sPos();
    double T = getT();

    if(Ippl::myNode() == 0) {
        for(int i = idx[2].min(); i <= idx[2].max(); ++i) {
            double localquantsum = 0.0;
            elem[2] = Index(i, i);
            for(int j = idx[1].min(); j <= idx[1].max(); ++j) {
                elem[1] = Index(j, j);
                for(int k = idx[0].min(); k <= idx[0].max(); ++k) {
                    elem[0] = Index(k, k);
                    localquantsum += rho_m.localElement(elem) / hz;
                }
            }
            listz.push_back(ListElem(spos, T, i, i, localquantsum));
        }
        // wait for msg from all processors (EXEPT NODE 0)
        int notReceived = Ippl::getNodes() - 1;
        int dataBlocks = 0;
        int coor;
        double projVal;
        while(notReceived > 0) {
            int node = COMM_ANY_NODE;
            std::unique_ptr<Message> rmsg(Ippl::Comm->receive_block(node, tag));
            if(!bool(rmsg)) {
                ERRORMSG("Could not receive from client nodes in main." << endl);
            }
            notReceived--;
            rmsg->get(&dataBlocks);
            for(int i = 0; i < dataBlocks; i++) {
                rmsg->get(&coor);
                rmsg->get(&projVal);
                listz.push_back(ListElem(spos, T, coor, coor, projVal));
            }
        }
        listz.sort();
        /// copy line density in listz to array of double
        fillArray(lineDensity_m.get(), listz);
    } else {
        Message *smsg = new Message();
        smsg->put(idx[2].max() - idx[2].min() + 1);
        for(int i = idx[2].min(); i <= idx[2].max(); ++i) {
            double localquantsum = 0.0;
            elem[2] = Index(i, i);
            for(int j = idx[1].min(); j <= idx[1].max(); ++j) {
                elem[1] = Index(j, j);
                for(int k = idx[0].min(); k <= idx[0].max(); ++k) {
                    elem[0] = Index(k, k);
                    localquantsum +=  rho_m.localElement(elem) / hz;
                }
            }
            smsg->put(i);
            smsg->put(localquantsum);
        }
        bool res = Ippl::Comm->send(smsg, 0, tag);
        if(! res)
            ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
    }
    reduce(&(lineDensity_m[0]), &(lineDensity_m[0]) + nBinsLineDensity_m, &(lineDensity_m[0]), OpAddAssign());
}

void PartBunch::fillArray(double *lineDensity, const list<ListElem> &l) {
    unsigned int mmax = 0;
    unsigned int nmax = 0;
    unsigned int count = 0;

    for(list<ListElem>::const_iterator it = l.begin(); it != l.end() ; ++it)  {
        if(it->m > mmax) mmax = it->m;
        if(it->n > nmax) nmax = it->n;
    }
    for(list<ListElem>::const_iterator it = l.begin(); it != l.end(); ++it)
        if((it->m < mmax) && (it->n < nmax)) {
            lineDensity[count] = it->den;
            count++;
        }
}

void PartBunch::getLineDensity(vector<double> &lineDensity) {
    if(bool(lineDensity_m)) {
        if(lineDensity.size() != nBinsLineDensity_m)
            lineDensity.resize(nBinsLineDensity_m, 0.0);
        for(unsigned int i  = 0; i < nBinsLineDensity_m; ++i)
            lineDensity[i] = lineDensity_m[i];
    }
}

void PartBunch::updateBinStructure()
{ }

void PartBunch::calcGammas() {

    const int emittedBins = dist_m->GetNumberOfEnergyBins();

    size_t s = 0;

    for(int i = 0; i < emittedBins; i++)
        bingamma_m[i] = 0.0;

    for(unsigned int n = 0; n < getLocalNum(); n++)
        bingamma_m[this->Bin[n]] += sqrt(1.0 + dot(this->P[n], this->P[n]));

    for(int i = 0; i < emittedBins; i++) {
        reduce(bingamma_m[i], bingamma_m[i], OpAddAssign());

        size_t pInBin = (binemitted_m[i]);
        reduce(pInBin, pInBin, OpAddAssign());
        if(pInBin != 0) {
            bingamma_m[i] /= pInBin;
            INFOMSG("Bin " << i << " gamma = " << setw(8) << scientific << setprecision(5) << bingamma_m[i] << "; NpInBin= " << setw(8) << setfill(' ') << pInBin << endl);
        } else {
            bingamma_m[i] = 1.0;
            INFOMSG("Bin " << i << " has no particles " << endl);
        }
        s += pInBin;
    }

    
    
    if(s != getTotalNum() && !OpalData::getInstance()->hasGlobalGeometry())
        ERRORMSG("sum(Bins)= " << s << " != sum(R)= " << getTotalNum() << endl;);

    if(emittedBins >= 2) {
        for(int i = 1; i < emittedBins; i++) {
            if(binemitted_m[i - 1] != 0 && binemitted_m[i] != 0)
                INFOMSG("dE= " << getM() * 1.0E-3 * (bingamma_m[i - 1] - bingamma_m[i]) << " [keV] of Bin " << i - 1 << " and " << i << endl);
        }
    }
}


void PartBunch::calcGammas_cycl() {

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


double PartBunch::getMaxdEBins() {

    const int emittedBins = pbin_m->getLastemittedBin();

    double maxdE = DBL_MIN;
    double maxdEGlobal = DBL_MIN;
    if(emittedBins >= 1) {
        for(int i = 1; i < emittedBins; i++) {
            const size_t pInBin1 = (binemitted_m[i]);
            const size_t pInBin2 = (binemitted_m[i - 1]);
            if(pInBin1 != 0 && pInBin2 != 0) {
                double de = fabs(getM() * 1.0E-3 * (bingamma_m[i - 1] - bingamma_m[i]));
                if(de > maxdE)
                    maxdE = de;
            }
        }

        reduce(maxdE, maxdEGlobal, OpMaxAssign());

        return maxdEGlobal;
    } else
        return DBL_MAX;
}


void PartBunch::computeSelfFields(int binNumber) {
    IpplTimings::startTimer(selfFieldTimer_m);

    /// Set initial charge density to zero. Create image charge
    /// potential field.
    rho_m = 0.0;
    Field_t imagePotential = rho_m;

    /// Set initial E field to zero.
    eg_m = Vector_t(0.0);

    if(fs_m->hasValidSolver()) {
        /// Scatter charge onto space charge grid.
        this->Q *= this->dt;
        if(!interpolationCacheSet_m) {
            if(interpolationCache_m.size() < getLocalNum()) {
                interpolationCache_m.create(getLocalNum() - interpolationCache_m.size());
            } else {
                interpolationCache_m.destroy(interpolationCache_m.size() - getLocalNum(),
                                             getLocalNum(),
                                             true);
            }
            interpolationCacheSet_m = true;

            this->Q.scatter(this->rho_m, this->R, IntrplCIC_t(), interpolationCache_m);
        } else {
            this->Q.scatter(this->rho_m, IntrplCIC_t(), interpolationCache_m);
        }
        this->Q /= this->dt;
        this->rho_m /= getdT();

        /// Calculate mesh-scale factor and get gamma for this specific slice (bin).
        double scaleFactor = 1;
        // double scaleFactor = Physics::c * getdT();
        double gammaz = getBinGamma(binNumber);

        /// Scale charge density to get charge density in real units. Account for
        /// Lorentz transformation in longitudinal direction.
        double tmp2 = 1 / hr_m[0] * 1 / hr_m[1] * 1 / hr_m[2] / (scaleFactor * scaleFactor * scaleFactor) / gammaz;
        rho_m *= tmp2;

        /// Scale mesh spacing to real units (meters). Lorentz transform the
        /// longitudinal direction.
        Vector_t hr_scaled = hr_m * Vector_t(scaleFactor);
        hr_scaled[2] *= gammaz;

        /// Find potential from charge in this bin (no image yet) using Poisson's
        /// equation (without coefficient: -1/(eps)).
        IpplTimings::startTimer(compPotenTimer_m);
        imagePotential = rho_m;

        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        IpplTimings::stopTimer(compPotenTimer_m);

        /// Scale mesh back (to same units as particle locations.)
        rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// The scalar potential is given back in rho_m
        /// and must be converted to the right units.
        rho_m *= getCouplingConstant();

        /// IPPL Grad numerical computes gradient to find the
        /// electric field (in bin rest frame).
        eg_m = -Grad(rho_m, eg_m);

        /// Scale field. Combine Lorentz transform with conversion to proper field
        /// units.
        eg_m *= Vector_t(gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz));

        /// Interpolate electric field at particle positions.  We reuse the
        /// cached information about where the particles are relative to the
        /// field, since the particles have not moved since this the most recent
        /// scatter operation.
        Eftmp.gather(eg_m, IntrplCIC_t(), interpolationCache_m);
        //Eftmp.gather(eg_m, this->R, IntrplCIC_t());

        /** Magnetic field in x and y direction induced by the electric field.
         *
         *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
         *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
         *  \f[ B_z = B_z^{'} = 0 \f]
         *
         */
        double betaC = sqrt(gammaz * gammaz - 1.0) / gammaz / Physics::c;

        Bf(0) = Bf(0) - betaC * Eftmp(1);
        Bf(1) = Bf(1) + betaC * Eftmp(0);

        Ef += Eftmp;


        /// Now compute field due to image charge. This is done separately as the image charge
        /// is moving to -infinity (instead of +infinity) so the Lorentz transform is different.

        /// Find z shift for shifted Green's function.
        Vector_t rmax, rmin;
        get_bounds(rmin, rmax);
        double zshift = - (rmax(2) + rmin(2)) * gammaz * scaleFactor;

        /// Find potential from image charge in this bin using Poisson's
        /// equation (without coefficient: -1/(eps)).
        IpplTimings::startTimer(compPotenTimer_m);
        fs_m->solver_m->computePotential(imagePotential, hr_scaled, zshift);
        IpplTimings::stopTimer(compPotenTimer_m);

        /// Scale mesh back (to same units as particle locations.)
        imagePotential *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// The scalar potential is given back in rho_m
        /// and must be converted to the right units.
        imagePotential *= getCouplingConstant();

        /// IPPL Grad numerical computes gradient to find the
        /// electric field (in rest frame of this bin's image
        /// charge).
        eg_m = -Grad(imagePotential, eg_m);

        /// Scale field. Combine Lorentz transform with conversion to proper field
        /// units.
        eg_m *= Vector_t(gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz));

        /// Interpolate electric field at particle positions.  We reuse the
        /// cached information about where the particles are relative to the
        /// field, since the particles have not moved since this the most recent
        /// scatter operation.
        Eftmp.gather(eg_m, IntrplCIC_t(), interpolationCache_m);
        //Eftmp.gather(eg_m, this->R, IntrplCIC_t());

        /** Magnetic field in x and y direction induced by the image charge electric field. Note that beta will have
         *  the opposite sign from the bunch charge field, as the image charge is moving in the opposite direction.
         *
         *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
         *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
         *  \f[ B_z = B_z^{'} = 0 \f]
         *
         */
        Bf(0) = Bf(0) + betaC * Eftmp(1);
        Bf(1) = Bf(1) - betaC * Eftmp(0);

        Ef += Eftmp;
    }
    IpplTimings::stopTimer(selfFieldTimer_m);
}

void PartBunch::resizeMesh() {
    double xmin = fs_m->solver_m->getXRangeMin();
    double xmax = fs_m->solver_m->getXRangeMax();
    double ymin = fs_m->solver_m->getYRangeMin();
    double ymax = fs_m->solver_m->getYRangeMax();
    double zmin = fs_m->solver_m->getZRangeMin();
    double zmax = fs_m->solver_m->getZRangeMax();
    
    if(xmin > rmin_m[0] || xmax < rmax_m[0] ||
       ymin > rmin_m[1] || ymax < rmax_m[1]) {

        for(unsigned int n = 0; n < getLocalNum(); n++) {

            if(R[n](0) < xmin || R[n](0) > xmax ||
               R[n](1) < ymin || R[n](1) > ymax) {

                // delete the particle
                INFOMSG("destroyed particle with id=" << n << endl;);
                destroy(1, n);
            }

        }

        update();
        boundp();
        get_bounds(rmin_m, rmax_m);
    }
 
    Vector_t mymin = Vector_t(xmin, ymin , zmin);
    Vector_t mymax = Vector_t(xmax, ymax , zmax);

    for(int i = 0; i < 3; i++) 
        hr_m[i]   = (mymax[i] - mymin[i]) / nr_m[i];

    getMesh().set_meshSpacing(&(hr_m[0]));
    getMesh().set_origin(mymin); 

    rho_m.initialize(getMesh(),
                     getFieldLayout(),
                     GuardCellSizes<Dim>(1),
                     bc_m);
    eg_m.initialize(getMesh(),
                    getFieldLayout(),
                    GuardCellSizes<Dim>(1),
                    vbc_m);

    update();

//    setGridIsFixed();
}

void PartBunch::computeSelfFields() {
    IpplTimings::startTimer(selfFieldTimer_m);
    rho_m = 0.0;
    eg_m = Vector_t(0.0);

    if(fs_m->hasValidSolver()) {

//        if(fs_m->getFieldSolverType() == "SAAMG" && !isGridFixed())
//            resizeMesh();
	INFOMSG("after resizeMesh" << hr_m << endl);

        //scatter charges onto grid
        this->Q *= this->dt;
        this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());
        this->Q /= this->dt;
        this->rho_m /= getdT();

        //calculating mesh-scale factor
        double gammaz = sum(this->P)[2] / getTotalNum();
        gammaz *= gammaz;
        gammaz = sqrt(gammaz + 1.0);
        double scaleFactor = 1;
        // double scaleFactor = Physics::c * getdT();
        //and get meshspacings in real units [m]
        Vector_t hr_scaled = hr_m * Vector_t(scaleFactor);
        hr_scaled[2] *= gammaz;

        //double tmp2 = 1/hr_m[0] * 1/hr_m[1] * 1/hr_m[2] / (scaleFactor*scaleFactor*scaleFactor) / gammaz;
        double tmp2 = 1 / hr_scaled[0] * 1 / hr_scaled[1] * 1 / hr_scaled[2];
        //divide charge by a 'grid-cube' volume to get [C/m^3]
        rho_m *= tmp2;

#ifdef DBG_SCALARFIELD
        INFOMSG("*** START DUMPING SCALAR FIELD ***" << endl);
        ofstream fstr1;
        fstr1.precision(9);

        std::ostringstream istr;
        istr << fieldDBGStep_m;

        string SfileName = OpalData::getInstance()->getInputBasename();

        string rho_fn = string("data/") + SfileName + string("-rho_scalar-") + string(istr.str());
        fstr1.open(rho_fn.c_str(), ios::out);
        NDIndex<3> myidx1 = getFieldLayout().getLocalNDIndex();
        for(int x = myidx1[0].first(); x <= myidx1[0].last(); x++) {
            for(int y = myidx1[1].first(); y <= myidx1[1].last(); y++) {
                for(int z = myidx1[2].first(); z <= myidx1[2].last(); z++) {
                    fstr1 << x + 1 << " " << y + 1 << " " << z + 1 << " " <<  rho_m[x][y][z].get() << endl;
                }
            }
        }
        fstr1.close();
        INFOMSG("*** FINISHED DUMPING SCALAR FIELD ***" << endl);
#endif
        // charge density is in rho_m
        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        //do the multiplication of the grid-cube volume coming
        //from the discretization of the convolution integral.
        //this is only necessary for the FFT solver
        //FIXME: later move this scaling into FFTPoissonSolver
        if(fs_m->getFieldSolverType() == "FFT" || fs_m->getFieldSolverType() == "FFTBOX")
            rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        // the scalar potential is given back in rho_m in units
        // [C/m] = [F*V/m] and must be divided by
        // 4*pi*\epsilon_0 [F/m] resulting in [V]
        rho_m *= getCouplingConstant();

        //write out rho


#ifdef DBG_SCALARFIELD
        INFOMSG("*** START DUMPING SCALAR FIELD ***" << endl);
        ostringstream oss;

        ofstream fstr2;
        fstr2.precision(9);

        string phi_fn = string("data/") + SfileName + string("-phi_scalar-") + string(istr.str());
        fstr2.open(phi_fn.c_str(), ios::out);
        NDIndex<3> myidx = getFieldLayout().getLocalNDIndex();
        for(int x = myidx[0].first(); x <= myidx[0].last(); x++) {
            for(int y = myidx[1].first(); y <= myidx[1].last(); y++) {
                for(int z = myidx[2].first(); z <= myidx[2].last(); z++) {
                    fstr2 << x + 1 << " " << y + 1 << " " << z + 1 << " " <<  rho_m[x][y][z].get() << endl;
                }
            }
        }
        fstr2.close();

        INFOMSG("*** FINISHED DUMPING SCALAR FIELD ***" << endl);
#endif

        // IPPL Grad divides by hr_m [m] resulting in
        // [V/m] for the electric field
        eg_m = -Grad(rho_m, eg_m);

        eg_m *= Vector_t(gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz));

        //write out e field
#ifdef DBG_SCALARFIELD
        INFOMSG("*** START DUMPING E FIELD ***" << endl);
        //ostringstream oss;
        //MPI_File file;
        //MPI_Status status;
        //MPI_Info fileinfo;
        //MPI_File_open(Ippl::getComm(), "rho_scalar", MPI_MODE_WRONLY | MPI_MODE_CREATE, fileinfo, &file);
        ofstream fstr;
        fstr.precision(9);

        string e_field = string("data/") + SfileName + string("-e_field-") + string(istr.str());
        fstr.open(e_field.c_str(), ios::out);
        NDIndex<3> myidxx = getFieldLayout().getLocalNDIndex();
        for(int x = myidxx[0].first(); x <= myidxx[0].last(); x++) {
            for(int y = myidxx[1].first(); y <= myidxx[1].last(); y++) {
                for(int z = myidxx[2].first(); z <= myidxx[2].last(); z++) {
                    fstr << x + 1 << " " << y + 1 << " " << z + 1 << " " <<  eg_m[x][y][z].get() << endl;
                }
            }
        }

        fstr.close();
        fieldDBGStep_m++;

        //MPI_File_write_shared(file, (char*)oss.str().c_str(), oss.str().length(), MPI_CHAR, &status);
        //MPI_File_close(&file);

        INFOMSG("*** FINISHED DUMPING E FIELD ***" << endl);
#endif

        // interpolate electric field at particle positions.  We reuse the
        // cached information about where the particles are relative to the
        // field, since the particles have not moved since this the most recent
        // scatter operation.
        Ef.gather(eg_m, this->R,  IntrplCIC_t());

        /** Magnetic field in x and y direction induced by the eletric field
         *
         *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
         *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
         *  \f[ B_z = B_z^{'} = 0 \f]
         *
         */
        double betaC = sqrt(gammaz * gammaz - 1.0) / gammaz / Physics::c;

        Bf(0) = Bf(0) - betaC * Ef(1);
        Bf(1) = Bf(1) + betaC * Ef(0);
    }
    IpplTimings::stopTimer(selfFieldTimer_m);
}

void PartBunch::computeSelfFields_cycl(double gamma) {
    IpplTimings::startTimer(selfFieldTimer_m);

    /// set initial charge density to zero.
    rho_m = 0.0;

    /// set initial E field to zero
    eg_m = Vector_t(0.0);

    if(fs_m->hasValidSolver()) {

        /// scatter particles charge onto grid.
        this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

        /// from charge to charge density.
        double tmp2 = 1.0 / gamma / (hr_m[0] * hr_m[1] * hr_m[2]);
        rho_m *= tmp2;

        /// Lorentz transformation
        /// In particle rest frame, the longitudinal length enlarged
        Vector_t hr_scaled = hr_m ;
        hr_scaled[1] *= gamma;

        /// now charge density is in rho_m
        /// calculate Possion equation (without coefficient: -1/(eps))
        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        /// additional work of FFT solver
        /// now the scalar potential is given back in rho_m
        rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// retrive coefficient: -1/(eps)
        rho_m *= getCouplingConstant();

        /// calculate electric field vectors from field potential
        eg_m = -Grad(rho_m, eg_m);

        /// back Lorentz transformation
        eg_m *= Vector_t(gamma, 1.0, gamma);

        /*
        //debug
        // output field on the grid points

        int m1 = (int)nr_m[0]-1;
        int m2 = (int)nr_m[0]/2;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along x axis E = " << eg_m[i][m2][m2] << " Pot = " << rho_m[i][m2][m2]  << endl;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along y axis E = " << eg_m[m2][i][m2] << " Pot = " << rho_m[m2][i][m2]  << endl;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along z axis E = " << eg_m[m2][m2][i] << " Pot = " << rho_m[m2][m2][i]  << endl;
        // end debug
         */

        /// interpolate electric field at particle positions.
        Ef.gather(eg_m, this->R,  IntrplCIC_t());

        /// calculate coefficient
        double betaC = sqrt(gamma * gamma - 1.0) / gamma / Physics::c;

        /// calculate B field from E field
        Bf(0) =  betaC * Ef(2);
        Bf(2) = -betaC * Ef(0);

    }
    // *gmsg<<"gamma ="<<gamma<<endl;
    // *gmsg<<"dx,dy,dz =("<<hr_m[0]<<", "<<hr_m[1]<<", "<<hr_m[2]<<") [m] "<<endl;
    // *gmsg<<"max of bunch is ("<<rmax_m(0)<<", "<<rmax_m(1)<<", "<<rmax_m(2)<<") [m] "<<endl;
    // *gmsg<<"min of bunch is ("<<rmin_m(0)<<", "<<rmin_m(1)<<", "<<rmin_m(2)<<") [m] "<<endl;
    IpplTimings::stopTimer(selfFieldTimer_m);
}

void PartBunch::computeSelfFields_cycl(double gamma, Vector_t const meanR, Vektor<double, 4> const quaternion) {
    IpplTimings::startTimer(selfFieldTimer_m);

    globalMeanR = meanR;
    globalToLocalQuaternion = quaternion;

    /// set initial charge density to zero.
    rho_m = 0.0;

    /// set initial E field to zero
    eg_m = Vector_t(0.0);

    if(fs_m->hasValidSolver()) { 
        /// scatter particles charge onto grid.
        this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

        /// from charge to charge density.
        double tmp2 = 1.0 / gamma / (hr_m[0] * hr_m[1] * hr_m[2]);
        rho_m *= tmp2;

        /// Lorentz transformation
        /// In particle rest frame, the longitudinal length enlarged
        Vector_t hr_scaled = hr_m ;
        hr_scaled[1] *= gamma;

        /// now charge density is in rho_m
        /// calculate Possion equation (without coefficient: -1/(eps))  
        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        /// additional work of FFT solver
        /// now the scalar potential is given back in rho_m
        rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// retrive coefficient: -1/(eps)
        rho_m *= getCouplingConstant();

        /// calculate electric field vectors from field potential
        eg_m = -Grad(rho_m, eg_m);

        /// back Lorentz transformation
        eg_m *= Vector_t(gamma, 1.0, gamma);

        /*
        //debug
        // output field on the grid points

        int m1 = (int)nr_m[0]-1;
        int m2 = (int)nr_m[0]/2;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along x axis E = " << eg_m[i][m2][m2] << " Pot = " << rho_m[i][m2][m2]  << endl;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along y axis E = " << eg_m[m2][i][m2] << " Pot = " << rho_m[m2][i][m2]  << endl;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along z axis E = " << eg_m[m2][m2][i] << " Pot = " << rho_m[m2][m2][i]  << endl;
        // end debug
         */
        /// interpolate electric field at particle positions.
        Ef.gather(eg_m, this->R,  IntrplCIC_t());

        /// calculate coefficient
        double betaC = sqrt(gamma * gamma - 1.0) / gamma / Physics::c;

        /// calculate B field from E field
        Bf(0) =  betaC * Ef(2);
        Bf(2) = -betaC * Ef(0);

    }
    // *gmsg<<"gamma ="<<gamma<<endl;
    // *gmsg<<"dx,dy,dz =("<<hr_m[0]<<", "<<hr_m[1]<<", "<<hr_m[2]<<") [m] "<<endl;
    // *gmsg<<"max of bunch is ("<<rmax_m(0)<<", "<<rmax_m(1)<<", "<<rmax_m(2)<<") [m] "<<endl;
    // *gmsg<<"min of bunch is ("<<rmin_m(0)<<", "<<rmin_m(1)<<", "<<rmin_m(2)<<") [m] "<<endl;
    IpplTimings::stopTimer(selfFieldTimer_m);
}

void PartBunch::computeSelfFields_cycl(int bin, Vector_t const meanR, Vektor<double, 4> const quaternion) {
    IpplTimings::startTimer(selfFieldTimer_m);

    globalMeanR = meanR;
    globalToLocalQuaternion = quaternion;

    /// set initial charge dentsity to zero.
    rho_m = 0.0;

    /// set initial E field to zero
    eg_m = Vector_t(0.0);

    /// get gamma of this bin
    double gamma = getBinGamma(bin);

    if(fs_m->hasValidSolver()) {

        /// scatter particles charge onto grid.
        this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

        /// from charge to charge density.
        double tmp2 = 1.0 / gamma / (hr_m[0] * hr_m[1] * hr_m[2]);
        rho_m *= tmp2;

        /// Lorentz transformation
        /// In particle rest frame, the longitudinal length enlarged
        Vector_t hr_scaled = hr_m ;
        hr_scaled[1] *= gamma;

        /// now charge density is in rho_m
        /// calculate Possion equation (without coefficient: -1/(eps))
        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        /// additional work of FFT solver
        /// now the scalar potential is given back in rho_m
        rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// retrive coefficient: -1/(eps)
        rho_m *= getCouplingConstant();

        /// calculate electric field vectors from field potential
        eg_m = -Grad(rho_m, eg_m);

        /// back Lorentz transformation
        eg_m *= Vector_t(gamma, 1.0, gamma);

        /*
        //debug
        // output field on the grid points

        int m1 = (int)nr_m[0]-1;
        int m2 = (int)nr_m[0]/2;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along x axis E = " << eg_m[i][m2][m2] << " Pot = " << rho_m[i][m2][m2]  << endl;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along y axis E = " << eg_m[m2][i][m2] << " Pot = " << rho_m[m2][i][m2]  << endl;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along z axis E = " << eg_m[m2][m2][i] << " Pot = " << rho_m[m2][m2][i]  << endl;
        // end debug
         */

        /// interpolate electric field at particle positions.
        Eftmp.gather(eg_m, this->R,  IntrplCIC_t());

        /// calculate coefficient
        double betaC = sqrt(gamma * gamma - 1.0) / gamma / Physics::c;

        /// calculate B_bin field from E_bin field accumulate B and E field
        Bf(0) = Bf(0) + betaC * Eftmp(2);
        Bf(2) = Bf(2) - betaC * Eftmp(0);

        Ef += Eftmp;
    }
    // *gmsg<<"gamma ="<<gamma<<endl;
    // *gmsg<<"dx,dy,dz =("<<hr_m[0]<<", "<<hr_m[1]<<", "<<hr_m[2]<<") [m] "<<endl;
    // *gmsg<<"max of bunch is ("<<rmax_m(0)<<", "<<rmax_m(1)<<", "<<rmax_m(2)<<") [m] "<<endl;
    // *gmsg<<"min of bunch is ("<<rmin_m(0)<<", "<<rmin_m(1)<<", "<<rmin_m(2)<<") [m] "<<endl;
    IpplTimings::stopTimer(selfFieldTimer_m);
}

void PartBunch::computeSelfFields_cycl(int bin) {
    IpplTimings::startTimer(selfFieldTimer_m);

    /// set initial charge dentsity to zero.
    rho_m = 0.0;

    /// set initial E field to zero
    eg_m = Vector_t(0.0);

    /// get gamma of this bin
    double gamma = getBinGamma(bin);

    if(fs_m->hasValidSolver()) {

        /// scatter particles charge onto grid.
        this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

        /// from charge to charge density.
        double tmp2 = 1.0 / gamma / (hr_m[0] * hr_m[1] * hr_m[2]);
        rho_m *= tmp2;

        /// Lorentz transformation
        /// In particle rest frame, the longitudinal length enlarged
        Vector_t hr_scaled = hr_m ;
        hr_scaled[1] *= gamma;

        /// now charge density is in rho_m
        /// calculate Possion equation (without coefficient: -1/(eps))
        fs_m->solver_m->computePotential(rho_m, hr_scaled);

        /// additional work of FFT solver
        /// now the scalar potential is given back in rho_m
        rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// retrive coefficient: -1/(eps)
        rho_m *= getCouplingConstant();

        /// calculate electric field vectors from field potential
        eg_m = -Grad(rho_m, eg_m);

        /// back Lorentz transformation
        eg_m *= Vector_t(gamma, 1.0, gamma);

        /*
        //debug
        // output field on the grid points

        int m1 = (int)nr_m[0]-1;
        int m2 = (int)nr_m[0]/2;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along x axis E = " << eg_m[i][m2][m2] << " Pot = " << rho_m[i][m2][m2]  << endl;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along y axis E = " << eg_m[m2][i][m2] << " Pot = " << rho_m[m2][i][m2]  << endl;

        for (int i=0; i<m1; i++ )
         *gmsg << "Field along z axis E = " << eg_m[m2][m2][i] << " Pot = " << rho_m[m2][m2][i]  << endl;
        // end debug
         */

        /// interpolate electric field at particle positions.
        Eftmp.gather(eg_m, this->R,  IntrplCIC_t());

        /// calculate coefficient
        double betaC = sqrt(gamma * gamma - 1.0) / gamma / Physics::c;

        /// calculate B_bin field from E_bin field accumulate B and E field
        Bf(0) = Bf(0) + betaC * Eftmp(2);
        Bf(2) = Bf(2) - betaC * Eftmp(0);

        Ef += Eftmp;
    }
    // *gmsg<<"gamma ="<<gamma<<endl;
    // *gmsg<<"dx,dy,dz =("<<hr_m[0]<<", "<<hr_m[1]<<", "<<hr_m[2]<<") [m] "<<endl;
    // *gmsg<<"max of bunch is ("<<rmax_m(0)<<", "<<rmax_m(1)<<", "<<rmax_m(2)<<") [m] "<<endl;
    // *gmsg<<"min of bunch is ("<<rmin_m(0)<<", "<<rmin_m(1)<<", "<<rmin_m(2)<<") [m] "<<endl;
    IpplTimings::stopTimer(selfFieldTimer_m);
}

void PartBunch::setBCAllOpen() {
    for(int i = 0; i < 2 * 3; ++i) {
        bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
        vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
        getBConds()[i] = ParticleNoBCond;
    }
    dcBeam_m=false;
    INFOMSG("BC set for normal Beam" << endl);
}

void PartBunch::setBCForDCBeam() {
    for(int i = 0; i < 2 * 3; ++i) {
        bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
        vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
        getBConds()[i] = ParticleNoBCond;
    }
    
    // z-direction
    bc_m[4] = new ParallelPeriodicFace<double,3,Mesh_t,Center_t>(4);
    this->getBConds()[4] = ParticlePeriodicBCond;
    bc_m[5] = new ParallelPeriodicFace<double,3,Mesh_t,Center_t>(5);
    this->getBConds()[5] = ParticlePeriodicBCond;
    dcBeam_m=true;
    INFOMSG("BC set for DC-Beam" << endl);
}

void PartBunch::boundp() {
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
	double hzSave;

	NDIndex<3> domain = getFieldLayout().getDomain();
	for(int i = 0; i < Dim; i++)
	  nr_m[i] = domain[i].length();
	get_bounds(rmin_m, rmax_m);
	Vector_t len = rmax_m - rmin_m;
	
	if (!fullUpdate) {
	  hzSave = hr_m[2];
	}
	else {
	  for(int i = 0; i < dimIdx; i++) {
	    rmax_m[i] += dh_m * abs(rmax_m[i] - rmin_m[i]);
	    rmin_m[i] -= dh_m * abs(rmax_m[i] - rmin_m[i]);
	    hr_m[i]    = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
	  }
	  //INFOMSG("It is a full boundp hz= " << hr_m << " rmax= " << rmax_m << " rmin= " << rmin_m << endl);
	}

	if (!fullUpdate) {
	  hr_m[2] = hzSave;
	  //INFOMSG("It is not a full boundp hz= " << hr_m << " rmax= " << rmax_m << " rmin= " << rmin_m << endl);
	}

	if(hr_m[0] * hr_m[1] * hr_m[2] > 0) {
	  getMesh().set_meshSpacing(&(hr_m[0]));
	  getMesh().set_origin(rmin_m - Vector_t(hr_m[0] / 2.0, hr_m[1] / 2.0, hr_m[2] / 2.0));
	  rho_m.initialize(getMesh(),
			   getFieldLayout(),
			   GuardCellSizes<Dim>(1),
			   bc_m);
	  eg_m.initialize(getMesh(),
			  getFieldLayout(),
			  GuardCellSizes<Dim>(1),
			  vbc_m);
	}	
    }
    update();
    R.resetDirtyFlag();

    IpplTimings::stopTimer(boundpTimer_m);
}

void PartBunch::calcWeightedAverages(Vector_t &CentroidPosition, Vector_t &CentroidMomentum) const {
    double gamma;
    double cent[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double N =  static_cast<double>(this->getTotalNum());

    for(unsigned int i = 0; i < this->getLocalNum(); i++) {
        gamma = sqrt(1.0 + dot(this->P[i], this->P[i]));
        cent[0] += this->R[i](0);
        cent[1] += this->R[i](1);
        cent[2] += this->R[i](2);
        cent[3] += this->P[i](0) / gamma;
        cent[4] += this->P[i](1) / gamma;
        cent[5] += this->P[i](2) / gamma;

    }
    reduce(&(cent[0]), &(cent[0]) + 6, &(cent[0]), OpAddAssign());

    CentroidPosition(0) = cent[0] / N;
    CentroidPosition(1) = cent[1] / N;
    CentroidPosition(2) = cent[2] / N;
    CentroidMomentum(0) = cent[3] / N;
    CentroidMomentum(1) = cent[4] / N;
    CentroidMomentum(2) = cent[5] / N;
}

void PartBunch::rotateAbout(const Vector_t &Center, const Vector_t &Momentum) {
    double AbsMomentumProj = sqrt(Momentum(0) * Momentum(0) + Momentum(2) * Momentum(2));
    double AbsMomentum = sqrt(dot(Momentum, Momentum));
    double cos0 = AbsMomentumProj / AbsMomentum;
    double sin0 = -Momentum(1) / AbsMomentum;
    double cos1 = Momentum(2) / AbsMomentumProj;
    double sin1 = -Momentum(0) / AbsMomentumProj;
    double sin2 = 0.0;
    double cos2 = sqrt(1.0 - sin2 * sin2);


    Tenzor<double, 3> T1(1.0,   0.0,  0.0,
                         0.0,  cos0, sin0,
                         0.0, -sin0, cos0);
    Tenzor<double, 3> T2(cos1, 0.0, sin1,
                         0.0, 1.0,  0.0,
                         -sin1, 0.0, cos1);
    Tenzor<double, 3> T3(cos2, 0.0, sin2,
                         0.0, 1.0,  0.0,
                         -sin2, 0.0, cos2);
    Tenzor<double, 3> TM = dot(T1, dot(T2, T3));
    for(unsigned int i = 0; i < this->getLocalNum(); i++) {
        R[i] = dot(TM, R[i] - Center) + Center;
        P[i] = dot(TM, P[i]);
    }
}

void PartBunch::moveBy(const Vector_t &Center) {
    for(unsigned int i = 0; i < this->getLocalNum(); i++) {
        R[i] += Center;
    }
}

void PartBunch::ResetLocalCoordinateSystem(const int &i, const Vector_t &Orientation, const double &origin) {

    Vector_t temp(R[i](0), R[i](1), R[i](2) - origin);

    if(fabs(Orientation(0)) > 1e-6 || fabs(Orientation(1)) > 1e-6 || fabs(Orientation(2)) > 1e-6) {

        // Rotate to the the element's local coordinate system.
        //
        // 1) Rotate about the z axis by angle negative Orientation(2).
        // 2) Rotate about the y axis by angle negative Orientation(0).
        // 3) Rotate about the x axis by angle Orientation(1).

        double cosa = cos(Orientation(0));
        double sina = sin(Orientation(0));
        double cosb = cos(Orientation(1));
        double sinb = sin(Orientation(1));
        double cosc = cos(Orientation(2));
        double sinc = sin(Orientation(2));

        X[i](0) = (cosa * cosc) * temp(0) + (cosa * sinc) * temp(1) - sina *        temp(2);
        X[i](1) = (-cosb * sinc - sina * sinb * cosc) * temp(0) + (cosb * cosc - sina * sinb * sinc) * temp(1) - cosa * sinb * temp(2);
        X[i](2) = (-sinb * sinc + sina * cosb * cosc) * temp(0) + (sinb * cosc + sina * cosb * sinc) * temp(1) + cosa * cosb * temp(2);

    } else
        X[i] = temp;
}


void PartBunch::beamEllipsoid(FVector<double, 6>   &centroid,
                              FMatrix<double, 6, 6> &moment) {
    for(int i = 0; i < 6; ++i) {
        centroid(i) = 0.0;
        for(int j = 0; j <= i; ++j) {
            moment(i, j) = 0.0;
        }
    }

    //  PartBunch::const_iterator last = end();
    // for (PartBunch::const_iterator part = begin(); part != last; ++part) {

    Particle part;

    for(unsigned int ii = 0; ii < this->getLocalNum(); ii++) {
        part = get_part(ii);
        for(int i = 0; i < 6; ++i) {
            centroid(i) += part[i];
            for(int j = 0; j <= i; ++j) {
                moment(i, j) += part[i] * part[j];
            }
        }
    }

    double factor = 1.0 / double(this->getTotalNum());
    for(int i = 0; i < 6; ++i) {
        centroid(i) *= factor;
        for(int j = 0; j <= i; ++j) {
            moment(j, i) = moment(i, j) *= factor;
        }
    }
}


void PartBunch::gatherLoadBalanceStatistics() {
    for(int i = 0; i < Ippl::getNodes(); i++)
        partPerNode_m[i] = globalPartPerNode_m[i] = 0.0;

    partPerNode_m[Ippl::myNode()] = this->getLocalNum();

    reduce(partPerNode_m.get(), partPerNode_m.get() + Ippl::getNodes(), globalPartPerNode_m.get(), OpAddAssign());

}

void PartBunch::calcMoments() {

    double part[2 * Dim];

    double loc_centroid[2 * Dim];
    double loc_moment[2 * Dim][2 * Dim];
    double moments[2 * Dim][2 * Dim];

    for(int i = 0; i < 2 * Dim; i++) {
        loc_centroid[i] = 0.0;
        for(int j = 0; j <= i; j++) {
            loc_moment[i][j] = 0.0;
            loc_moment[j][i] = loc_moment[i][j];
        }
    }

    for(unsigned long k = 0; k < this->getLocalNum(); ++k) {
        part[1] = this->P[k](0);
        part[3] = this->P[k](1);
        part[5] = this->P[k](2);
        part[0] = this->R[k](0);
        part[2] = this->R[k](1);
        part[4] = this->R[k](2);

        for(int i = 0; i < 2 * Dim; i++) {
            loc_centroid[i]   += part[i];
            for(int j = 0; j <= i; j++) {
                loc_moment[i][j]   += part[i] * part[j];
            }
        }
    }

    for(int i = 0; i < 2 * Dim; i++) {
        for(int j = 0; j < i; j++) {
            loc_moment[j][i] = loc_moment[i][j];
        }
    }

    reduce(&(loc_moment[0][0]), &(loc_moment[0][0]) + 2 * Dim * 2 * Dim,
           &(moments[0][0]), OpAddAssign());

    reduce(&(loc_centroid[0]), &(loc_centroid[0]) + 2 * Dim,
           &(centroid_m[0]), OpAddAssign());

    for(int i = 0; i < 2 * Dim; i++) {
        for(int j = 0; j <= i; j++) {
            moments_m(i, j) = moments[i][j];
            moments_m(j, i) = moments_m(i, j);
        }
    }
}

void PartBunch::calcMomentsInitial() {

    double part[2 * Dim];

    for(int i = 0; i < 2 * Dim; ++i) {
        centroid_m[i] = 0.0;
        for(int j = 0; j <= i; ++j) {
            moments_m(i, j) = 0.0;
            moments_m(j, i) = moments_m(i, j);
        }
    }

    for(size_t k = 0; k < pbin_m->getNp(); k++) {
        for(int binNumber = 0; binNumber < pbin_m->getNBins(); binNumber++) {
            vector<double> p;

            if(pbin_m->getPart(k, binNumber, p)) {
                part[0] = p.at(0);
                part[1] = p.at(3);
                part[2] = p.at(1);
                part[3] = p.at(4);
                part[4] = p.at(2);
                part[5] = p.at(5);

                for(int i = 0; i < 2 * Dim; ++i) {
                    centroid_m[i] += part[i];
                    for(int j = 0; j <= i; ++j) {
                        moments_m(i, j) += part[i] * part[j];
                    }
                }
            }
        }
    }

    for(int i = 0; i < 2 * Dim; ++i) {
        for(int j = 0; j < i; ++j) {
            moments_m(j, i) = moments_m(i, j);
        }
    }
}

void PartBunch::calcBeamParameters() {
    using Physics::c;

    Vector_t eps2, fac, rsqsum, psqsum, rpsum;

    const double m0 = getM() * 1.E-6;

    IpplTimings::startTimer(statParamTimer_m);

    const size_t locNp = this->getLocalNum();
    const double N =  static_cast<double>(this->getTotalNum());
    const double zero = 0.0;

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

    /*
      double rmax = sqrt(dot(rmax_m,rmax_m));
      fplasma_m = sqrt(2.0*get_perverance())*get_beta()*c/rmax;
      budkerp_m = (get_perverance()/2.0)*pow(get_gamma(),3.0)*pow(get_beta(),2.0);
     */

    // Find unnormalized emittance.
    double gamma = 0.0;
    for(size_t i = 0; i < locNp; i++)
        gamma += sqrt(1.0 + dot(P[i], P[i]));

    reduce(gamma, gamma, OpAddAssign());
    gamma /= N;
    eKin_m = (gamma - 1.0) * m0;
    // calculate energy spread
    dE_m = prms_m(2) * sqrt(eKin_m * (eKin_m + 2.*m0) / (1. + eKin_m * (eKin_m + 2.*m0) / (m0 * m0)));

    eps_m = eps_norm_m / Vector_t(gamma * sqrt(1.0 - 1.0 / (gamma * gamma)));
    IpplTimings::stopTimer(statParamTimer_m);
}

void PartBunch::calcBeamParametersLight() {
    // for Autophase, avoids communication

    const double m0 = getM() * 1.E-6;

    const size_t locNp = this->getLocalNum();

    // Find unnormalized emittance.
    double gamma = 0.0;
    for(size_t i = 0; i < locNp; i++)
        gamma += sqrt(1.0 + dot(P[i], P[i]));

    eKin_m = (gamma - 1.0) * m0;

    IpplTimings::stopTimer(statParamTimer_m);
}

void PartBunch::calcBeamParametersInitial() {
    using Physics::c;

    const double N =  static_cast<double>(this->getTotalNum());

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
                eps_m(i)  = sqrt(max(eps2(i), zero));
                double tmp = rrms_m(i) * prms_m(i);
                fac(i) = (tmp == 0) ? zero : 1.0 / tmp;
            }
            rprms_m = rpsum * fac;
        }
    }
}

void PartBunch::setSolver(FieldSolver *fs) {
    fs_m = fs;
    fs_m->initSolver(*this);
    /**
       CAN not re-inizialize ParticleLayout
       this is an IPPL issue
     */
    if(!OpalData::getInstance()->hasBunchAllocated())
        initialize(&(fs_m->getParticleLayout()));
}

void PartBunch::maximumAmplitudes(const FMatrix<double, 6, 6> &D,
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

/**
 * Here we emit particles from the cathode. All particles in a new simulation (not a restart) initially reside in the bin
 container "pbin_m" and are not part of the beam bunch (so they cannot "see" fields, space charge etc.). In pbin_m, particles
 are sorted into the bins of a time histogram that describes the longitudinal time distribution of the beam, where the number
 of bins is given by \f$NBIN \times SBIN\f$. \f$NBIN\f$ and \f$SBIN\f$ are parameters given when defining the initial beam
 distribution. During emission, the time step of the simulation is set so that an integral number of these bins are emitted each step.
 Once all of the particles have been emitted, the simulation time step is reset to the value defined in the input file.

 A typical integration time step, \f$\Delta t\f$, is broken down into 3 sub-steps:

 1) Drift particles for \f$\frac{\Delta t}{2}\f$.

 2) Calculate fields and advance momentum.

 3) Drift particles for \f$\frac{\Delta t}{2}\f$ at the new momentum to complete the
 full time step.

 The difficulty for emission is that at the cathode position there is a step function discontinuity in the  fields. If we
 apply the typical integration time step across this boundary, we get an artificial numerical bunching of the beam, especially
 at very high accelerating fields. This function takes the cathode position boundary into account in order to achieve
 smoother particle emission.

 During an emission step, an integral number of time bins from the distribution histogram are emitted. However, each particle
 contained in those time bins will actually be emitted from the cathode at a different time, so will only spend some fraction
 of the total time step, \f$\Delta t_{full-timestep}\f$, in the simulation. The trick to emission is to give each particle
 a unique time step, \f$Delta t_{temp}\f$, that is equal to the actual time during the emission step that the particle
 exists in the simulation. For the next integration time step, the particle's time step is set back to the global time step,
 \f$\Delta t_{full-timestep}\f$.
  */

double PartBunch::getTBin() {
    if(dist_m)
        return dist_m->GetEnergyBinDeltaT();
    else
        return 0.0;
}

double PartBunch::GetEmissionDeltaT() {
    return dist_m->GetEmissionDeltaT();
}

size_t PartBunch::EmitParticles(double eZ) {

    return dist_m->EmitParticles(*this, eZ);

}


double PartBunch::calcTimeDelay(const double &jifactor) {
    double gamma = pbin_m->getGamma();
    double beta = sqrt(1. - (1. / (gamma * gamma)));
    double xmin, xmax;
    pbin_m->getExtrema(xmin, xmax);
    double length = xmax - xmin;

    return length * jifactor / (Physics::c * beta);
}

void PartBunch::moveBunchToCathode(double &t) {
    double avrg_betagamma = 0.0;
    double maxspos = -9999999.99;

    for(int bin = 0; bin < getNumBins(); ++bin) {
        for(size_t i = 0; i < pbin_m->getNp(); ++i) {
            vector<double> p;
            if(pbin_m->getPart(i, bin, p)) {
                avrg_betagamma += sqrt(1.0 + p[3] * p[3] + p[4] * p[4] + p[5] * p[5]);
                if(p[2] > maxspos) maxspos = p[2];
            }
        }
    }
    avrg_betagamma /= pbin_m->getNp();
    double dist_per_step = sqrt(1.0 - (1.0 / (avrg_betagamma * avrg_betagamma))) * Physics::c * getdT();
    if(maxspos < 0.0) {
        double num_steps = floor(-maxspos / dist_per_step);
        t += num_steps * getdT();

        Inform gmsg("PartBunch");
        gmsg << "move bunch by " << num_steps *dist_per_step << "; DT = " << num_steps *getdT() << endl;
        for(int bin = 0; bin < getNumBins(); ++bin) {
            for(size_t i = 0; i < pbin_m->getNp(); ++i) {
                vector<double> p;
                if(pbin_m->getPart(i, bin, p)) {
                    pbin_m->updatePartPos(i, bin, p[2] + num_steps * dist_per_step);
                }
            }
        }
    }
}

void PartBunch::printBinHist() {
    if(weHaveBins()) {
        std::unique_ptr<int[]> binhisto(new int[getNumBins()]);
        double maxz = -999999999.99, minz = 999999999.99;

        for(int bin = 0; bin < getNumBins(); ++bin) {
            binhisto[bin] = 0;
        }
        for(size_t i = 0; i < getLocalNum(); ++i) {
            ++binhisto[this->Bin[i]];
        }
        reduce(&(binhisto[0]), &(binhisto[0]) + getNumBins(), &(binhisto[0]), OpAddAssign());

        if(Ippl::myNode() == 0) {
            for(int bin = 0; bin < getNumBins(); ++bin) {
                for(size_t i = 0; i < pbin_m->getNp(); ++i) {
                    vector<double> p;
                    if(pbin_m->getPart(i, bin, p)) {
                        if(p[2] > maxz) maxz = p[2];
                        if(p[2] < minz) minz = p[2];
                    }
                }
            }
            double dz = (maxz - minz) / 20;
            int minihist[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            for(int bin = 0; bin < getNumBins(); ++bin) {
                for(size_t i = 0; i < pbin_m->getNp(); ++i) {
                    vector<double> p;
                    if(pbin_m->getPart(i, bin, p)) {
                        ++minihist[(int)floor((p[2] - minz) / dz)];
                    }
                }
            }
            INFOMSG("particles are between " << minz << " and " << maxz << endl;);
            for(int i = 0; i < 20; ++i) {
                INFOMSG(i << ": " << minihist[i] << endl;);
            }

            INFOMSG("\n\n"
                    << "effective histogram:" << endl;);
            for(int bin = 0; bin < getNumBins(); ++bin)
                INFOMSG(bin << ": " << binhisto[bin] << endl;);
        }
    }
}


Inform &PartBunch::print(Inform &os) {
    if(this->getTotalNum() != 0) {  // to suppress Nan's
        Inform::FmtFlags_t ff = os.flags();
        os << scientific;
        os << "* ************** B U N C H ********************************************************* " << endl;
        os << "* NP              =   " << this->getTotalNum() << "\n";
        os << "* Qtot            =   " << setw(12) << setprecision(5) << abs(sum(Q)) * 1.0e9 << " [nC]       "
           << "Qi    = " << setw(12) << std::abs(qi_m) * 1e9 << " [nC]" << "\n";
        os << "* Ekin            =   " << setw(12) << setprecision(5) << eKin_m << " [MeV]      "
           << "dEkin = " << setw(12) << dE_m << " [MeV]" << endl;
        os << "* rmax            = " << setw(12) << setprecision(5) << rmax_m << " [m]" << endl;
        os << "* rmin            = " << setw(12) << setprecision(5) << rmin_m << " [m]" << endl;
        os << "* rms beam size   = " << setw(12) << setprecision(5) << rrms_m << " [m]" << endl;
        os << "* rms momenta     = " << setw(12) << setprecision(5) << prms_m << " [beta gamma]" << endl;
        os << "* mean position   = " << setw(12) << setprecision(5) << rmean_m << " [m]" << endl;
        os << "* mean momenta    = " << setw(12) << setprecision(5) << pmean_m << " [beta gamma]" << endl;
        os << "* rms emittance   = " << setw(12) << setprecision(5) << eps_m << " (not normalized)" << endl;
        os << "* rms correlation = " << setw(12) << setprecision(5) << rprms_m << endl;
        os << "* hr              = " << setw(12) << setprecision(5) << hr_m << " [m]" << endl;
        os << "* dh              =   " << setw(12) << setprecision(5) << dh_m << " [m]" << endl;
        os << "* t               =   " << setw(12) << setprecision(5) << getT() << " [s]        "
           << "dT    = " << setw(12) << getdT() << " [s]" << endl;
        os << "* spos            =   " << setw(12) << setprecision(5) << get_sPos() << " [m]" << endl;
        os << "* ********************************************************************************** " << endl;
        os.flags(ff);
    }
    return os;
}


void PartBunch::calcBeamParameters_cycl() {
    using Physics::c;

    Vector_t eps2, fac, rsqsum, psqsum, rpsum;

    const size_t locNp = this->getLocalNum();
    double localeKin = 0.0;

    const double zero = 0.0;
    const double TotalNp =  static_cast<double>(this->getTotalNum());

    // calculate centroid_m and moments_m
    calcMoments();

    for(unsigned int i = 0 ; i < Dim; i++) {
        rmean_m(i) = centroid_m[2 * i] / TotalNp;
        pmean_m(i) = centroid_m[(2 * i) + 1] / TotalNp;
        rsqsum(i) = moments_m(2 * i, 2 * i) - TotalNp * rmean_m(i) * rmean_m(i);
        psqsum(i) = moments_m((2 * i) + 1, (2 * i) + 1) - TotalNp * pmean_m(i) * pmean_m(i);
        rpsum(i) =  moments_m((2 * i), (2 * i) + 1) - TotalNp * rmean_m(i) * pmean_m(i);
    }
    eps2      = (rsqsum * psqsum - rpsum * rpsum) / (TotalNp * TotalNp);
    rpsum /= TotalNp;

    for(unsigned int i = 0 ; i < Dim; i++) {
        rrms_m(i) = sqrt(rsqsum(i) / TotalNp);
        prms_m(i) = sqrt(psqsum(i) / TotalNp);
        //eps_m(i)  = sqrt( max( eps2(i), zero ) );
        eps_norm_m(i)  = sqrt(max(eps2(i), zero));
        double tmp    = rrms_m(i) * prms_m(i);
        fac(i)  = (tmp == 0) ? zero : 1.0 / tmp;
    }

    rprms_m = rpsum * fac;

    // y: longitudinal direction; z: vertical direction.
    Dx_m = moments_m(0, 3) / TotalNp;
    DDx_m = moments_m(1, 3) / TotalNp;

    Dy_m = moments_m(4, 3) / TotalNp;
    DDy_m = moments_m(5, 3) / TotalNp;

    // calculate mean energy
    eKin_m = 0.0;
    for(unsigned int k = 0; k < locNp; k++)
        eKin_m += (sqrt(dot(P[k], P[k]) + 1.0) - 1.0) * getM() * 1e-6;

    localeKin = eKin_m / locNp;
    // sum energy of all nodes
    reduce(eKin_m, eKin_m, OpAddAssign());
    eKin_m /= TotalNp;
    
    double meanLocalBetaGamma = sqrt(pow(1 + localeKin / (getM() * 1.0e-6), 2.0) - 1);

    double betagamma = meanLocalBetaGamma * locNp;
    // sum the betagamma of all nodes
    reduce(betagamma, betagamma, OpAddAssign());
    betagamma /= TotalNp;

    // obtain the global RMS emmitance, it make no sense for multi-bunch simulation
    eps_m = eps_norm_m / Vector_t(betagamma);
}



size_t PartBunch::boundp_destroyT() {

    NDIndex<Dim> domain = getFieldLayout().getDomain();
    for(int i = 0; i < Dim; i++)
        nr_m[i] = domain[i].length();

    std::unique_ptr<size_t[]> tmpbinemitted;

    size_t ne = 0;
    boundp();

    if(WeHaveEnergyBins()) {
        tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[GetNumberOfEnergyBins()]);
        for(int i = 0; i < GetNumberOfEnergyBins(); i++) {
            tmpbinemitted[i] = 0;
        }
        for(unsigned int i = 0; i < this->getLocalNum(); i++) {
            if(Bin[i] < 0) {
                ne++;
                this->destroy(1, i);
            } else
                tmpbinemitted[Bin[i]]++;
        }
    } else {
        for(unsigned int i = 0; i < this->getLocalNum(); i++) {
            if(Bin[i] < 0) {
                ne++;
                this->destroy(1, i);
            }
        }
    }
    update();
    boundp();
    calcBeamParameters();

    if(WeHaveEnergyBins()) {
        const int lastBin = dist_m->GetLastEmittedEnergyBin() + 1;
        for(int i = 0; i < lastBin; i++) {
            binemitted_m[i] = tmpbinemitted[i];
        }
    }
    reduce(ne, ne, OpAddAssign());
    return ne;
}


void PartBunch::boundp_destroy() {

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

    NDIndex<3> domain = getFieldLayout().getDomain();
    for(int i = 0; i < Dim; i++)
        nr_m[i] = domain[i].length();

    get_bounds(rmin_m, rmax_m);
    len = rmax_m - rmin_m;

    calcBeamParameters_cycl();

    const int checkfactor = Options::remotePartDel;
    // check the bunch if its full size is larger than checkfactor times of its rms size
    if(checkfactor != -1) {
        if(len[0] > checkfactor * rrms_m[0] || len[1] > checkfactor * rrms_m[1] || len[2] > checkfactor * rrms_m[2]) {
            for(unsigned int ii = 0; ii < this->getLocalNum(); ii++) {
                // delete the particle if the ditance to the beam center is larger than 8 times of beam's rms size
                if(abs(R[ii](0) - rmean_m(0)) > checkfactor * rrms_m[0] || abs(R[ii](1) - rmean_m(1)) > checkfactor * rrms_m[1] || abs(R[ii](2) - rmean_m(2)) > checkfactor * rrms_m[2]) {
                    // put particle onto deletion list
                    destroy(1, ii);
                    //update bin parameter
                    if(weHaveBins()) countLost[Bin[ii]] += 1 ;

                    gmsgAll << "REMOTE PARTICLE DELETION: ID = " << ID[ii] << ", R = " << R[ii] << ", beam rms = " << rrms_m << endl;
                }
            }
        }
    }
    for(int i = 0; i < dimIdx; i++) {
        rmax_m[i] += dh_m * abs(rmax_m[i] - rmin_m[i]);
        rmin_m[i] -= dh_m * abs(rmax_m[i] - rmin_m[i]);
        hr_m[i]    = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
    }

    // rescale mesh
    getMesh().set_meshSpacing(&(hr_m[0]));
    getMesh().set_origin(rmin_m);

    rho_m.initialize(getMesh(),
                     getFieldLayout(),
                     GuardCellSizes<Dim>(1),
                     bc_m);
    eg_m.initialize(getMesh(),
                    getFieldLayout(),
                    GuardCellSizes<Dim>(1),
                    vbc_m);


    if(weHaveBins()) {
        pbin_m->updatePartInBin_cyc(countLost.get());
    }

    update();

    IpplTimings::stopTimer(boundpTimer_m);

}

// angle range [0~2PI) degree
double PartBunch::calculateAngle(double x, double y) {
    double thetaXY = atan2(y, x);

    // if(x < 0)                   thetaXY = pi + atan(y / x);
    // else if((x > 0) && (y >= 0))  thetaXY = atan(y / x);
    // else if((x > 0) && (y < 0))   thetaXY = 2.0 * pi + atan(y / x);
    // else if((x == 0) && (y > 0)) thetaXY = pi / 2.0;
    // else if((x == 0) && (y < 0)) thetaXY = 3.0 / 2.0 * pi;

    return thetaXY >= 0 ? thetaXY : thetaXY + Physics::two_pi;

}

// angle range [-PI~PI) degree
double PartBunch::calculateAngle2(double x, double y) {

    // double thetaXY = atan2(y, x);
    // if(x > 0)              thetaXY = atan(y / x);
    // else if((x < 0)  && (y > 0)) thetaXY = pi + atan(y / x);
    // else if((x < 0)  && (y <= 0)) thetaXY = -pi + atan(y / x);
    // else if((x == 0) && (y > 0)) thetaXY = pi / 2.0;
    // else if((x == 0) && (y < 0)) thetaXY = -pi / 2.0;

    return atan2(y, x);

}

double PartBunch::calcMeanPhi() {

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

    for(int ii = 0; ii < emittedBins; ii++) {
        reduce(px[ii], px[ii], OpAddAssign());
        reduce(py[ii], py[ii], OpAddAssign());

        phi[ii] = calculateAngle(px[ii], py[ii]);
        meanPhi += phi[ii];
        INFOMSG("Bin " << ii  << "  mean phi = " << phi[ii] * 180.0 / pi - 90.0 << "[degree]" << endl);
    }

    meanPhi /= emittedBins;

    INFOMSG("mean phi of all particles " <<  meanPhi * 180.0 / pi - 90.0 << "[degree]" << endl);

    return meanPhi;
}

size_t PartBunch::getNumPartInBin(int BinID) const {
    if(weHaveBins())
        return pbin_m->getGlobalBinCount(BinID);
    else
        return this->getTotalNum();
}

// this function reset the BinID for each particles according to its current gamma
// for the time it is designed for cyclotron where energy gain per turn may be changing.
bool PartBunch::resetPartBinID() {

    size_t partInBin[numBunch_m];

    if(numBunch_m != pbin_m->getLastemittedBin()) {
        ERRORMSG("Bunch number does NOT equal to bin number! Not implemented yet!!!" << endl);
        return false;
    }

    for(int ii = 0; ii < numBunch_m; ii++) partInBin[ii] = 0 ;

    // update mass gamma for each bin first
    INFOMSG("Before reset Bin: " << endl);
    calcGammas_cycl();

    // reset bin index for each particle and
    // calculate total particles number for each bin
    for(unsigned long int n = 0; n < this->getLocalNum(); n++) {
        double deltgamma[numBunch_m];
        double gamma = sqrt(1.0 + dot(P[n], P[n]));

        int index = 0;
        for(int ii = 0; ii < numBunch_m; ii++)
            deltgamma[ii] = abs(bingamma_m[ii] - gamma);

        for(int ii = 0; ii < numBunch_m; ii++)
            if(*(deltgamma + index) > *(deltgamma + ii))
                index = ii;

        Bin[n] = index;
        partInBin[index]++;
    }

    pbin_m->resetPartInBin(partInBin);

    // after reset Particle Bin ID, update mass gamma of each bin again
    INFOMSG("After reset Bin: " << endl);
    calcGammas_cycl();

    return true;

}

// this function reset the BinID for each particles according to its current beta*gamma
// it is for multi-turn extraction cyclotron with small energy gain
// the bin number can be different with the bunch number

bool PartBunch::resetPartBinID2(const double eta) {


    INFOMSG("Before reset Bin: " << endl);
    calcGammas_cycl();
    int maxbin = pbin_m->getNBins();
    size_t partInBin[maxbin];
    for(int ii = 0; ii < maxbin; ii++) partInBin[ii] = 0;

    double pMin0 = 1.0e9;
    double pMin = 0.0;
    double maxbinIndex = 0;

    for(unsigned long int n = 0; n < this->getLocalNum(); n++) {
        double temp_betagamma = sqrt(pow(P[n](0), 2) + pow(P[n](1), 2));
        if(pMin0 > temp_betagamma)
            pMin0 = temp_betagamma;
    }
    reduce(pMin0, pMin, OpMinAssign());
    INFOMSG("minimal beta*gamma = " << pMin << endl);

    double asinh0 = asinh(pMin);
    for(unsigned long int n = 0; n < this->getLocalNum(); n++) {

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

void PartBunch::setPBins(PartBins *pbin) {
    pbin_m = pbin;
    *gmsg << *pbin_m << endl;
    bingamma_m = std::unique_ptr<double[]>(new double[pbin_m->getNBins()]);
    binemitted_m = std::unique_ptr<size_t[]>(new size_t[pbin_m->getNBins()]);
    for(int i = 0; i < pbin_m->getNBins(); i++)
        binemitted_m[i] = 0;
}


void PartBunch::setPBins(PartBinsCyc *pbin) {

    pbin_m = pbin;
    bingamma_m = std::unique_ptr<double[]>(new double[pbin_m->getNBins()]);
    binemitted_m = std::unique_ptr<size_t[]>(new size_t[pbin_m->getNBins()]);
    for(int i = 0; i < pbin_m->getNBins(); i++)
        binemitted_m[i] = 0;

}

void PartBunch::stash() {

    size_t Nloc = getLocalNum();
    if(Nloc == 0) return;
    if(bunchStashed_m) {
        *gmsg << "ERROR: bunch already stashed, call pop() first" << endl;
        return;
    }

    // save all particles
    stash_Nloc_m = Nloc;
    stash_iniR_m = get_rmean();
    stash_iniP_m = get_pmean();

    stash_id_m.create(Nloc);
    stash_r_m.create(Nloc);
    stash_p_m.create(Nloc);
    stash_x_m.create(Nloc);
    stash_q_m.create(Nloc);
    stash_bin_m.create(Nloc);
    stash_dt_m.create(Nloc);
    stash_ls_m.create(Nloc);
    stash_ptype_m.create(Nloc);

    stash_id_m    = this->ID;
    stash_r_m     = this->R;
    stash_p_m     = this->P;
    stash_x_m     = this->X;
    stash_q_m     = this->Q;
    stash_bin_m   = this->Bin;
    stash_dt_m    = this->dt;
    stash_ls_m    = this->LastSection;
    stash_ptype_m = this->PType;

    // and destroy all particles in bunch
    destroy(Nloc, 0);
    update();

    bunchStashed_m = true;
}

void PartBunch::pop() {

    if(!bunchStashed_m) return;

    size_t Nloc = getLocalNum();
    if(getTotalNum() > 0) {
        destroy(Nloc, 0);
    }
    update();

    this->create(stash_Nloc_m);

    this->ID          = stash_id_m;
    this->R           = stash_r_m;
    this->P           = stash_p_m;
    this->X           = stash_x_m;
    this->Q           = stash_q_m;
    this->Bin         = stash_bin_m;
    this->dt          = stash_dt_m;
    this->LastSection = stash_ls_m;
    this->PType       = stash_ptype_m;

    stash_iniR_m = Vector_t(0.0);
    stash_iniP_m = Vector_t (0.0, 0.0, 1E-6);

    stash_id_m.destroy(stash_Nloc_m, 0);
    stash_r_m.destroy(stash_Nloc_m, 0);
    stash_p_m.destroy(stash_Nloc_m, 0);
    stash_x_m.destroy(stash_Nloc_m, 0);
    stash_q_m.destroy(stash_Nloc_m, 0);
    stash_dt_m.destroy(stash_Nloc_m, 0);
    stash_bin_m.destroy(stash_Nloc_m, 0);
    stash_ls_m.destroy(stash_Nloc_m, 0);
    stash_ptype_m.destroy(stash_Nloc_m, 0);

    bunchStashed_m = false;

    update();
}

double PartBunch::getZPos() {

    if(sum(PType != 0)) {
        size_t numberOfPrimaryParticles = 0;
        double zAverage = 0.0;
        if(getLocalNum() != 0) {
            for(size_t partIndex = 0; partIndex < getLocalNum(); partIndex++) {
                if(PType[partIndex] == 0) {
                    zAverage += X[partIndex](2);
                    numberOfPrimaryParticles++;
                }
            }
            if(numberOfPrimaryParticles != 0)
                zAverage /= numberOfPrimaryParticles;
        }
        reduce(zAverage, zAverage, OpAddAssign());
        zAverage /= Ippl::getNodes();

        return zAverage;
    } else {
        if(getTotalNum() > 0)
            return sum(X(2)) / getTotalNum();
        else
            return 0.0;
    }
}

void PartBunch::getXBounds(Vector_t &xMin, Vector_t &xMax) {
    bounds(X, xMin, xMax);
}

void PartBunch::iterateEmittedBin(int binNumber) {
    binemitted_m[binNumber]++;
}

bool PartBunch::doEmission() {
    if (dist_m != NULL)
        return dist_m->GetIfDistEmitting();
    else
        return false;
}

bool PartBunch::GetIfBeamEmitting() {

    if (dist_m != NULL) {
        size_t isBeamEmitted = dist_m->GetIfDistEmitting();
        reduce(isBeamEmitted, isBeamEmitted, OpAddAssign());
        if (isBeamEmitted > 0)
            return true;
        else
            return false;
    } else
        return false;

}

int PartBunch::GetLastEmittedEnergyBin() {
    /*
     * Get maximum last emitted energy bin.
     */
    int lastEmittedBin = dist_m->GetLastEmittedEnergyBin();
    std::vector<int> emittedBins;
    for (int nodeIndex = 0; nodeIndex < Ippl::getNodes(); nodeIndex++) {
        if (nodeIndex == Ippl::myNode())
            emittedBins.push_back(lastEmittedBin);
        else
            emittedBins.push_back(0);

        reduce(emittedBins.at(nodeIndex), emittedBins.at(nodeIndex), OpAddAssign());
    }

    int maxEmittedBin = 0;
    std::vector<int>::iterator binIt;
    for (binIt = emittedBins.begin(); binIt != emittedBins.end(); binIt++) {
        if (binIt == emittedBins.begin())
            maxEmittedBin = *binIt;
        else if (maxEmittedBin < *binIt)
            maxEmittedBin = *binIt;
    }

    return maxEmittedBin;
}

size_t PartBunch::GetNumberOfEmissionSteps() {
    return dist_m->GetNumberOfEmissionSteps();
}

bool PartBunch::weHaveBins() const {
    if(pbin_m != NULL)
        return pbin_m->weHaveBins();
    else
        return false;
}

int PartBunch::GetNumberOfEnergyBins() {
    return dist_m->GetNumberOfEnergyBins();
}

void PartBunch::Rebin() {

    size_t isBeamEmitting = dist_m->GetIfDistEmitting();
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

void PartBunch::SetEnergyBins(int numberOfEnergyBins) {
    bingamma_m = std::unique_ptr<double[]>(new double[numberOfEnergyBins]);
    binemitted_m = std::unique_ptr<size_t[]>(new size_t[numberOfEnergyBins]);
    for(int i = 0; i < numberOfEnergyBins; i++)
        binemitted_m[i] = 0;
}

bool PartBunch::WeHaveEnergyBins() {

    if (dist_m != NULL)
        return dist_m->GetNumberOfEnergyBins() > 0;
    else
        return false;
}
