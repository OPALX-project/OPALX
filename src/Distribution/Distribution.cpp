// ------------------------------------------------------------------------
// $RCSfile: Distribution.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Distribution
//   The class for the OPAL Distribution command.
//
// ------------------------------------------------------------------------

#include "Distribution/Distribution.h"
#include "Distribution/SigmaGenerator.h"
#include "AbsBeamline/SpecificElementVisitor.h"

#include <cmath>
#include <cfloat>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "AbstractObjects/Expressions.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include "Utilities/OpalException.h"
#include "halton1d_sequence.hh"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include "Structure/Beam.h"
#include "Structure/BoundaryGeometry.h"
#include "Algorithms/PartBinsCyc.h"
#include "BasicActions/Option.h"
#include "Distribution/LaserProfile.h"
#include "Elements/OpalBeamline.h"
#include "AbstractObjects/BeamSequence.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPC.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>


extern Inform *gmsg;

#define DISTDBG1
#define noDISTDBG2

#define SMALLESTCUTOFF 1e-12

SymTenzor<double, 6> getUnit6x6() {
    SymTenzor<double, 6> unit6x6;
    for (unsigned int i = 0; i < 6u; ++ i) {
        unit6x6(i,i) = 1.0;
    }
    return unit6x6;
}

//
// Class Distribution
// ------------------------------------------------------------------------

namespace AttributesT
{
    enum AttributesT {
        DISTRIBUTION,
        FNAME,
        WRITETOFILE,
        WEIGHT,
        INPUTMOUNITS,
        EMITTED,
        EMISSIONSTEPS,
        EMISSIONMODEL,
        EKIN,
        ELASER,
        W,
        FE,
        CATHTEMP,
        NBIN,
        XMULT,
        YMULT,
        ZMULT,
        TMULT,
        PXMULT,
        PYMULT,
        PZMULT,
        OFFSETX,
        OFFSETY,
        OFFSETZ,
        OFFSETT,
        OFFSETPX,
        OFFSETPY,
        OFFSETPZ,
        SIGMAX,
        SIGMAY,
        SIGMAR,
        SIGMAZ,
        SIGMAT,
        TPULSEFWHM,
        TRISE,
        TFALL,
        SIGMAPX,
        SIGMAPY,
        SIGMAPZ,
        MX,
        MY,
        MZ,
        MT,
        CUTOFFX,
        CUTOFFY,
        CUTOFFR,
        CUTOFFLONG,
        CUTOFFPX,
        CUTOFFPY,
        CUTOFFPZ,
        FTOSCAMPLITUDE,
        FTOSCPERIODS,
        R,                          // the correlation matrix (a la transport)
        CORRX,
        CORRY,
        CORRZ,
        CORRT,
        R51,
        R52,
        R61,
        R62,
        LASERPROFFN,
        IMAGENAME,
        INTENSITYCUT,
        FLIPX,
        FLIPY,
        ROTATE90,
        ROTATE180,
        ROTATE270,
        NPDARKCUR,
        INWARDMARGIN,
        EINITHR,
        FNA,
        FNB,
        FNY,
        FNVYZERO,
        FNVYSECOND,
        FNPHIW,
        FNBETA,
        FNFIELDTHR,
        FNMAXEMI,
        SECONDARYFLAG,
        NEMISSIONMODE,
        VSEYZERO,                   // sey_0 in Vaughn's model.
        VEZERO,                     // Energy related to sey_0 in Vaughan's model.
        VSEYMAX,                    // sey max in Vaughan's model.
        VEMAX,                      // Emax in Vaughan's model.
        VKENERGY,                   // Fitting parameter denotes the roughness of
        // surface for impact energy in Vaughn's model.
        VKTHETA,                    // Fitting parameter denotes the roughness of
        // surface for impact angle in Vaughn's model.
        VVTHERMAL,                  // Thermal velocity of Maxwellian distribution
        // of secondaries in Vaughan's model.
        VW,
        SURFMATERIAL,               // Add material type, currently 0 for copper
        // and 1 for stainless steel.
        EX,                         // below is for the matched distribution
        EY,
        ET,
        MAGSYM,                     // number of sector magnets
        LINE,
        FMAPFN,
        FMTYPE,                     // field map type used in matched gauss distribution
        RESIDUUM,
        MAXSTEPSCO,
        MAXSTEPSSI,
        ORDERMAPS,
        E2,
	RGUESS,
        ID1,                       // special particle that the user can set
        ID2,                       // special particle that the user can set
        SIZE
    };
}

namespace LegacyAttributesT
{
    enum LegacyAttributesT {
        // DESCRIPTION OF THE DISTRIBUTION:
        DEBIN = AttributesT::SIZE,
        SBIN,
        TEMISSION,
        SIGLASER,
        AG,
        SIGMAPT,
        TRANSVCUTOFF,
        CUTOFF,
        Z,
        T,
        PT,
        ALPHAX,
        ALPHAY,
        BETAX,
        BETAY,
        DX,
        DDX,
        DY,
        DDY,
        SIZE
    };
}

Distribution::Distribution():
    Definition( LegacyAttributesT::SIZE, "DISTRIBUTION",
                "The DISTRIBUTION statement defines data for the 6D particle distribution."),
    distrTypeT_m(DistrTypeT::NODIST),
    numberOfDistributions_m(1),
    emitting_m(false),
    scan_m(false),
    emissionModel_m(EmissionModelT::NONE),
    tEmission_m(0.0),
    tBin_m(0.0),
    currentEmissionTime_m(0.0),
    currentEnergyBin_m(1),
    currentSampleBin_m(0),
    numberOfEnergyBins_m(0),
    numberOfSampleBins_m(0),
    energyBins_m(NULL),
    energyBinHist_m(NULL),
    randGenEmit_m(NULL),
    pTotThermal_m(0.0),
    pmean_m(0.0),
    cathodeWorkFunc_m(0.0),
    laserEnergy_m(0.0),
    cathodeFermiEnergy_m(0.0),
    cathodeTemp_m(0.0),
    emitEnergyUpperLimit_m(0.0),
    inputMoUnits_m(InputMomentumUnitsT::NONE),
    sigmaTRise_m(0.0),
    sigmaTFall_m(0.0),
    tPulseLengthFWHM_m(0.0),
    correlationMatrix_m(getUnit6x6()),
    laserProfileFileName_m(""),
    laserImageName_m(""),
    laserIntensityCut_m(0.0),
    laserProfile_m(NULL),
    darkCurrentParts_m(0),
    darkInwardMargin_m(0.0),
    eInitThreshold_m(0.0),
    workFunction_m(0.0),
    fieldEnhancement_m(0.0),
    fieldThrFN_m(0.0),
    maxFN_m(0),
    paraFNA_m(0.0),
    paraFNB_m(0.0),
    paraFNY_m(0.0),
    paraFNVYSe_m(0.0),
    paraFNVYZe_m(0.0),
    secondaryFlag_m(0),
    ppVw_m(0.0),
    vVThermal_m(0.0),
    referencePz_m(0.0),
    referenceZ_m(0.0),
    avrgpz_m(0.0),
    I_m(0.0),
    E_m(0.0),
    bega_m(0.0),
    M_m(0.0)
{
    SetAttributes();

    Distribution *defaultDistribution = clone("UNNAMED_Distribution");
    defaultDistribution->builtin = true;

    try {
        OpalData::getInstance()->define(defaultDistribution);
    } catch(...) {
        delete defaultDistribution;
    }

    SetFieldEmissionParameters();
}
/**
 *
 *
 * @param name
 * @param parent
 */
Distribution::Distribution(const std::string &name, Distribution *parent):
    Definition(name, parent),
    distT_m(parent->distT_m),
    distrTypeT_m(DistrTypeT::NODIST),
    numberOfDistributions_m(parent->numberOfDistributions_m),
    emitting_m(parent->emitting_m),
    scan_m(parent->scan_m),
    particleRefData_m(parent->particleRefData_m),
    addedDistributions_m(parent->addedDistributions_m),
    particlesPerDist_m(parent->particlesPerDist_m),
    emissionModel_m(parent->emissionModel_m),
    tEmission_m(parent->tEmission_m),
    tBin_m(parent->tBin_m),
    currentEmissionTime_m(parent->currentEmissionTime_m),
    currentEnergyBin_m(parent->currentEmissionTime_m),
    currentSampleBin_m(parent->currentSampleBin_m),
    numberOfEnergyBins_m(parent->numberOfEnergyBins_m),
    numberOfSampleBins_m(parent->numberOfSampleBins_m),
    energyBins_m(NULL),
    energyBinHist_m(NULL),
    randGenEmit_m(parent->randGenEmit_m),
    pTotThermal_m(parent->pTotThermal_m),
    pmean_m(parent->pmean_m),
    cathodeWorkFunc_m(parent->cathodeWorkFunc_m),
    laserEnergy_m(parent->laserEnergy_m),
    cathodeFermiEnergy_m(parent->cathodeFermiEnergy_m),
    cathodeTemp_m(parent->cathodeTemp_m),
    emitEnergyUpperLimit_m(parent->emitEnergyUpperLimit_m),
    xDist_m(parent->xDist_m),
    pxDist_m(parent->pxDist_m),
    yDist_m(parent->yDist_m),
    pyDist_m(parent->pyDist_m),
    tOrZDist_m(parent->tOrZDist_m),
    pzDist_m(parent->pzDist_m),
    xWrite_m(parent->xWrite_m),
    pxWrite_m(parent->pxWrite_m),
    yWrite_m(parent->yWrite_m),
    pyWrite_m(parent->pyWrite_m),
    tOrZWrite_m(parent->tOrZWrite_m),
    pzWrite_m(parent->pzWrite_m),
    avrgpz_m(parent->avrgpz_m),
    inputMoUnits_m(parent->inputMoUnits_m),
    sigmaTRise_m(parent->sigmaTRise_m),
    sigmaTFall_m(parent->sigmaTFall_m),
    tPulseLengthFWHM_m(parent->tPulseLengthFWHM_m),
    sigmaR_m(parent->sigmaR_m),
    sigmaP_m(parent->sigmaP_m),
    cutoffR_m(parent->cutoffR_m),
    cutoffP_m(parent->cutoffP_m),
    correlationMatrix_m(parent->correlationMatrix_m),
    laserProfileFileName_m(parent->laserProfileFileName_m),
    laserImageName_m(parent->laserImageName_m),
    laserIntensityCut_m(parent->laserIntensityCut_m),
    laserProfile_m(NULL),
    darkCurrentParts_m(parent->darkCurrentParts_m),
    darkInwardMargin_m(parent->darkInwardMargin_m),
    eInitThreshold_m(parent->eInitThreshold_m),
    workFunction_m(parent->workFunction_m),
    fieldEnhancement_m(parent->fieldEnhancement_m),
    fieldThrFN_m(parent->fieldThrFN_m),
    maxFN_m(parent->maxFN_m),
    paraFNA_m(parent-> paraFNA_m),
    paraFNB_m(parent-> paraFNB_m),
    paraFNY_m(parent-> paraFNY_m),
    paraFNVYSe_m(parent-> paraFNVYSe_m),
    paraFNVYZe_m(parent-> paraFNVYZe_m),
    secondaryFlag_m(parent->secondaryFlag_m),
    ppVw_m(parent->ppVw_m),
    vVThermal_m(parent->vVThermal_m),
    tRise_m(parent->tRise_m),
    tFall_m(parent->tFall_m),
    sigmaRise_m(parent->sigmaRise_m),
    sigmaFall_m(parent->sigmaFall_m),
    cutoff_m(parent->cutoff_m),
    I_m(parent->I_m),
    E_m(parent->E_m),
    bega_m(parent->bega_m),
    M_m(parent->M_m)
{
}

Distribution::~Distribution() {

    if((Ippl::getNodes() == 1) && (os_m.is_open()))
        os_m.close();

    if (energyBins_m != NULL) {
        delete energyBins_m;
        energyBins_m = NULL;
    }

    if (energyBinHist_m != NULL) {
        gsl_histogram_free(energyBinHist_m);
        energyBinHist_m = NULL;
    }

    if (randGenEmit_m != NULL) {
        gsl_rng_free(randGenEmit_m);
        randGenEmit_m = NULL;
    }

    if(laserProfile_m) {
        delete laserProfile_m;
        laserProfile_m = NULL;
    }

}
/*
  void Distribution::printSigma(SigmaGenerator<double,unsigned int>::matrix_type& M, Inform& out) {
  for(int i=0; i<M.size1(); ++i) {
  for(int j=0; j<M.size2(); ++j) {
  *gmsg  << M(i,j) << " ";
  }
  *gmsg << endl;
  }
  }
*/

/**
 * At the moment only write the header into the file dist.dat
 * PartBunch will then append (very uggly)
 * @param
 * @param
 * @param
 */
void Distribution::WriteToFile() {
    /*
      if(Ippl::getNodes() == 1) {
      if(os_m.is_open()) {
      ;
      } else {
      *gmsg << " Write distribution to file data/dist.dat" << endl;
      std::string file("data/dist.dat");
      os_m.open(file.c_str());
      if(os_m.bad()) {
      *gmsg << "Unable to open output file " <<  file << endl;
      }
      os_m << "# x y ti px py pz "  << std::endl;
      os_m.close();
      }
      }
    */
}

/// Distribution can only be replaced by another distribution.
bool Distribution::canReplaceBy(Object *object) {
    return dynamic_cast<Distribution *>(object) != 0;
}

Distribution *Distribution::clone(const std::string &name) {
    return new Distribution(name, this);
}

void Distribution::execute() {
}

void Distribution::update() {
}

void Distribution::Create(size_t &numberOfParticles, double massIneV) {

    SetFieldEmissionParameters();

    switch (distrTypeT_m) {

    case DistrTypeT::MATCHEDGAUSS:
        CreateMatchedGaussDistribution(numberOfParticles, massIneV);
        break;
    case DistrTypeT::FROMFILE:
        CreateDistributionFromFile(numberOfParticles, massIneV);
        break;
    case DistrTypeT::GAUSS:
        CreateDistributionGauss(numberOfParticles, massIneV);
        break;
    case DistrTypeT::BINOMIAL:
        CreateDistributionBinomial(numberOfParticles, massIneV);
        break;
    case DistrTypeT::FLATTOP:
        CreateDistributionFlattop(numberOfParticles, massIneV);
        break;
    case DistrTypeT::GUNGAUSSFLATTOPTH:
        CreateDistributionFlattop(numberOfParticles, massIneV);
        break;
    case DistrTypeT::ASTRAFLATTOPTH:
        CreateDistributionFlattop(numberOfParticles, massIneV);
        break;
    default:
        INFOMSG("Distribution unknown." << endl;);
        break;
    }

    if (emitting_m) {
        unsigned int numAdditionalRNsPerParticle = 0;
        if (emissionModel_m == EmissionModelT::ASTRA) {
            numAdditionalRNsPerParticle = 2;
        }

        if (Options::cZero) {
            numAdditionalRNsPerParticle *= 2;
        }

        int saveProcessor = -1;

        for (size_t partIndex = 0; partIndex < numberOfParticles; ++ partIndex) {

            // Save to each processor in turn.
            ++ saveProcessor;
            if (saveProcessor >= Ippl::getNodes())
                saveProcessor = 0;

            if (Ippl::myNode() == saveProcessor) {
                std::vector<double> rns;
                for (unsigned int i = 0; i < numAdditionalRNsPerParticle; ++ i) {
                    double x = gsl_rng_uniform(randGenEmit_m);
                    rns.push_back(x);
                }
                additionalRNs_m.push_back(rns);
            } else {
                for (unsigned int i = 0; i < numAdditionalRNsPerParticle; ++ i) {
                    gsl_rng_uniform(randGenEmit_m);
                }
            }
        }
    }

    // AAA Scale and shift coordinates according to distribution input.
    ScaleDistCoordinates();
}

void  Distribution::CreatePriPart(PartBunch *beam, BoundaryGeometry &bg) {

    if(Options::ppdebug) {  // This is Parallel Plate Benchmark.
        int pc = 0;
        size_t lowMark = beam->getLocalNum();
        double vw = this->GetVw();
        double vt = this->GetvVThermal();
        double f_max = vw / vt * exp(-0.5);
        double test_a = vt / vw;
        double test_asq = test_a * test_a;
        size_t count = 0;
        size_t N_mean = static_cast<size_t>(floor(bg.getN() / Ippl::getNodes()));
        size_t N_extra = static_cast<size_t>(bg.getN() - N_mean * Ippl::getNodes());
        if(Ippl::myNode() == 0)
            N_mean += N_extra;
        if(bg.getN() != 0) {
            for(size_t i = 0; i < bg.getN(); i++) {
                if(pc == Ippl::myNode()) {
                    if(count < N_mean) {
                        /*==============Parallel Plate Benchmark=====================================*/
                        double test_s = 1;
                        double f_x = 0;
                        double test_x = 0;
                        while(test_s > f_x) {
                            test_s = IpplRandom();
                            test_s *= f_max;
                            test_x = IpplRandom();
                            test_x *= 10 * test_a; //range for normalized emission speed(0,10*test_a);
                            f_x = test_x / test_asq * exp(-test_x * test_x / 2 / test_asq);
                        }
                        double v_emi = test_x * vw;

                        double betaemit = v_emi / Physics::c;
                        double betagamma = betaemit / sqrt(1 - betaemit * betaemit);
                        /*============================================================================ */
                        beam->create(1);
                        if(pc != 0) {
                            beam->R[lowMark + count] = bg.getCooridinate(Ippl::myNode() * N_mean + count + N_extra);
                            beam->P[lowMark + count] = betagamma * bg.getMomenta(Ippl::myNode() * N_mean + count);
                        } else {
                            beam->R[lowMark + count] = bg.getCooridinate(count);
                            beam->P[lowMark + count] = betagamma * bg.getMomenta(count);
                        }
                        beam->Bin[lowMark + count] = 0;
                        beam->PType[lowMark + count] = ParticleType::REGULAR;
                        beam->TriID[lowMark + count] = 0;
                        beam->Q[lowMark + count] = beam->getChargePerParticle();
                        beam->LastSection[lowMark + count] = 0;
                        beam->Ef[lowMark + count] = Vector_t(0.0);
                        beam->Bf[lowMark + count] = Vector_t(0.0);
                        beam->dt[lowMark + count] = beam->getdT();
                        count ++;
                    }
                }
                pc++;
                if(pc == Ippl::getNodes())
                    pc = 0;
            }
            bg.clearCooridinateArray();
            bg.clearMomentaArray();
            beam->boundp();
        }
        *gmsg << *beam << endl;

    } else {// Normal procedure to create primary particles

        int pc = 0;
        size_t lowMark = beam->getLocalNum();
        size_t count = 0;
        size_t N_mean = static_cast<size_t>(floor(bg.getN() / Ippl::getNodes()));
        size_t N_extra = static_cast<size_t>(bg.getN() - N_mean * Ippl::getNodes());

        if(Ippl::myNode() == 0)
            N_mean += N_extra;
        if(bg.getN() != 0) {
            for(size_t i = 0; i < bg.getN(); i++) {

                if(pc == Ippl::myNode()) {
                    if(count < N_mean) {
                        beam->create(1);
                        if(pc != 0)
                            beam->R[lowMark + count] = bg.getCooridinate(Ippl::myNode() * N_mean + count + N_extra); // node 0 will emit the particle with coordinate ID from 0 to N_mean+N_extra, so other nodes should shift to node_number*N_mean+N_extra
                        else
                            beam->R[lowMark + count] = bg.getCooridinate(count); // for node0 the particle number N_mean =  N_mean + N_extra
                        beam->P[lowMark + count] = Vector_t(0.0);
                        beam->Bin[lowMark + count] = 0;
                        beam->PType[lowMark + count] = ParticleType::REGULAR;
                        beam->TriID[lowMark + count] = 0;
                        beam->Q[lowMark + count] = beam->getChargePerParticle();
                        beam->LastSection[lowMark + count] = 0;
                        beam->Ef[lowMark + count] = Vector_t(0.0);
                        beam->Bf[lowMark + count] = Vector_t(0.0);
                        beam->dt[lowMark + count] = beam->getdT();
                        count++;

                    }

                }
                pc++;
                if(pc == Ippl::getNodes())
                    pc = 0;

            }

        }
        bg.clearCooridinateArray();
        beam->boundp();//fixme if bg.getN()==0?
    }
    *gmsg << *beam << endl;
}

void Distribution::DoRestartOpalT(PartBunch &beam, size_t Np, int restartStep, H5PartWrapper *dataSource) {

    IpplTimings::startTimer(beam.distrReload_m);

    long numParticles = dataSource->getNumParticles();
    size_t numParticlesPerNode = numParticles / Ippl::getNodes();

    size_t firstParticle = numParticlesPerNode * Ippl::myNode();
    size_t lastParticle = firstParticle + numParticlesPerNode - 1;
    if (Ippl::myNode() == Ippl::getNodes() - 1)
        lastParticle = numParticles - 1;

    numParticles = lastParticle - firstParticle + 1;
    PAssert(numParticles >= 0);

    beam.create(numParticles);

    dataSource->readHeader();
    dataSource->readStep(beam, firstParticle, lastParticle);

    beam.boundp();

    double actualT = beam.getT();
    long long ltstep = beam.getLocalTrackStep();
    long long gtstep = beam.getGlobalTrackStep();

    IpplTimings::stopTimer(beam.distrReload_m);

    *gmsg << "Total number of particles in the h5 file= " << beam.getTotalNum() << "\n"
          << "Global step= " << gtstep << "; Local step= " << ltstep << "; "
          << "restart step= " << restartStep << "\n"
          << "time of restart= " << actualT << "; phishift= " << OpalData::getInstance()->getGlobalPhaseShift() << endl;
}

void Distribution::DoRestartOpalCycl(PartBunch &beam,
                                     size_t Np,
                                     int restartStep,
                                     const int specifiedNumBunch,
                                     H5PartWrapper *dataSource) {

    // h5_int64_t rc;
    IpplTimings::startTimer(beam.distrReload_m);
    INFOMSG("---------------- Start reading hdf5 file----------------" << endl);

    long numParticles = dataSource->getNumParticles();
    size_t numParticlesPerNode = numParticles / Ippl::getNodes();

    size_t firstParticle = numParticlesPerNode * Ippl::myNode();
    size_t lastParticle = firstParticle + numParticlesPerNode - 1;
    if (Ippl::myNode() == Ippl::getNodes() - 1)
        lastParticle = numParticles - 1;

    numParticles = lastParticle - firstParticle + 1;
    PAssert(numParticles >= 0);

    beam.create(numParticles);

    dataSource->readHeader();
    dataSource->readStep(beam, firstParticle, lastParticle);

    beam.Q = beam.getChargePerParticle();

    beam.boundp();

    double meanE = static_cast<H5PartWrapperForPC*>(dataSource)->getMeanKineticEnergy();

    const int globalN = beam.getTotalNum();
    INFOMSG("Restart from hdf5 format file " << OpalData::getInstance()->getInputBasename() << ".h5" << endl);
    INFOMSG("total number of particles = " << globalN << endl);
    INFOMSG("* Restart Energy = " << meanE << " (MeV), Path lenght = " << beam.getLPath() << " (m)" <<  endl);
    INFOMSG("Tracking Step since last bunch injection is " << beam.getSteptoLastInj() << endl);
    INFOMSG(beam.getNumBunch() << " Bunches(bins) exist in this file" << endl);

    double gamma = 1 + meanE / beam.getM() * 1.0e6;
    double beta = sqrt(1.0 - (1.0 / std::pow(gamma, 2.0)));

    INFOMSG("* Gamma = " << gamma << ", Beta = " << beta << endl);

    if(dataSource->predecessorIsSameFlavour()) {
        INFOMSG("Restart from hdf5 file generated by OPAL-cycl" << endl);
        if(specifiedNumBunch > 1) {
            // the allowed maximal bin number is set to 1000
            beam.setPBins(new PartBinsCyc(1000, beam.getNumBunch()));
        }

    } else {
        INFOMSG("Restart from hdf5 file generated by OPAL-t" << endl);

        Vector_t meanR(0.0, 0.0, 0.0);
        Vector_t meanP(0.0, 0.0, 0.0);
        unsigned long int newLocalN = beam.getLocalNum();
        for(unsigned int i = 0; i < newLocalN; ++i) {
            for(int d = 0; d < 3; ++d) {
                meanR(d) += beam.R[i](d);
                meanP(d) += beam.P[i](d);
            }
        }
        reduce(meanR, meanR, OpAddAssign());
        meanR /= Vector_t(globalN);
        reduce(meanP, meanP, OpAddAssign());
        meanP /= Vector_t(globalN);
        INFOMSG("Rmean = " << meanR << "[m], Pmean=" << meanP << endl);

        for(unsigned int i = 0; i < newLocalN; ++i) {
            beam.R[i] -= meanR;
            beam.P[i] -= meanP;
        }
    }

    INFOMSG("---------------Finished reading hdf5 file---------------" << endl);
    IpplTimings::stopTimer(beam.distrReload_m);
}

void Distribution::DoRestartOpalE(EnvelopeBunch &beam, size_t Np, int restartStep,
                                  H5PartWrapper *dataSource) {
    IpplTimings::startTimer(beam.distrReload_m);
    int N = dataSource->getNumParticles();
    *gmsg << "total number of slices = " << N << endl;

    beam.distributeSlices(N);
    beam.createBunch();
    long long starti = beam.mySliceStartOffset();
    long long endi = beam.mySliceEndOffset();

    dataSource->readHeader();
    dataSource->readStep(beam, starti, endi);

    beam.setCharge(beam.getChargePerParticle());
    IpplTimings::stopTimer(beam.distrReload_m);
}

Distribution *Distribution::find(const std::string &name) {
    Distribution *dist = dynamic_cast<Distribution *>(OpalData::getInstance()->find(name));

    if(dist == 0) {
        throw OpalException("Distribution::find()", "Distribution \"" + name + "\" not found.");
    }

    return dist;
}

double Distribution::GetTEmission() {
    if(tEmission_m > 0.0) {
        return tEmission_m;
    }

    distT_m = Attributes::getString(itsAttr[AttributesT::DISTRIBUTION]);
    if(distT_m == "GAUSS")
        distrTypeT_m = DistrTypeT::GAUSS;
    else if(distT_m == "GUNGAUSSFLATTOPTH")
        distrTypeT_m = DistrTypeT::GUNGAUSSFLATTOPTH;
    else if(distT_m == "FROMFILE")
        distrTypeT_m = DistrTypeT::FROMFILE;
    else if(distT_m == "BINOMIAL")
        distrTypeT_m = DistrTypeT::BINOMIAL;

    tPulseLengthFWHM_m = Attributes::getReal(itsAttr[AttributesT::TPULSEFWHM]);
    cutoff_m = Attributes::getReal(itsAttr[ LegacyAttributesT::CUTOFF]);
    tRise_m = Attributes::getReal(itsAttr[AttributesT::TRISE]);
    tFall_m = Attributes::getReal(itsAttr[AttributesT::TFALL]);
    double tratio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));
    sigmaRise_m = tRise_m / tratio;
    sigmaFall_m = tFall_m / tratio;

    switch(distrTypeT_m) {
    case DistrTypeT::ASTRAFLATTOPTH: {
        double a = tPulseLengthFWHM_m / 2;
        double sig = tRise_m / 2;
        double inv_erf08 = 0.906193802436823; // erfinv(0.8)
        double sqr2 = sqrt(2.);
        double t = a - sqr2 * sig * inv_erf08;
        double tmps = sig;
        double tmpt = t;
        for(int i = 0; i < 10; ++ i) {
            sig = (t + tRise_m - a) / (sqr2 * inv_erf08);
            t = a - 0.5 * sqr2 * (sig + tmps) * inv_erf08;
            sig = (0.5 * (t + tmpt) + tRise_m - a) / (sqr2 * inv_erf08);
            tmps = sig;
            tmpt = t;
        }
        tEmission_m = tPulseLengthFWHM_m + 10 * sig;
        break;
    }
    case DistrTypeT::GUNGAUSSFLATTOPTH: {
        tEmission_m = tPulseLengthFWHM_m + (cutoff_m - sqrt(2.0 * log(2.0))) * (sigmaRise_m + sigmaFall_m);
        break;
    }
    default:
        tEmission_m = 0.0;
    }
    return tEmission_m;
}

double Distribution::GetEkin() const {return Attributes::getReal(itsAttr[AttributesT::EKIN]);}
double Distribution::GetLaserEnergy() const {return Attributes::getReal(itsAttr[AttributesT::ELASER]);}
double Distribution::GetWorkFunctionRf() const {return Attributes::getReal(itsAttr[AttributesT::W]);}

size_t Distribution::GetNumberOfDarkCurrentParticles() { return (size_t) Attributes::getReal(itsAttr[AttributesT::NPDARKCUR]);}
double Distribution::GetDarkCurrentParticlesInwardMargin() { return Attributes::getReal(itsAttr[AttributesT::INWARDMARGIN]);}
double Distribution::GetEInitThreshold() { return Attributes::getReal(itsAttr[AttributesT::EINITHR]);}
double Distribution::GetWorkFunction() { return Attributes::getReal(itsAttr[AttributesT::FNPHIW]); }
double Distribution::GetFieldEnhancement() { return Attributes::getReal(itsAttr[AttributesT::FNBETA]); }
size_t Distribution::GetMaxFNemissionPartPerTri() { return (size_t) Attributes::getReal(itsAttr[AttributesT::FNMAXEMI]);}
double Distribution::GetFieldFNThreshold() { return Attributes::getReal(itsAttr[AttributesT::FNFIELDTHR]);}
double Distribution::GetFNParameterA() { return Attributes::getReal(itsAttr[AttributesT::FNA]);}
double Distribution::GetFNParameterB() { return Attributes::getReal(itsAttr[AttributesT::FNB]);}
double Distribution::GetFNParameterY() { return Attributes::getReal(itsAttr[AttributesT::FNY]);}
double Distribution::GetFNParameterVYZero() { return Attributes::getReal(itsAttr[AttributesT::FNVYZERO]);}
double Distribution::GetFNParameterVYSecond() { return Attributes::getReal(itsAttr[AttributesT::FNVYSECOND]);}
int    Distribution::GetSecondaryEmissionFlag() { return Attributes::getReal(itsAttr[AttributesT::SECONDARYFLAG]);}
bool   Distribution::GetEmissionMode() { return Attributes::getBool(itsAttr[AttributesT::NEMISSIONMODE]);}

std::string Distribution::GetTypeofDistribution() { return (std::string) Attributes::getString(itsAttr[AttributesT::DISTRIBUTION]);}

double Distribution::GetvSeyZero() {return Attributes::getReal(itsAttr[AttributesT::VSEYZERO]);}// return sey_0 in Vaughan's model
double Distribution::GetvEZero() {return Attributes::getReal(itsAttr[AttributesT::VEZERO]);}// return the energy related to sey_0 in Vaughan's model
double Distribution::GetvSeyMax() {return Attributes::getReal(itsAttr[AttributesT::VSEYMAX]);}// return sey max in Vaughan's model
double Distribution::GetvEmax() {return Attributes::getReal(itsAttr[AttributesT::VEMAX]);}// return Emax in Vaughan's model
double Distribution::GetvKenergy() {return Attributes::getReal(itsAttr[AttributesT::VKENERGY]);}// return fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
double Distribution::GetvKtheta() {return Attributes::getReal(itsAttr[AttributesT::VKTHETA]);}// return fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
double Distribution::GetvVThermal() {return Attributes::getReal(itsAttr[AttributesT::VVTHERMAL]);}// thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
double Distribution::GetVw() {return Attributes::getReal(itsAttr[AttributesT::VW]);}// velocity scalar for parallel plate benchmark;

int Distribution::GetSurfMaterial() {return (int)Attributes::getReal(itsAttr[AttributesT::SURFMATERIAL]);}// Surface material number for Furman-Pivi's Model;

Inform &Distribution::printInfo(Inform &os) const {

    os << "* ************* D I S T R I B U T I O N ********************************************" << endl;
    os << "* " << endl;
    if (OpalData::getInstance()->inRestartRun()) {
        os << "* In restart. Distribution read in from .h5 file." << endl;
    } else {
        if (addedDistributions_m.size() > 0)
            os << "* Main Distribution" << endl
               << "-----------------" << endl;

        if (particlesPerDist_m.empty())
            PrintDist(os, 0);
        else
            PrintDist(os, particlesPerDist_m.at(0));

        size_t distCount = 1;
        for (unsigned distIndex = 0; distIndex < addedDistributions_m.size(); distIndex++) {
            os << "* " << endl;
            os << "* Added Distribution #" << distCount << endl;
            os << "* ----------------------" << endl;
            addedDistributions_m.at(distIndex)->PrintDist(os, particlesPerDist_m.at(distCount));
            distCount++;
        }

        os << "* " << endl;
        if (numberOfEnergyBins_m > 0) {
            os << "* Number of energy bins    = " << numberOfEnergyBins_m << endl;

            //            if (numberOfEnergyBins_m > 1)
            //    PrintEnergyBins(os);
        }

        if (emitting_m) {
            os << "* Distribution is emitted. " << endl;
            os << "* Emission time            = " << tEmission_m << " [sec]" << endl;
            os << "* Time per bin             = " << tEmission_m / numberOfEnergyBins_m << " [sec]" << endl;
            os << "* Delta t during emission  = " << tBin_m / numberOfSampleBins_m << " [sec]" << endl;
            os << "* " << endl;
            PrintEmissionModel(os);
            os << "* " << endl;
        } else
            os << "* Distribution is injected." << endl;
    }
    os << "* " << endl;
    os << "* **********************************************************************************" << endl;

    return os;
}

const PartData &Distribution::GetReference() const {
    // Cast away const, to allow logically constant Distribution to update.
    const_cast<Distribution *>(this)->update();
    return particleRefData_m;
}

void Distribution::AddDistributions() {
    /*
     * Move particle coordinates from added distributions to main distribution.
     */

    std::vector<Distribution *>::iterator addedDistIt;
    for (addedDistIt = addedDistributions_m.begin();
         addedDistIt != addedDistributions_m.end(); addedDistIt++) {

        std::vector<double>::iterator distIt;
        for (distIt = (*addedDistIt)->GetXDist().begin();
             distIt != (*addedDistIt)->GetXDist().end();
             distIt++) {
            xDist_m.push_back(*distIt);
        }
        (*addedDistIt)->EraseXDist();

        for (distIt = (*addedDistIt)->GetBGxDist().begin();
             distIt != (*addedDistIt)->GetBGxDist().end();
             distIt++) {
            pxDist_m.push_back(*distIt);
        }
        (*addedDistIt)->EraseBGxDist();

        for (distIt = (*addedDistIt)->GetYDist().begin();
             distIt != (*addedDistIt)->GetYDist().end();
             distIt++) {
            yDist_m.push_back(*distIt);
        }
        (*addedDistIt)->EraseYDist();

        for (distIt = (*addedDistIt)->GetBGyDist().begin();
             distIt != (*addedDistIt)->GetBGyDist().end();
             distIt++) {
            pyDist_m.push_back(*distIt);
        }
        (*addedDistIt)->EraseBGyDist();

        for (distIt = (*addedDistIt)->GetTOrZDist().begin();
             distIt != (*addedDistIt)->GetTOrZDist().end();
             distIt++) {
            tOrZDist_m.push_back(*distIt);
        }
        (*addedDistIt)->EraseTOrZDist();

        for (distIt = (*addedDistIt)->GetBGzDist().begin();
             distIt != (*addedDistIt)->GetBGzDist().end();
             distIt++) {
            pzDist_m.push_back(*distIt);
        }
        (*addedDistIt)->EraseBGzDist();
    }
}

void Distribution::ApplyEmissionModel(double eZ, double &px, double &py, double &pz, const std::vector<double> &additionalRNs) {

    switch (emissionModel_m) {

    case EmissionModelT::NONE:
        ApplyEmissModelNone(pz);
        break;
    case EmissionModelT::ASTRA:
        ApplyEmissModelAstra(px, py, pz, additionalRNs);
        break;
    case EmissionModelT::NONEQUIL:
        ApplyEmissModelNonEquil(eZ, px, py, pz);
        break;
    default:
        break;
    }
}

void Distribution::ApplyEmissModelAstra(double &px, double &py, double &pz, const std::vector<double> &additionalRNs) {

    double phi = 2.0 * acos(sqrt(additionalRNs[0]));
    double theta = 2.0 * Physics::pi * additionalRNs[1];

    px = pTotThermal_m * sin(phi) * cos(theta);
    py = pTotThermal_m * sin(phi) * sin(theta);
    pz = pTotThermal_m * std::abs(cos(phi));

}

void Distribution::ApplyEmissModelNone(double &pz) {
    pz += pTotThermal_m;
}

void Distribution::ApplyEmissModelNonEquil(double eZ,
                                           double &bgx,
                                           double &bgy,
                                           double &bgz) {

    double phiEffective = cathodeWorkFunc_m
        - sqrt(Physics::q_e * eZ /
               (4.0 * Physics::pi * Physics::epsilon_0));
    double lowEnergyLimit = cathodeFermiEnergy_m + phiEffective - laserEnergy_m;

    // Generate emission energy.
    double energy = 0.0;
    bool allow = false;
    while (!allow) {
        energy = lowEnergyLimit + (gsl_rng_uniform(randGenEmit_m)*emitEnergyUpperLimit_m);
        double randFuncValue = gsl_rng_uniform(randGenEmit_m);
        double funcValue = (1.0
                            - 1.0
                            / (1.0
                               + exp((energy + laserEnergy_m - cathodeFermiEnergy_m)
                                     / (Physics::kB * cathodeTemp_m))))
            / (1.0
               + exp((energy - cathodeFermiEnergy_m)
                     / (Physics::kB * cathodeTemp_m)));
        if (randFuncValue <= funcValue)
            allow = true;
    }

    // Compute emission angles.
    double energyInternal = energy + laserEnergy_m;
    double energyExternal = energy + laserEnergy_m
        - cathodeFermiEnergy_m - phiEffective;

    double thetaInMax = acos(sqrt((cathodeFermiEnergy_m + phiEffective)
                                  / (energy + laserEnergy_m)));
    double thetaIn = gsl_rng_uniform(randGenEmit_m)*thetaInMax;
    double sinThetaOut = sin(thetaIn) * sqrt(energyInternal / energyExternal);
    double phi = Physics::two_pi * gsl_rng_uniform(randGenEmit_m);

    // Compute emission momenta.
    double betaGammaExternal
        = sqrt(pow(energyExternal / (Physics::m_e * 1.0e9) + 1.0, 2.0) - 1.0);

    bgx = betaGammaExternal * sinThetaOut * cos(phi);
    bgy = betaGammaExternal * sinThetaOut * sin(phi);
    bgz = betaGammaExternal * sqrt(1.0 - pow(sinThetaOut, 2.0));

}

void Distribution::CalcPartPerDist(size_t numberOfParticles) {

    typedef std::vector<Distribution *>::iterator iterator;

    if (numberOfDistributions_m == 1)
        particlesPerDist_m.push_back(numberOfParticles);
    else {
        double totalWeight = GetWeight();
        for (iterator it = addedDistributions_m.begin(); it != addedDistributions_m.end(); it++) {
            totalWeight += (*it)->GetWeight();
        }

        particlesPerDist_m.push_back(0);
        size_t numberOfCommittedPart = 0;
        for (iterator it = addedDistributions_m.begin(); it != addedDistributions_m.end(); it++) {
            size_t particlesCurrentDist = numberOfParticles * (*it)->GetWeight() / totalWeight;
            particlesPerDist_m.push_back(particlesCurrentDist);
            numberOfCommittedPart += particlesCurrentDist;
        }

        // Remaining particles go into main distribution.
        particlesPerDist_m.at(0) = numberOfParticles - numberOfCommittedPart;

    }

}

void Distribution::CheckEmissionParameters() {

    // There must be at least on energy bin for an emitted beam.
    numberOfEnergyBins_m
        = std::abs(static_cast<int> (Attributes::getReal(itsAttr[AttributesT::NBIN])));
    if (numberOfEnergyBins_m == 0)
        numberOfEnergyBins_m = 1;

    int emissionSteps = std::abs(static_cast<int> (Attributes::getReal(itsAttr[AttributesT::EMISSIONSTEPS])));

    // There must be at least one emission step.
    if (emissionSteps == 0)
        emissionSteps = 1;

    // Set number of sample bins per energy bin from the number of emission steps.
    numberOfSampleBins_m = static_cast<int> (std::ceil(emissionSteps / numberOfEnergyBins_m));
    if (numberOfSampleBins_m <= 0)
        numberOfSampleBins_m = 1;

    // Initialize emission counters.
    currentEnergyBin_m = 1;
    currentSampleBin_m = 0;

}

void Distribution::CheckIfEmitted() {

    emitting_m = Attributes::getBool(itsAttr[AttributesT::EMITTED]);

    switch (distrTypeT_m) {

    case DistrTypeT::ASTRAFLATTOPTH:
        emitting_m = true;
        break;
    case DistrTypeT::GUNGAUSSFLATTOPTH:
        emitting_m = true;
        break;
    default:
        break;
    }
}

void Distribution::CheckParticleNumber(size_t &numberOfParticles) {

    size_t numberOfDistParticles = tOrZDist_m.size();
    reduce(numberOfDistParticles, numberOfDistParticles, OpAddAssign());

    if (numberOfDistParticles != numberOfParticles) {
        *gmsg << "\n--------------------------------------------------" << endl
              << "Error!! The number of particles in the initial" << endl
              << "distribution is " << numberOfDistParticles << "." << endl << endl
              << "This is different from the number of particles" << endl
              << "defined by the BEAM command: " << numberOfParticles << endl << endl
              << "This is often happens when using a FROMFILE type" << endl
              << "distribution and not matching the number of" << endl
              << "particles in the particles file(s) with the number" << endl
              << "given in the BEAM command." << endl << endl
              << "The number of particles in the initial distribution" << endl
              << "(" << numberOfDistParticles << ") "
              << "will take precedence." << endl
              << "---------------------------------------------------\n" << endl;
        throw OpalException("Distribution::CheckParticleNumber",
                            "Number of macro particles and NPART on BEAM are not equal");
    }
    numberOfParticles = numberOfDistParticles;
}

void Distribution::ChooseInputMomentumUnits(InputMomentumUnitsT::InputMomentumUnitsT inputMoUnits) {

    /*
     * Toggle what units to use for inputing momentum.
     */
    std::string inputUnits = Attributes::getString(itsAttr[AttributesT::INPUTMOUNITS]);
    if (inputUnits == "NONE")
        inputMoUnits_m = InputMomentumUnitsT::NONE;
    else if (inputUnits == "EV")
        inputMoUnits_m = InputMomentumUnitsT::EV;
    else
        inputMoUnits_m = inputMoUnits;

}

double Distribution::ConvertBetaGammaToeV(double valueInBetaGamma, double massIneV) {
    if (valueInBetaGamma < 0)
        return -1.0 * (sqrt(pow(valueInBetaGamma, 2.0) + 1.0) - 1.0) * massIneV;
    else
        return (sqrt(pow(valueInBetaGamma, 2.0) + 1.0) - 1.0) * massIneV;
}

double Distribution::ConverteVToBetaGamma(double valueIneV, double massIneV) {
    if (valueIneV < 0)
        return -1.0 * sqrt( pow( -1.0 * valueIneV / massIneV + 1.0, 2.0) - 1.0);
    else
        return sqrt( pow( valueIneV / massIneV + 1.0, 2.0) - 1.0);
}

double Distribution::ConvertMeVPerCToBetaGamma(double valueInMeVPerC, double massIneV) {
    return sqrt(pow(valueInMeVPerC * 1.0e6 * Physics::c / massIneV + 1.0, 2.0) - 1.0);
}

void Distribution::CreateDistributionBinomial(size_t numberOfParticles, double massIneV) {

    SetDistParametersBinomial(massIneV);
    GenerateBinomial(numberOfParticles);
}

void Distribution::CreateDistributionFlattop(size_t numberOfParticles, double massIneV) {

    SetDistParametersFlattop(massIneV);

    if (emitting_m) {
        if (laserProfile_m == NULL)
            GenerateFlattopT(numberOfParticles);
        else
            GenerateFlattopLaserProfile(numberOfParticles);
    } else
        GenerateFlattopZ(numberOfParticles);
}

void Distribution::CreateDistributionFromFile(size_t numberOfParticles, double massIneV) {

    *gmsg << level3 << "\n"
          << "------------------------------------------------------------------------------------\n";
    *gmsg << "READ INITIAL DISTRIBUTION FROM FILE \""
          << Attributes::getString(itsAttr[AttributesT::FNAME])
          << "\"\n";
    *gmsg << "------------------------------------------------------------------------------------\n" << endl;

    // Data input file is only read by node 0.
    std::ifstream inputFile;
    size_t numberOfParticlesRead = 0;
    if (Ippl::myNode() == 0) {
        std::string fileName = Attributes::getString(itsAttr[AttributesT::FNAME]);
        inputFile.open(fileName.c_str());
        if (inputFile.fail())
            throw OpalException("Distribution::Create()",
                                "Open file operation failed, please check if \""
                                + fileName
                                + "\" really exists.");
        else {
            int tempInt = 0;
            inputFile >> tempInt;
            numberOfParticlesRead = static_cast<size_t>(tempInt);
        }
    }
    reduce(numberOfParticlesRead, numberOfParticlesRead, OpAddAssign());

    /*
     * We read in the particle information using node zero, but save the particle
     * data to each processor in turn.
     */
    int saveProcessor = -1;

    for (size_t particleIndex = 0; particleIndex < numberOfParticlesRead; particleIndex++) {

        // Read particle.
        double x = 0.0;
        double px =0.0;
        double y = 0.0;
        double py = 0.0;
        double tOrZ = 0.0;
        double pz = 0.0;
        if (Ippl::myNode() == 0) {
            inputFile >> x
                      >> px
                      >> y
                      >> py
                      >> tOrZ
                      >> pz;

            switch (inputMoUnits_m) {

            case InputMomentumUnitsT::EV:
                px = ConverteVToBetaGamma(px, massIneV);
                py = ConverteVToBetaGamma(py, massIneV);
                pz = ConverteVToBetaGamma(pz, massIneV);
                break;

            default:
                break;
            }

            saveProcessor++;
            if (saveProcessor >= Ippl::getNodes())
                saveProcessor = 0;
        } else
            saveProcessor = 0;

        // Share data.
        reduce(saveProcessor, saveProcessor, OpAddAssign());
        reduce(x, x, OpAddAssign());
        reduce(px, px, OpAddAssign());
        reduce(y, y, OpAddAssign());
        reduce(py, py, OpAddAssign());
        reduce(tOrZ, tOrZ, OpAddAssign());
        reduce(pz, pz, OpAddAssign());

        // Store particles on proper processor.
        if (Ippl::myNode() == saveProcessor) {
            xDist_m.push_back(x);
            pxDist_m.push_back(px);
            yDist_m.push_back(y);
            pyDist_m.push_back(py);
            tOrZDist_m.push_back(tOrZ);
            pzDist_m.push_back(pz);
        }
    }
    if (Ippl::myNode() == 0)
        inputFile.close();
}


void Distribution::CreateMatchedGaussDistribution(size_t numberOfParticles, double massIneV) {

  /*
    ToDo:
    - add flag in order to calculate tunes from FMLOWE to FMHIGHE
    - eliminate physics and error
  */

  std::string LineName = Attributes::getString(itsAttr[AttributesT::LINE]);
  if (LineName != "") {
    const BeamSequence* LineSequence = BeamSequence::find(LineName);
    if (LineSequence != NULL) {
      SpecificElementVisitor<Cyclotron> CyclotronVisitor(*LineSequence->fetchLine());
      CyclotronVisitor.execute();
      size_t NumberOfCyclotrons = CyclotronVisitor.size();

      if (NumberOfCyclotrons > 1) {
	throw OpalException("Distribution::CreateMatchedGaussDistribution",
			    "I am confused, found more than one Cyclotron element in line");
      }
      if (NumberOfCyclotrons == 0) {
	throw OpalException("Distribution::CreateMatchedGaussDistribution",
			    "didn't find any Cyclotron element in line");
      }
      const Cyclotron* CyclotronElement = CyclotronVisitor.front();

      *gmsg << "* ----------------------------------------------------" << endl;
      *gmsg << "* About to find closed orbit and matched distribution " << endl;
      *gmsg << "* I= " << I_m*1E3 << " (mA)  E= " << E_m*1E-6 << " (MeV)" << endl;
      *gmsg << "* EX= "  << Attributes::getReal(itsAttr[AttributesT::EX])
	    << "* EY= " << Attributes::getReal(itsAttr[AttributesT::EY])
	    << "* ET= " << Attributes::getReal(itsAttr[AttributesT::ET])
	    << "* FMAPFN " << Attributes::getString(itsAttr[AttributesT::FMAPFN]) << endl //CyclotronElement->getFieldMapFN() << endl
	    << "* FMSYM= " << (int)Attributes::getReal(itsAttr[AttributesT::MAGSYM])
	    << "* HN= "   << CyclotronElement->getCyclHarm()
	    << "* PHIINIT= " << CyclotronElement->getPHIinit() << endl;
      *gmsg << "* ----------------------------------------------------" << endl;

      const double wo = CyclotronElement->getRfFrequ()*1E6/CyclotronElement->getCyclHarm()*2.0*Physics::pi;

      const double fmLowE  = CyclotronElement->getFMLowE();
      const double fmHighE = CyclotronElement->getFMHighE();

      double lE,hE;
      lE = fmLowE;
      hE = fmHighE;

      if ((lE<0) || (hE<0)) {
	lE = E_m*1E-6;
	hE = E_m*1E-6;
      }

      int Nint = 1000;
      bool writeMap = true;

      SigmaGenerator<double,unsigned int> *siggen = new SigmaGenerator<double,unsigned int>(I_m,
						 Attributes::getReal(itsAttr[AttributesT::EX])*1E6,
						 Attributes::getReal(itsAttr[AttributesT::EY])*1E6,
						 Attributes::getReal(itsAttr[AttributesT::ET])*1E6,
						 wo,
						 E_m*1E-6,
						 CyclotronElement->getCyclHarm(),
						 massIneV*1E-6,
						 lE,
						 hE,
						 (int)Attributes::getReal(itsAttr[AttributesT::MAGSYM]),
						 Nint,
						 Attributes::getString(itsAttr[AttributesT::FMAPFN]),
						 Attributes::getReal(itsAttr[AttributesT::ORDERMAPS]),
						 writeMap);


      if(siggen->match(Attributes::getReal(itsAttr[AttributesT::RESIDUUM]),
                       Attributes::getReal(itsAttr[AttributesT::MAXSTEPSSI]),
                       Attributes::getReal(itsAttr[AttributesT::MAXSTEPSCO]),
                       CyclotronElement->getPHIinit(),
                       Attributes::getReal(itsAttr[AttributesT::RGUESS]),
                       Attributes::getString(itsAttr[AttributesT::FMTYPE]),
                       false))  {

	std::array<double,3> Emit = siggen->getEmittances();

	if (Attributes::getReal(itsAttr[AttributesT::RGUESS]) > 0)
	  *gmsg << "* RGUESS " << Attributes::getReal(itsAttr[AttributesT::RGUESS])/1000.0 << " (m) " << endl;

	*gmsg << "* Converged (Ex, Ey, Ez) = (" << Emit[0] << ", " << Emit[1] << ", "
              << Emit[2] << ") pi mm mrad for E= " << E_m*1E-6 << " (MeV)" << endl;
	*gmsg << "* Sigma-Matrix " << endl;

	for(unsigned int i = 0; i < siggen->getSigma().size1(); ++ i) {
	  *gmsg << std::setprecision(4)  << std::setw(11) << siggen->getSigma()(i,0);
	  for(unsigned int j = 1; j < siggen->getSigma().size2(); ++ j) {
	    *gmsg << " & " <<  std::setprecision(4)  << std::setw(11) << siggen->getSigma()(i,j);
	  }
	  *gmsg << " \\\\" << endl;
	}

	/*

	  Now setup the distribution generator
	  Units of the Sigma Matrix:

	  spatial: mm
	  moment:  rad

	*/

	if(Options::cloTuneOnly)
	  throw OpalException("Do only CLO and tune calculation","... ");


	auto sigma = siggen->getSigma();
	// change units from mm to m
	for (unsigned int i = 0; i < 3; ++ i) {
	  for (unsigned int j = 0; j < 6; ++ j) {
	    sigma(2 * i, j) *= 1e-3;
	    sigma(j, 2 * i) *= 1e-3;
	  }
	}

	for (unsigned int i = 0; i < 3; ++ i) {
          if ( sigma(2 * i, 2 * i) < 0 || sigma(2 * i + 1, 2 * i + 1) < 0 )
              throw OpalException("Distribution::CreateMatchedGaussDistribution()",
                                  "Negative value on the diagonal of the sigma matrix.");

	  sigmaR_m[i] = std::sqrt(sigma(2 * i, 2 * i));
	  sigmaP_m[i] = std::sqrt(sigma(2 * i + 1, 2 * i + 1));


	}

	if (inputMoUnits_m == InputMomentumUnitsT::EV) {
	  for (unsigned int i = 0; i < 3; ++ i) {
	    sigmaP_m[i] = ConverteVToBetaGamma(sigmaP_m[i], massIneV);
	  }
	}

	correlationMatrix_m(1, 0) = sigma(0, 1) / (sqrt(sigma(0, 0) * sigma(1, 1)));
	correlationMatrix_m(3, 2) = sigma(2, 3) / (sqrt(sigma(2, 2) * sigma(3, 3)));
	correlationMatrix_m(5, 4) = sigma(4, 5) / (sqrt(sigma(4, 4) * sigma(5, 5)));
	correlationMatrix_m(4, 0) = sigma(0, 4) / (sqrt(sigma(0, 0) * sigma(4, 4)));
	correlationMatrix_m(4, 1) = sigma(1, 4) / (sqrt(sigma(1, 1) * sigma(4, 4)));
	correlationMatrix_m(5, 0) = sigma(0, 5) / (sqrt(sigma(0, 0) * sigma(5, 5)));
	correlationMatrix_m(5, 1) = sigma(1, 5) / (sqrt(sigma(1, 1) * sigma(5, 5)));

	CreateDistributionGauss(numberOfParticles, massIneV);
      }
      else {
	*gmsg << "* Not converged for " << E_m*1E-6 << " MeV" << endl;

        if (siggen)
            delete siggen;

        throw OpalException("Distribution::CreateMatchedGaussDistribution", "didn't find any matched distribution.");
      }

      if (siggen)
	delete siggen;
    }
    else
      throw OpalException("Distribution::CreateMatchedGaussDistribution", "didn't find any Cyclotron element in line");
  }
}

void Distribution::CreateDistributionGauss(size_t numberOfParticles, double massIneV) {

    SetDistParametersGauss(massIneV);

    if (emitting_m) {
        GenerateTransverseGauss(numberOfParticles);
        GenerateLongFlattopT(numberOfParticles);
    } else {
        GenerateGaussZ(numberOfParticles);
    }
}

void  Distribution::CreateBoundaryGeometry(PartBunch &beam, BoundaryGeometry &bg) {

    int pc = 0;
    size_t N_mean = static_cast<size_t>(floor(bg.getN() / Ippl::getNodes()));
    size_t N_extra = static_cast<size_t>(bg.getN() - N_mean * Ippl::getNodes());
    if(Ippl::myNode() == 0)
        N_mean += N_extra;
    size_t count = 0;
    size_t lowMark = beam.getLocalNum();
    if(bg.getN() != 0) {

        for(size_t i = 0; i < bg.getN(); i++) {
            if(pc == Ippl::myNode()) {
                if(count < N_mean) {
                    beam.create(1);
                    if(pc != 0)
                        beam.R[lowMark + count] = bg.getCooridinate(Ippl::myNode() * N_mean + count + N_extra);
                    else
                        beam.R[lowMark + count] = bg.getCooridinate(count);
                    beam.P[lowMark + count] = Vector_t(0.0);
                    beam.Bin[lowMark + count] = 0;
                    beam.PType[lowMark + count] = ParticleType::FIELDEMISSION;
                    beam.TriID[lowMark + count] = 0;
                    beam.Q[lowMark + count] = beam.getChargePerParticle();
                    beam.LastSection[lowMark + count] = 0;
                    beam.Ef[lowMark + count] = Vector_t(0.0);
                    beam.Bf[lowMark + count] = Vector_t(0.0);
                    beam.dt[lowMark + count] = beam.getdT();
                    count ++;
                }
            }
            pc++;
            if(pc == Ippl::getNodes())
                pc = 0;
        }
    }
    bg.clearCooridinateArray();
    beam.boundp();
    *gmsg << &beam << endl;
}

void Distribution::CreateOpalCycl(PartBunch &beam,
                                  size_t numberOfParticles,
                                  double current, const Beamline &bl,
                                  bool scan) {

    /*
     *  setup data for matched distribution generation
     */

    E_m = (beam.getGamma()-1.0)*beam.getM();
    I_m = current;
    bega_m = beam.getGamma()*beam.getBeta();

    /*
     * When scan mode is true, we need to destroy particles except
     * for the first pass.
     */

    /*
      Fixme:

      avrgpz_m = beam.getP()/beam.getM();
    */
    size_t numberOfPartToCreate = numberOfParticles;
    if (beam.getTotalNum() != 0) {
        scan_m = scan;
        numberOfPartToCreate = beam.getLocalNum();
    } else
        scan_m = false;

    /*
     * If in scan mode, destroy existing beam so we
     * can create new particles.
     */
    if (scan_m)
        DestroyBeam(beam);

    // Setup particle bin structure.
    SetupParticleBins(beam.getM(),beam);

    /*
     * Set what units to use for input momentum units. Default is
     * eV.
     */
    ChooseInputMomentumUnits(InputMomentumUnitsT::EV);

    // Set distribution type.
    SetDistType();

    // Emitting particles in not supported in OpalCyclT.
    emitting_m = false;

    /*
     * If Options::cZero is true than we reflect generated distribution
     * to ensure that the transverse averages are 0.0.
     *
     * For now we just cut the number of generated particles in half.
     */
    if (Options::cZero && !(distrTypeT_m == DistrTypeT::FROMFILE))
        numberOfPartToCreate /= 2;

    // Create distribution.
    Create(numberOfPartToCreate, beam.getM());

    // Now reflect particles if Options::cZero is true.
    if (Options::cZero && !(distrTypeT_m == DistrTypeT::FROMFILE))
        ReflectDistribution(numberOfPartToCreate);

    ShiftDistCoordinates(beam.getM());

    // Setup energy bins.
    if (numberOfEnergyBins_m > 0) {
        FillParticleBins();
        beam.setPBins(energyBins_m);
    }

    // Check number of particles in distribution.
    CheckParticleNumber(numberOfParticles);

    InitializeBeam(beam);
    WriteOutFileHeader();

    InjectBeam(beam);

}

void Distribution::CreateOpalE(Beam *beam,
                               std::vector<Distribution *> addedDistributions,
                               EnvelopeBunch *envelopeBunch,
                               double distCenter,
                               double Bz0) {

    IpplTimings::startTimer(envelopeBunch->distrCreate_m);

    double beamEnergy = beam->getMass() * (beam->getGamma() - 1.0) * 1.0e9;
    numberOfEnergyBins_m = static_cast<int>(fabs(Attributes::getReal(itsAttr[AttributesT::NBIN])));

    /*
     * Set what units to use for input momentum units. Default is
     * unitless (i.e. BetaXGamma, BetaYGamma, BetaZGamma).
     */
    ChooseInputMomentumUnits(InputMomentumUnitsT::NONE);

    SetDistType();

    // Check if this is to be an emitted beam.
    CheckIfEmitted();

    switch (distrTypeT_m) {

    case DistrTypeT::FLATTOP:
        SetDistParametersFlattop(beam->getMass());
        beamEnergy = Attributes::getReal(itsAttr[AttributesT::EKIN]);
        break;
    case DistrTypeT::GAUSS:
        SetDistParametersGauss(beam->getMass());
        break;
    case DistrTypeT::GUNGAUSSFLATTOPTH:
        INFOMSG("GUNGAUSSFLATTOPTH"<<endl);
        SetDistParametersFlattop(beam->getMass());
        beamEnergy = Attributes::getReal(itsAttr[AttributesT::EKIN]);
        break;
    default:
        *gmsg << "Only GAUSS and GUNGAUSSFLATTOPTH distribution types supported " << endl
              << "in envelope mode. Assuming FLATTOP." << endl;
        distrTypeT_m = DistrTypeT::FLATTOP;
        SetDistParametersFlattop(beam->getMass());
        break;
    }

    tEmission_m = tPulseLengthFWHM_m + (cutoffR_m[2] - sqrt(2.0 * log(2.0)))
        * (sigmaTRise_m + sigmaTFall_m);
    double beamWidth = tEmission_m * Physics::c * sqrt(1.0 - (1.0 / pow(beam->getGamma(), 2)));
    double beamCenter = -1.0 * beamWidth / 2.0;

    envelopeBunch->initialize(beam->getNumberOfSlices(),
                              beam->getCharge(),
                              beamEnergy,
                              beamWidth,
                              tEmission_m,
                              0.9,
                              beam->getCurrent(),
                              beamCenter,
                              sigmaR_m[0],
                              sigmaR_m[1],
                              0.0,
                              0.0,
                              Bz0,
                              numberOfEnergyBins_m);

    IpplTimings::stopTimer(envelopeBunch->distrCreate_m);
}

void Distribution::CreateOpalT(PartBunch &beam,
                               std::vector<Distribution *> addedDistributions,
                               size_t &numberOfParticles,
                               bool scan) {

    addedDistributions_m = addedDistributions;
    CreateOpalT(beam, numberOfParticles, scan);
}

void Distribution::CreateOpalT(PartBunch &beam,
                               size_t &numberOfParticles,
                               bool scan) {

    IpplTimings::startTimer(beam.distrCreate_m);

    // This is PC from BEAM
    avrgpz_m = beam.getP()/beam.getM();

    /*
     * If in scan mode, destroy existing beam so we
     * can create new particles.
     */
    scan_m = scan;
    if (scan_m)
        DestroyBeam(beam);

    /*
     * Set what units to use for input momentum units. Default is
     * unitless (i.e. BetaXGamma, BetaYGamma, BetaZGamma).
     */
    ChooseInputMomentumUnits(InputMomentumUnitsT::NONE);

    // Set distribution type(s).
    SetDistType();
    std::vector<Distribution *>::iterator addedDistIt;
    for (addedDistIt = addedDistributions_m.begin();
         addedDistIt != addedDistributions_m.end(); addedDistIt++)
        (*addedDistIt)->SetDistType();

    // Check if this is to be an emitted beam.
    CheckIfEmitted();

    if (emitting_m) {
        CheckEmissionParameters();
        SetupEmissionModel(beam);
    } else {
        pmean_m = Vector_t(0.0, 0.0, ConverteVToBetaGamma(GetEkin(), beam.getM()));
    }

    /*
     * If Options::cZero is true than we reflect generated distribution
     * to ensure that the transverse averages are 0.0.
     *
     * For now we just cut the number of generated particles in half.
     */
    if (Options::cZero && !(distrTypeT_m == DistrTypeT::FROMFILE))
        numberOfParticles /= 2;

    /*
     * Determine the number of particles for each distribution. Note
     * that if a distribution is generated from an input file, then
     * the number of particles in that file will override what is
     * determined here.
     */
    CalcPartPerDist(numberOfParticles);

    /*
     * Force added distributions to the same emission condition as the main
     * distribution.
     */
    for (addedDistIt = addedDistributions_m.begin();
         addedDistIt != addedDistributions_m.end(); addedDistIt++)
        (*addedDistIt)->SetDistToEmitted(emitting_m);

    // Create distributions.
    Create(particlesPerDist_m.at(0), beam.getM());

    size_t distCount = 1;
    for (addedDistIt = addedDistributions_m.begin();
         addedDistIt != addedDistributions_m.end(); addedDistIt++) {
        (*addedDistIt)->Create(particlesPerDist_m.at(distCount), beam.getM());
        distCount++;
    }

    // Move added distribution particles to main distribution.
    AddDistributions();

    // Now reflect particles if Options::cZero is true
    if (Options::cZero && !(distrTypeT_m == DistrTypeT::FROMFILE))
        ReflectDistribution(numberOfParticles);

    ShiftDistCoordinates(beam.getM());

    // Check number of particles in distribution.
    CheckParticleNumber(numberOfParticles);

    /*
     * Find max. and min. particle positions in the bunch. For the
     * case of an emitted beam these will be in time (seconds). For
     * an injected beam in z (meters).
     */
    double maxTOrZ = GetMaxTOrZ();
    double minTOrZ = GetMinTOrZ();

    /*
     * Set emission time and/or shift distributions if necessary.
     * For an emitted beam we place all particles at negative time.
     * For an injected beam we just ensure that there are no
     * particles at z < 0.
     */
    if (emitting_m)
        SetEmissionTime(maxTOrZ, minTOrZ);
    ShiftBeam(maxTOrZ, minTOrZ);

    if (numberOfEnergyBins_m > 0) {
        SetupEnergyBins(maxTOrZ, minTOrZ);
        FillEBinHistogram();
    } else
        energyBins_m = NULL;

    InitializeBeam(beam);
    WriteOutFileHeader();

    if (emitting_m && Options::cZero) {
        std::vector<std::vector<double> > mirrored;
        const auto end = additionalRNs_m.end();
        for (auto it = additionalRNs_m.begin(); it != end; ++ it) {
            std::vector<double> tmp((*it).begin() + 2, (*it).end());
            mirrored.push_back(tmp);
            (*it).erase((*it).begin() + 2, (*it).end());
        }
        additionalRNs_m.insert(additionalRNs_m.end(), mirrored.begin(), mirrored.end());
    }
    /*
     * If this is an injected beam, we create particles right away.
     * Emitted beams get created during the course of the simulation.
     */
    if (!emitting_m)
        InjectBeam(beam);

    IpplTimings::stopTimer(beam.distrCreate_m);
}

void Distribution::DestroyBeam(PartBunch &beam) {

    particlesPerDist_m.clear();
    xDist_m.clear();
    pxDist_m.clear();
    yDist_m.clear();
    pyDist_m.clear();
    tOrZDist_m.clear();
    pzDist_m.clear();
    xWrite_m.clear();
    pxWrite_m.clear();
    yWrite_m.clear();
    pyWrite_m.clear();
    tOrZWrite_m.clear();
    pzWrite_m.clear();
    binWrite_m.clear();

    beam.destroy(beam.getLocalNum(), 0);
    beam.update();
    INFOMSG("In scan mode: delete all particles in the bunch before next run.");

}

/**
 * Here we emit particles from the cathode.

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
size_t Distribution::EmitParticles(PartBunch &beam, double eZ) {

    // Number of particles that have already been emitted and are on this processor.
    size_t numberOfEmittedParticles = beam.getLocalNum();
    size_t oldNumberOfEmittedParticles = numberOfEmittedParticles;

    if (!tOrZDist_m.empty() && emitting_m) {

        /*
         * Loop through emission beam coordinate containers and emit particles with
         * the appropriate time coordinate. Once emitted, remove particle from the
         * "to be emitted" list.
         */
        std::vector<size_t> particlesToBeErased;
        for (size_t particleIndex = 0; particleIndex < tOrZDist_m.size(); particleIndex++) {

            // Advance particle time.
            tOrZDist_m.at(particleIndex) += beam.getdT();

            // If particle time is now greater than zero, we emit it.
            if (tOrZDist_m.at(particleIndex) >= 0.0) {

                particlesToBeErased.push_back(particleIndex);

                beam.create(1);
                double deltaT = tOrZDist_m.at(particleIndex);
                beam.dt[numberOfEmittedParticles] = deltaT;

                double oneOverCDt = 1.0 / (Physics::c * deltaT);

                double px = pxDist_m.at(particleIndex);
                double py = pyDist_m.at(particleIndex);
                double pz = pzDist_m.at(particleIndex);
                std::vector<double> additionalRNs;
                if (additionalRNs_m.size() > particleIndex) {
                    additionalRNs = additionalRNs_m[particleIndex];
                } else {
                    additionalRNs = std::vector<double>({gsl_rng_uniform(randGenEmit_m),
                                                         gsl_rng_uniform(randGenEmit_m)});
                }
                ApplyEmissionModel(eZ, px, py, pz, additionalRNs);

                double particleGamma
                    = std::sqrt(1.0
                                + std::pow(px, 2.0)
                                + std::pow(py, 2.0)
                                + std::pow(pz, 2.0));

                beam.R[numberOfEmittedParticles]
                    = Vector_t(oneOverCDt * (xDist_m.at(particleIndex)
                                             + px * deltaT * Physics::c / (2.0 * particleGamma)),
                               oneOverCDt * (yDist_m.at(particleIndex)
                                             + py * deltaT * Physics::c / (2.0 * particleGamma)),
                               pz / (2.0 * particleGamma));
                beam.P[numberOfEmittedParticles]
                    = Vector_t(px, py, pz);
                beam.Bin[numberOfEmittedParticles] = currentEnergyBin_m - 1;
                beam.Q[numberOfEmittedParticles] = beam.getChargePerParticle();
                beam.LastSection[numberOfEmittedParticles] = -1;
                beam.Ef[numberOfEmittedParticles] = Vector_t(0.0);
                beam.Bf[numberOfEmittedParticles] = Vector_t(0.0);
                beam.PType[numberOfEmittedParticles] = ParticleType::REGULAR;
                beam.TriID[numberOfEmittedParticles] = 0;
                numberOfEmittedParticles++;
                beam.iterateEmittedBin(currentEnergyBin_m - 1);

                // Save particles to vectors for writing initial distribution.
                xWrite_m.push_back(xDist_m.at(particleIndex));
                pxWrite_m.push_back(px);
                yWrite_m.push_back(yDist_m.at(particleIndex));
                pyWrite_m.push_back(py);
                tOrZWrite_m.push_back(-(beam.getdT() - deltaT + currentEmissionTime_m));
                pzWrite_m.push_back(pz);
                binWrite_m.push_back(currentEnergyBin_m);
            }
        }

        // Now erase particles that were emitted.
        std::vector<size_t>::reverse_iterator ptbErasedIt;
        for (ptbErasedIt = particlesToBeErased.rbegin();
             ptbErasedIt < particlesToBeErased.rend();
             ptbErasedIt++) {

            /*
             * We don't use the vector class function erase because it
             * is much slower than doing a swap and popping off the end
             * of the vector.
             */
            std::swap(xDist_m.at(*ptbErasedIt), xDist_m.back());
            std::swap(pxDist_m.at(*ptbErasedIt), pxDist_m.back());
            std::swap(yDist_m.at(*ptbErasedIt), yDist_m.back());
            std::swap(pyDist_m.at(*ptbErasedIt), pyDist_m.back());
            std::swap(tOrZDist_m.at(*ptbErasedIt), tOrZDist_m.back());
            std::swap(pzDist_m.at(*ptbErasedIt), pzDist_m.back());
            if (additionalRNs_m.size() == xDist_m.size()) {
                std::swap(additionalRNs_m.at(*ptbErasedIt), additionalRNs_m.back());

                additionalRNs_m.pop_back();
            }

            xDist_m.pop_back();
            pxDist_m.pop_back();
            yDist_m.pop_back();
            pyDist_m.pop_back();
            tOrZDist_m.pop_back();
            pzDist_m.pop_back();

        }

        /*
         * Set energy bin to emitted if all particles in the bin have been emitted.
         * However, be careful with the last energy bin. We cannot emit it until all
         * of the particles have been accounted for. So when on the last bin, keep it
         * open for the rest of the beam.
         */
        currentSampleBin_m++;
        if (currentSampleBin_m == numberOfSampleBins_m) {

            INFOMSG(level3 << "* Bin number: "
                    << currentEnergyBin_m
                    << " has emitted all particles (new emit)." << endl);
            currentSampleBin_m = 0;
            currentEnergyBin_m++;

        }

        /*
         * Set beam to emitted. Make sure temporary storage is cleared.
         */
        if (currentEnergyBin_m > numberOfEnergyBins_m || tOrZDist_m.empty()) {
            emitting_m = false;

            xDist_m.clear();
            pxDist_m.clear();
            yDist_m.clear();
            pyDist_m.clear();
            tOrZDist_m.clear();
            pzDist_m.clear();

            currentEnergyBin_m = numberOfEnergyBins_m;
        }

    }
    currentEmissionTime_m += beam.getdT();

    // Write emitted particles to file.
    WriteOutFileEmission();

    return numberOfEmittedParticles - oldNumberOfEmittedParticles;

}

void Distribution::EraseXDist() {
    xDist_m.erase(xDist_m.begin(), xDist_m.end());
}

void Distribution::EraseBGxDist() {
    pxDist_m.erase(pxDist_m.begin(), pxDist_m.end());
}

void Distribution::EraseYDist() {
    yDist_m.erase(yDist_m.begin(), yDist_m.end());
}

void Distribution::EraseBGyDist() {
    pyDist_m.erase(pyDist_m.begin(), pyDist_m.end());
}

void Distribution::EraseTOrZDist() {
    tOrZDist_m.erase(tOrZDist_m.begin(), tOrZDist_m.end());
}

void Distribution::EraseBGzDist() {
    pzDist_m.erase(pzDist_m.begin(), pzDist_m.end());
}

void Distribution::FillEBinHistogram() {

    /*
     * Loop over longitudinal beam coordinates and add them to the energy
     * bin histogram. Because of how we defined the bin structure the maximum
     * longitudinal coordinate will not be counted. So, we just force it in there
     * by artificially incrementing the last bin.
     */
    for (int processor = 0; processor < Ippl::getNodes(); processor++) {

        size_t numberOfParticles = 0;
        if (Ippl::myNode() == processor)
            numberOfParticles = tOrZDist_m.size();
        reduce(numberOfParticles, numberOfParticles, OpAddAssign());

        for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

            double tOrZ = 0.0;
            if (Ippl::myNode() == processor)
                tOrZ = tOrZDist_m.at(partIndex);
            reduce(tOrZ, tOrZ, OpAddAssign());

            if (gsl_histogram_increment(energyBinHist_m, tOrZ) == GSL_EDOM) {
                double upper = 0.0;
                double lower = 0.0;
                gsl_histogram_get_range(energyBinHist_m,
                                        gsl_histogram_bins(energyBinHist_m) - 1,
                                        &lower, &upper);
                gsl_histogram_increment(energyBinHist_m, lower);
            }

        }
    }
}

void Distribution::FillParticleBins() {

    for (size_t particleIndex = 0; particleIndex < tOrZDist_m.size(); particleIndex++) {

        std::vector<double> partCoord;

        partCoord.push_back(xDist_m.at(particleIndex));
        partCoord.push_back(yDist_m.at(particleIndex));
        partCoord.push_back(tOrZDist_m.at(particleIndex));
        partCoord.push_back(pxDist_m.at(particleIndex));
        partCoord.push_back(pyDist_m.at(particleIndex));
        partCoord.push_back(pzDist_m.at(particleIndex));
        partCoord.push_back(0.0);

        energyBins_m->fill(partCoord);

    }

    energyBins_m->sortArray();
}

size_t Distribution::FindEBin(double tOrZ) {

    if (tOrZ <= gsl_histogram_min(energyBinHist_m)) {
        return 0;
    } else if (tOrZ >= gsl_histogram_max(energyBinHist_m)) {
        return numberOfEnergyBins_m - 1;
    } else {
        gsl_histogram *energyBinHist = gsl_histogram_alloc(numberOfEnergyBins_m);
        gsl_histogram_set_ranges_uniform(energyBinHist,
                                         gsl_histogram_min(energyBinHist_m),
                                         gsl_histogram_max(energyBinHist_m));

        size_t binNumber;
        gsl_histogram_find(energyBinHist, tOrZ, &binNumber);
        gsl_histogram_free(energyBinHist);
        return binNumber;
    }
}

void Distribution::GenerateAstraFlattopT(size_t numberOfParticles) {

    /*
     * Legacy function to support the ASTRAFLATTOPOTH
     * distribution type.
     */
    CheckEmissionParameters();

    gsl_rng_env_setup();
    gsl_rng *randGen = gsl_rng_alloc(gsl_rng_1dhalton);
    gsl_qrng *quasiRandGen = gsl_qrng_alloc(gsl_qrng_halton, 2);

    int numberOfSampleBins
        = std::abs(static_cast<int> (Attributes::getReal(itsAttr[LegacyAttributesT::SBIN])));
    int numberOfEnergyBins
        = std::abs(static_cast<int> (Attributes::getReal(itsAttr[AttributesT::NBIN])));

    int binTotal = numberOfSampleBins * numberOfEnergyBins;

    double *distributionTable = new double[binTotal + 1];

    double a = tPulseLengthFWHM_m / 2.;
    double sig = tRise_m / 2.;
    double inv_erf08 = 0.906193802436823; // erfinv(0.8)
    double sqr2 = sqrt(2.0);
    double t = a - sqr2 * sig * inv_erf08;
    double tmps = sig;
    double tmpt = t;

    for(int i = 0; i < 10; ++ i) {
        sig = (t + tRise_m - a) / (sqr2 * inv_erf08);
        t = a - 0.5 * sqr2 * (sig + tmps) * inv_erf08;
        sig = (0.5 * (t + tmpt) + tRise_m - a) / (sqr2 * inv_erf08);
        tmps = sig;
        tmpt = t;
    }

    /*
     * Emission time is set here during distribution particle creation only for
     * the ASTRAFLATTOPTH distribution type.
     */
    tEmission_m = tPulseLengthFWHM_m + 10. * sig;
    tBin_m = tEmission_m / numberOfEnergyBins;

    double lo = -tBin_m / 2.0 * numberOfEnergyBins;
    double hi = tBin_m / 2.0 * numberOfEnergyBins;
    double dx = tBin_m / numberOfSampleBins;
    double x = lo;
    double tot = 0;
    double weight = 2.0;

    // sample the function that describes the histogram of the requested distribution
    for(int i = 0; i < binTotal + 1; ++ i, x += dx, weight = 6. - weight) {
        distributionTable[i] = gsl_sf_erf((x + a) / (sqrt(2) * sig)) - gsl_sf_erf((x - a) / (sqrt(2) * sig));
        tot += distributionTable[i] * weight;
    }
    tot -= distributionTable[binTotal] * (5. - weight);
    tot -= distributionTable[0];

    int saveProcessor = -1;
    double tCoord = 0.0;

    int effectiveNumParticles = 0;
    int largestBin = 0;
    std::vector<int> numParticlesInBin(numberOfEnergyBins,0);
    for(int k = 0; k < numberOfEnergyBins; ++ k) {
        double loc_fraction = -distributionTable[numberOfSampleBins * k] / tot;

        weight = 2.0;
        for(int i = numberOfSampleBins * k; i < numberOfSampleBins * (k + 1) + 1;
            ++ i, weight = 6. - weight) {
            loc_fraction += distributionTable[i] * weight / tot;
        }
        loc_fraction -= distributionTable[numberOfSampleBins * (k + 1)]
            * (5. - weight) / tot;
        numParticlesInBin[k] = static_cast<int>(std::floor(loc_fraction * numberOfParticles + 0.5));
        effectiveNumParticles += numParticlesInBin[k];
        if (numParticlesInBin[k] > numParticlesInBin[largestBin]) largestBin = k;
    }

    numParticlesInBin[largestBin] += (numberOfParticles - effectiveNumParticles);

    for(int k = 0; k < numberOfEnergyBins; ++ k) {
        gsl_ran_discrete_t *table
            = gsl_ran_discrete_preproc(numberOfSampleBins,
                                       &(distributionTable[numberOfSampleBins * k]));

        for(int i = 0; i < numParticlesInBin[k]; i++) {
            double xx[2];
            gsl_qrng_get(quasiRandGen, xx);
            tCoord = hi * (xx[1] + static_cast<int>(gsl_ran_discrete(randGen, table))
                           - binTotal / 2 + k * numberOfSampleBins) / (binTotal / 2);

            saveProcessor++;
            if (saveProcessor >= Ippl::getNodes())
                saveProcessor = 0;

            if (Ippl::myNode() == saveProcessor) {
                tOrZDist_m.push_back(tCoord);
                pzDist_m.push_back(0.0);
            }
        }
        gsl_ran_discrete_free(table);
    }

    gsl_rng_free(randGen);
    gsl_qrng_free(quasiRandGen);
    delete [] distributionTable;

}

void Distribution::GenerateBinomial(size_t numberOfParticles) {

    /*!
     *
     * \brief Following W. Johos for his report  <a href="http://gfa.web.psi.ch/publications/presentations/WernerJoho/TM-11-14.pdf"> TM-11-14 </a>
     *
     * For the \f$x,p_x\f$ phase space we have:
     * \f[
     *  \epsilon_x = \sigma_x \sigma_{p_x} \cos{( \arcsin{(\sigma_{12}) }) }
     *  \f]
     *
     * \f{eqnarray*}{
     * \beta_x &=& \frac{\sigma_x^2}{\epsilon_x} \\
     * \gamma_x &=& \frac{\sigma_{p_x}^2}{\epsilon_x} \\
     * \alpha_x &=& -\sigma_{12} \sqrt{(\beta_x \gamma_x)} \\
     * \f}
     *
     * This holds similar for the other dimensions.
     */


    // Calculate emittance and Twiss parameters.
    Vector_t emittance;
    for (unsigned int index = 0; index < 3; index++) {
        emittance(index) = sigmaR_m[index] * sigmaP_m[index]
            * cos(asin(correlationMatrix_m(2 * index + 1, 2 * index)));
    }

    Vector_t beta;
    Vector_t gamma;
    Vector_t alpha;
    for (unsigned int index = 0; index < 3; index++) {
        beta(index) = pow(sigmaR_m[index], 2.0) / emittance(index);
        gamma(index) = pow(sigmaP_m[index], 2.0) / emittance(index);
        alpha(index) = -correlationMatrix_m(2 * index + 1, 2 * index)
                        * sqrt(beta(index) * gamma(index));
    }

    Vector_t M = Vector_t(0.0);
    Vector_t PM = Vector_t(0.0);
    Vector_t COSCHI = Vector_t(0.0);
    Vector_t SINCHI = Vector_t(0.0);
    Vector_t CHI = Vector_t(0.0);
    Vector_t AMI = Vector_t(0.0);
    Vector_t L = Vector_t(0.0);
    Vector_t PL = Vector_t(0.0);

    for(unsigned int index = 0; index < 3; index++) {
        gamma(index) *= cutoffR_m[index];
        beta(index)  *= cutoffP_m[index];
        M[index]      =  sqrt(emittance(index) * beta(index));
        PM[index]     =  sqrt(emittance(index) * gamma(index));
        COSCHI[index] =  sqrt(1.0 / (1.0 + pow(alpha(index), 2.0)));
        SINCHI[index] = -alpha(index) * COSCHI[index];
        CHI[index]    =  atan2(SINCHI[index], COSCHI[index]);
        AMI[index]    =  1.0 / mBinomial_m[index];
        L[index]      =  sqrt((mBinomial_m[index] + 1.0) / 2.0) * M[index];
        PL[index]     =  sqrt((mBinomial_m[index] + 1.0) / 2.0) * PM[index];
    }

    int saveProcessor = -1;
    Vector_t x = Vector_t(0.0);
    Vector_t p = Vector_t(0.0);
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double S1 = 0.0;
        double S2 = 0.0;
        double A = 0.0;
        double AL = 0.0;
        double U = 0.0;
        double V = 0.0;

        S1 = IpplRandom();
        S2 = IpplRandom();

        if (mBinomial_m[0] <= 10000) {

            A = sqrt(1.0 - pow(S1, AMI[0]));
            AL = Physics::two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);
            double Ucp = U;
            double Vcp = V;
            x[0] = L[0] * U;
            p[0] = PL[0] * (U * SINCHI[0] + V * COSCHI[0]);

            S1 = IpplRandom();
            S2 = IpplRandom();
            A = sqrt(1.0 - pow(S1, AMI[1]));
            AL = Physics::two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);
            x[1] = L[1] * U;
            p[1] = PL[1] * (U * SINCHI[1] + V * COSCHI[1]);

            S1 = IpplRandom();
            S2 = IpplRandom();
            A = sqrt(1.0 - pow(S1, AMI[2]));
            AL = Physics::two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);

            double tempa = 1.0 - std::pow(correlationMatrix_m(1, 0), 2.0);
            const double l32 = (correlationMatrix_m(4, 1) - correlationMatrix_m(1, 0) * correlationMatrix_m(4, 0)) / sqrt(std::abs(tempa)) * tempa / std::abs(tempa);
            double tempb = 1 - std::pow(correlationMatrix_m(4, 0), 2.0) - l32 * l32;
            const double l33 = sqrt(std::abs(tempb)) * tempb / std::abs(tempb);
            x[2] = Ucp * correlationMatrix_m(4, 0) + Vcp * l32 + U * l33;
            const double l42 = (correlationMatrix_m(5, 1) - correlationMatrix_m(1, 0) * correlationMatrix_m(5, 0)) / sqrt(std::abs(tempa)) * tempa / std::abs(tempa);
            const double l43 = (correlationMatrix_m(5, 4) - correlationMatrix_m(4, 0) * correlationMatrix_m(5, 0) - l42 * l32) / l33;
            double tempc = 1 - std::pow(correlationMatrix_m(5, 0), 2.0) - l42 * l42 - l43 * l43;
            const double l44 = sqrt(std::abs(tempc)) * tempc / std::abs(tempc);

            p[2] = Ucp * correlationMatrix_m(5, 0) + Vcp * l42 + U * l43 + V * l44;
            x[2]  *= L[2];
            p[2] *= PL[2];

        } else {

            A = sqrt(2.0) / 2.0 * sqrt(-log(S1));
            AL = Physics::two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);
            double Ucp = U;
            double Vcp = V;
            x[0] = M[0] * U;
            p[0] = PM[0] * (U * SINCHI[0] + V * COSCHI[0]);


            S1 = IpplRandom();
            S2 = IpplRandom();
            A = sqrt(2.0) / 2.0 * sqrt(-log(S1));
            AL = Physics::two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);
            x[1] = M[1] * U;
            p[1] = PM[1] * (U * SINCHI[1] + V * COSCHI[1]);

            S1 = IpplRandom();
            S2 = IpplRandom();
            A = sqrt(2.0) / 2.0 * sqrt(-log(S1));
            AL = Physics::two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);

            double tempa = 1.0 - std::pow(correlationMatrix_m(1, 0), 2.0);
            const double l32 = std::copysign(1.0, tempa) * (correlationMatrix_m(4, 1) - correlationMatrix_m(1, 0) * correlationMatrix_m(4, 0)) / sqrt(std::abs(tempa));
            double tempb = 1 - std::pow(correlationMatrix_m(4, 0), 2.0) - l32 * l32;
            const double l33 = std::copysign(1.0, tempb) * sqrt(std::abs(tempb));
            x[2] = Ucp * correlationMatrix_m(4, 0) + Vcp * l32 + U * l33;
            const double l42 = std::copysign(1.0, tempa) * (correlationMatrix_m(5, 1) - correlationMatrix_m(1, 0) * correlationMatrix_m(5, 0)) / sqrt(std::abs(tempa));
            const double l43 = (correlationMatrix_m(5, 4) - correlationMatrix_m(4, 0) * correlationMatrix_m(5, 0) - l42 * l32) / l33;
            double tempc = 1 - std::pow(correlationMatrix_m(5, 0), 2.0) - l42 * l42 - l43 * l43;
            const double l44 = std::copysign(1.0, tempc) * sqrt(std::abs(tempc));

            p[2] = Ucp * correlationMatrix_m(5, 0) + Vcp * l42 + U * l43 + V * l44;
            x[2]  *= M[2];
            p[2] *= PM[2];

        }

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= Ippl::getNodes())
            saveProcessor = 0;

        if (Ippl::myNode() == saveProcessor) {
            xDist_m.push_back(x[0]);
            pxDist_m.push_back(p[0]);
            yDist_m.push_back(x[1]);
            pyDist_m.push_back(p[1]);
            tOrZDist_m.push_back(x[2]);
            pzDist_m.push_back(avrgpz_m * (1 + p[2]));
        }

    }
}

void Distribution::GenerateFlattopLaserProfile(size_t numberOfParticles) {

    int saveProcessor = -1;
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;

        laserProfile_m->getXY(x, y);

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= Ippl::getNodes())
            saveProcessor = 0;

        if (Ippl::myNode() == saveProcessor) {
            xDist_m.push_back(x * sigmaR_m[0]);
            pxDist_m.push_back(px);
            yDist_m.push_back(y * sigmaR_m[1]);
            pyDist_m.push_back(py);
        }
    }

    if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH)
        GenerateAstraFlattopT(numberOfParticles);
    else
        GenerateLongFlattopT(numberOfParticles);

}

void Distribution::GenerateFlattopT(size_t numberOfParticles) {

    gsl_rng_env_setup();
    gsl_rng  *randGenStandard = gsl_rng_alloc(gsl_rng_default);
    gsl_qrng *quasiRandGen2D = NULL;

    if(Options::rngtype != std::string("RANDOM")) {
        INFOMSG("RNGTYPE= " << Options::rngtype << endl);
        if(Options::rngtype == std::string("HALTON")) {
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_halton, 2);
        } else if(Options::rngtype == std::string("SOBOL")) {
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_sobol, 2);
        } else if(Options::rngtype == std::string("NIEDERREITER")) {
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);
        } else {
            INFOMSG("RNGTYPE= " << Options::rngtype << " not known, using HALTON" << endl);
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_halton, 2);
        }
    }

    int saveProcessor = -1;
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;

        bool allow = false;
        while (!allow) {

            if (quasiRandGen2D != NULL) {
                double randNums[2];
                gsl_qrng_get(quasiRandGen2D, randNums);
                x = -1.0 + 2.0 * randNums[0];
                y = -1.0 + 2.0 * randNums[1];
            } else {
                x = -1.0 + 2.0 * gsl_rng_uniform(randGenStandard);
                y = -1.0 + 2.0 * gsl_rng_uniform(randGenStandard);
            }

            allow = (pow(x, 2.0) + pow(y, 2.0) <= 1.0);

        }
        x *= sigmaR_m[0];
        y *= sigmaR_m[1];

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= Ippl::getNodes())
            saveProcessor = 0;

        if (Ippl::myNode() == saveProcessor) {
            xDist_m.push_back(x);
            pxDist_m.push_back(px);
            yDist_m.push_back(y);
            pyDist_m.push_back(py);
        }

    }

    gsl_rng_free(randGenStandard);

    if (quasiRandGen2D != NULL)
        gsl_qrng_free(quasiRandGen2D);

    if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH)
        GenerateAstraFlattopT(numberOfParticles);
    else
        GenerateLongFlattopT(numberOfParticles);

}

void Distribution::GenerateFlattopZ(size_t numberOfParticles) {

    gsl_rng_env_setup();
    gsl_rng *randGenStandard = gsl_rng_alloc(gsl_rng_default);
    gsl_qrng *quasiRandGen1D = NULL;
    gsl_qrng *quasiRandGen2D = NULL;
    if(Options::rngtype != std::string("RANDOM")) {
        INFOMSG("RNGTYPE= " << Options::rngtype << endl);
        if(Options::rngtype == std::string("HALTON")) {
            quasiRandGen1D = gsl_qrng_alloc(gsl_qrng_halton, 1);
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_halton, 2);
        } else if(Options::rngtype == std::string("SOBOL")) {
            quasiRandGen1D = gsl_qrng_alloc(gsl_qrng_sobol, 1);
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_sobol, 2);
        } else if(Options::rngtype == std::string("NIEDERREITER")) {
            quasiRandGen1D = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 1);
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);
        } else {
            INFOMSG("RNGTYPE= " << Options::rngtype << " not known, using HALTON" << endl);
            quasiRandGen1D = gsl_qrng_alloc(gsl_qrng_halton, 1);
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_halton, 2);
        }
    }

    int saveProcessor = -1;
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;
        double z = 0.0;
        double pz = 0.0;

        bool allow = false;
        while (!allow) {

            if (quasiRandGen2D != NULL) {
                double randNums[2];
                gsl_qrng_get(quasiRandGen2D, randNums);
                x = -1.0 + 2.0 * randNums[0];
                y = -1.0 + 2.0 * randNums[1];
            } else {
                x = -1.0 + 2.0 * gsl_rng_uniform(randGenStandard);
                y = -1.0 + 2.0 * gsl_rng_uniform(randGenStandard);
            }

            allow = (pow(x, 2.0) + pow(y, 2.0) <= 1.0);

        }
        x *= sigmaR_m[0];
        y *= sigmaR_m[1];

        if (quasiRandGen1D != NULL)
            gsl_qrng_get(quasiRandGen1D, &z);
        else
            z = gsl_rng_uniform(randGenStandard);

        z = (z - 0.5) * sigmaR_m[2];

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= Ippl::getNodes())
            saveProcessor = 0;

        if (Ippl::myNode() == saveProcessor) {
            xDist_m.push_back(x);
            pxDist_m.push_back(px);
            yDist_m.push_back(y);
            pyDist_m.push_back(py);
            tOrZDist_m.push_back(z);
            pzDist_m.push_back(pz);
        }
    }

    gsl_rng_free(randGenStandard);

    if (quasiRandGen1D != NULL) {
        gsl_qrng_free(quasiRandGen1D);
        gsl_qrng_free(quasiRandGen2D);
    }
}

void Distribution::GenerateGaussZ(size_t numberOfParticles) {

    gsl_rng_env_setup();
    gsl_rng *randGen = gsl_rng_alloc(gsl_rng_default);

    gsl_matrix * corMat  = gsl_matrix_alloc (6, 6);
    gsl_vector * rx = gsl_vector_alloc(6);
    gsl_vector * ry = gsl_vector_alloc(6);

    for (unsigned int i = 0; i < 6; ++ i) {
        gsl_matrix_set(corMat, i, i, correlationMatrix_m(i, i));
        for (unsigned int j = 0; j < i; ++ j) {
            gsl_matrix_set(corMat, i, j, correlationMatrix_m(i, j));
            gsl_matrix_set(corMat, j, i, correlationMatrix_m(i, j));
        }
    }

#define DISTDBG1
#ifdef DISTDBG1
    *gmsg << "* m before gsl_linalg_cholesky_decomp" << endl;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (j==0)
                *gmsg << "* r(" << std::setprecision(1) << i << ", : ) = "
                      << std::setprecision(4) << std::setw(10) << gsl_matrix_get (corMat, i, j);
            else
                *gmsg << " & " << std::setprecision(4) << std::setw(10) << gsl_matrix_get (corMat, i, j);
        }
        *gmsg << " \\\\" << endl;
    }
#endif

    int errcode = gsl_linalg_cholesky_decomp(corMat);

    if (errcode == GSL_EDOM) {
        throw OpalException("Distribution::GenerateGaussZ",
                            "gsl_linalg_cholesky_decomp failed");
    }
    // so make the upper part zero.
    for (int i = 0; i < 5; ++ i) {
        for (int j = i+1; j < 6 ; ++ j) {
            gsl_matrix_set (corMat, i, j, 0.0);
        }
    }

#define DISTDBG2
#ifdef DISTDBG2
    *gmsg << "* m after gsl_linalg_cholesky_decomp" << endl;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (j==0)
                *gmsg << "* r(" << std::setprecision(1) << i << ", : ) = "
                      << std::setprecision(4) << std::setw(10) << gsl_matrix_get (corMat, i, j);
            else
                *gmsg << " & " << std::setprecision(4) << std::setw(10) << gsl_matrix_get (corMat, i, j);
        }
        *gmsg << " \\\\" << endl;
    }
#endif

    int saveProcessor = -1;
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {
        bool allow = false;

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;
        double z = 0.0;
        double pz = 0.0;

        while (!allow) {
            gsl_vector_set(rx, 0, gsl_ran_gaussian(randGen, 1.0));
            gsl_vector_set(rx, 1, gsl_ran_gaussian(randGen, 1.0));
            gsl_vector_set(rx, 2, gsl_ran_gaussian(randGen, 1.0));
            gsl_vector_set(rx, 3, gsl_ran_gaussian(randGen, 1.0));
            gsl_vector_set(rx, 4, gsl_ran_gaussian(randGen, 1.0));
            gsl_vector_set(rx, 5, gsl_ran_gaussian(randGen, 1.0));

            gsl_blas_dgemv(CblasNoTrans, 1.0, corMat, rx, 0.0, ry);

            x = gsl_vector_get(ry, 0);
            px = gsl_vector_get(ry, 1);
            y = gsl_vector_get(ry, 2);
            py = gsl_vector_get(ry, 3);
            z = gsl_vector_get(ry, 4);
            pz = gsl_vector_get(ry, 5);

            bool xAndYOk = false;
            if (cutoffR_m[0] < SMALLESTCUTOFF && cutoffR_m[1] < SMALLESTCUTOFF)
                xAndYOk = true;
            else if (cutoffR_m[0] < SMALLESTCUTOFF)
                xAndYOk = (std::abs(y) <= cutoffR_m[1]);
            else if (cutoffR_m[1] < SMALLESTCUTOFF)
                xAndYOk = (std::abs(x) <= cutoffR_m[0]);
            else
                xAndYOk = (pow(x / cutoffR_m[0], 2.0) + pow(y / cutoffR_m[1], 2.0) <= 1.0);

            bool pxAndPyOk = false;
            if (cutoffP_m[0] < SMALLESTCUTOFF && cutoffP_m[1] < SMALLESTCUTOFF)
                pxAndPyOk = true;
            else if (cutoffP_m[0] < SMALLESTCUTOFF)
                pxAndPyOk = (std::abs(py) <= cutoffP_m[1]);
            else if (cutoffP_m[1] < SMALLESTCUTOFF)
                pxAndPyOk = (std::abs(px) <= cutoffP_m[0]);
            else
                pxAndPyOk = (pow(px / cutoffP_m[0], 2.0) + pow(py / cutoffP_m[1], 2.0) <= 1.0);

            allow = (xAndYOk && pxAndPyOk && std::abs(z) < cutoffR_m[2] && std::abs(pz) < cutoffP_m[2]);
        }

        x *= sigmaR_m[0];
        y *= sigmaR_m[1];
        z *= sigmaR_m[2];
        px *= sigmaP_m[0];
        py *= sigmaP_m[1];
        pz *= sigmaP_m[2];

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= Ippl::getNodes())
            saveProcessor = 0;

        if (Ippl::myNode() == saveProcessor) {
            xDist_m.push_back(x);
            pxDist_m.push_back(px);
            yDist_m.push_back(y);
            pyDist_m.push_back(py);
            tOrZDist_m.push_back(z);
            pzDist_m.push_back(avrgpz_m + pz);
        }
    }
    gsl_rng_free(randGen);
    gsl_vector_free(rx);
    gsl_vector_free(ry);
    gsl_matrix_free(corMat);
}

void Distribution::GenerateLongFlattopT(size_t numberOfParticles) {

    double flattopTime = tPulseLengthFWHM_m
        - sqrt(2.0 * log(2.0)) * (sigmaTRise_m + sigmaTFall_m);

    if (flattopTime < 0.0)
        flattopTime = 0.0;

    double normalizedFlankArea = 0.5 * sqrt(2.0 * Physics::pi) * gsl_sf_erf(cutoffR_m[2] / sqrt(2.0));
    double distArea = flattopTime
        + (sigmaTRise_m + sigmaTFall_m) * normalizedFlankArea;

    // Find number of particles in rise, fall and flat top.
    size_t numRise = numberOfParticles * sigmaTRise_m * normalizedFlankArea / distArea;
    size_t numFall = numberOfParticles * sigmaTFall_m * normalizedFlankArea / distArea;
    size_t numFlat = numberOfParticles - numRise - numFall;

    // Generate particles in tail.
    gsl_rng_env_setup();
    gsl_rng *randGenGSL = gsl_rng_alloc(gsl_rng_default);
    int saveProcessor = -1;

    for (size_t partIndex = 0; partIndex < numFall; partIndex++) {

        double t = 0.0;
        double pz = 0.0;

        bool allow = false;
        while (!allow) {
            t = gsl_ran_gaussian_tail(randGenGSL, 0, sigmaTFall_m);
            if (t <= sigmaTRise_m * cutoffR_m[2]) {
                t = -t + sigmaTFall_m * cutoffR_m[2];
                allow = true;
            }
        }

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= Ippl::getNodes())
            saveProcessor = 0;

        if (Ippl::myNode() == saveProcessor) {
            tOrZDist_m.push_back(t);
            pzDist_m.push_back(pz);
        }
    }

    /*
     * Generate particles in flat top. The flat top can also have sinusoidal
     * modulations.
     */
    double modulationAmp = Attributes::getReal(itsAttr[AttributesT::FTOSCAMPLITUDE])
        / 100.0;
    if (modulationAmp > 1.0)
        modulationAmp = 1.0;
    double numModulationPeriods
        = std::abs(Attributes::getReal(itsAttr[AttributesT::FTOSCPERIODS]));
    double modulationPeriod = 0.0;
    if (numModulationPeriods != 0.0)
        modulationPeriod = flattopTime / numModulationPeriods;

    gsl_rng_env_setup();
    gsl_rng *randGen = gsl_rng_alloc(gsl_rng_default);
    gsl_qrng *quasiRandGen1D = NULL;
    gsl_qrng *quasiRandGen2D = NULL;
    if(Options::rngtype != std::string("RANDOM")) {
        INFOMSG("RNGTYPE= " << Options::rngtype << endl);
        if(Options::rngtype == std::string("HALTON")) {
            quasiRandGen1D = gsl_qrng_alloc(gsl_qrng_halton, 1);
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_halton, 2);
        } else if(Options::rngtype == std::string("SOBOL")) {
            quasiRandGen1D = gsl_qrng_alloc(gsl_qrng_sobol, 1);
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_sobol, 2);
        } else if(Options::rngtype == std::string("NIEDERREITER")) {
            quasiRandGen1D = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 1);
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);
        } else {
            INFOMSG("RNGTYPE= " << Options::rngtype << " not known, using HALTON" << endl);
            quasiRandGen1D = gsl_qrng_alloc(gsl_qrng_halton, 1);
            quasiRandGen2D = gsl_qrng_alloc(gsl_qrng_halton, 2);
        }
    }

    for (size_t partIndex = 0; partIndex < numFlat; partIndex++) {

        double t = 0.0;
        double pz = 0.0;

        if (modulationAmp == 0.0 || numModulationPeriods == 0.0) {

            if (quasiRandGen1D != NULL)
                gsl_qrng_get(quasiRandGen1D, &t);
            else
                t = gsl_rng_uniform(randGen);

            t = flattopTime * t + sigmaTFall_m * cutoffR_m[2];

        } else {

            double funcAmp = 0.0;
            bool allow = false;
            while (!allow) {
                if (quasiRandGen2D != NULL) {
                    double randNums[2] = {0.0, 0.0};
                    gsl_qrng_get(quasiRandGen2D, randNums);
                    t = randNums[0] * flattopTime;
                    funcAmp = randNums[1];
                } else {
                    t = gsl_rng_uniform(randGen)*flattopTime;
                    funcAmp = gsl_rng_uniform(randGen);
                }

                double funcValue = (1.0 + modulationAmp
                                    * sin(Physics::two_pi * t / modulationPeriod))
                    / (1.0 + std::abs(modulationAmp));

                allow = (funcAmp <= funcValue);

            }
        }

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= Ippl::getNodes())
            saveProcessor = 0;

        if (Ippl::myNode() == saveProcessor) {
            tOrZDist_m.push_back(t);
            pzDist_m.push_back(pz);
        }
    }

    // Generate particles in rise.
    for (size_t partIndex = 0; partIndex < numRise; partIndex++) {

        double t = 0.0;
        double pz = 0.0;

        bool allow = false;
        while (!allow) {
            t = gsl_ran_gaussian_tail(randGenGSL, 0, sigmaTRise_m);
            if (t <= sigmaTRise_m * cutoffR_m[2]) {
                t += sigmaTFall_m * cutoffR_m[2] + flattopTime;
                allow = true;
            }
        }

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= Ippl::getNodes())
            saveProcessor = 0;

        if (Ippl::myNode() == saveProcessor) {
            tOrZDist_m.push_back(t);
            pzDist_m.push_back(pz);
        }
    }

    gsl_rng_free(randGenGSL);
    gsl_rng_free(randGen);

    if (quasiRandGen1D != NULL) {
        gsl_qrng_free(quasiRandGen1D);
        gsl_qrng_free(quasiRandGen2D);
    }
}

void Distribution::GenerateTransverseGauss(size_t numberOfParticles) {

    // Generate coordinates.
    gsl_rng_env_setup();
    gsl_rng *randGen = gsl_rng_alloc(gsl_rng_default);
    gsl_matrix * corMat  = gsl_matrix_alloc (4, 4);
    gsl_vector * rx = gsl_vector_alloc(4);
    gsl_vector * ry = gsl_vector_alloc(4);

    for (unsigned int i = 0; i < 4; ++ i) {
        gsl_matrix_set(corMat, i, i, correlationMatrix_m(i, i));
        for (unsigned int j = 0; j < i; ++ j) {
            gsl_matrix_set(corMat, i, j, correlationMatrix_m(i, j));
            gsl_matrix_set(corMat, j, i, correlationMatrix_m(i, j));
        }
    }

    int errcode = gsl_linalg_cholesky_decomp(corMat);

    if (errcode == GSL_EDOM) {
        throw OpalException("Distribution::GenerateTransverseGauss",
                            "gsl_linalg_cholesky_decomp failed");
    }

    for (int i = 0; i < 3; ++ i) {
        for (int j = i+1; j < 4 ; ++ j) {
            gsl_matrix_set (corMat, i, j, 0.0);
        }
    }

    int saveProcessor = -1;
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        double x = 0.0;
        double px = 0.0;
        double y = 0.0;
        double py = 0.0;

        bool allow = false;
        while (!allow) {

            gsl_vector_set(rx, 0, gsl_ran_gaussian (randGen,1.0));
            gsl_vector_set(rx, 1, gsl_ran_gaussian (randGen,1.0));
            gsl_vector_set(rx, 2, gsl_ran_gaussian (randGen,1.0));
            gsl_vector_set(rx, 3, gsl_ran_gaussian (randGen,1.0));

            gsl_blas_dgemv(CblasNoTrans, 1.0, corMat, rx, 0.0, ry);
            x = gsl_vector_get(ry, 0);
            px = gsl_vector_get(ry, 1);
            y = gsl_vector_get(ry, 2);
            py = gsl_vector_get(ry, 3);

            bool xAndYOk = false;
            if (cutoffR_m[0] < SMALLESTCUTOFF && cutoffR_m[1] < SMALLESTCUTOFF)
                xAndYOk = true;
            else if (cutoffR_m[0] < SMALLESTCUTOFF)
                xAndYOk = (std::abs(y) <= cutoffR_m[1]);
            else if (cutoffR_m[1] < SMALLESTCUTOFF)
                xAndYOk = (std::abs(x) <= cutoffR_m[0]);
            else
                xAndYOk = (pow(x / cutoffR_m[0], 2.0) + pow(y / cutoffR_m[1], 2.0) <= 1.0);

            bool pxAndPyOk = false;
            if (cutoffP_m[0] < SMALLESTCUTOFF && cutoffP_m[1] < SMALLESTCUTOFF)
                pxAndPyOk = true;
            else if (cutoffP_m[0] < SMALLESTCUTOFF)
                pxAndPyOk = (std::abs(py) <= cutoffP_m[1]);
            else if (cutoffP_m[1] < SMALLESTCUTOFF)
                pxAndPyOk = (std::abs(px) <= cutoffP_m[0]);
            else
                pxAndPyOk = (pow(px / cutoffP_m[0], 2.0) + pow(py / cutoffP_m[1], 2.0) <= 1.0);

            allow = (xAndYOk && pxAndPyOk);

        }
        x *= sigmaR_m[0];
        y *= sigmaR_m[1];
        px *= sigmaP_m[0];
        py *= sigmaP_m[1];

        // Save to each processor in turn.
        saveProcessor++;
        if (saveProcessor >= Ippl::getNodes())
            saveProcessor = 0;

        if (Ippl::myNode() == saveProcessor) {
            xDist_m.push_back(x);
            pxDist_m.push_back(px);
            yDist_m.push_back(y);
            pyDist_m.push_back(py);
        }
    }

    gsl_rng_free(randGen);
    gsl_matrix_free(corMat);
    gsl_vector_free(rx);
    gsl_vector_free(ry);
}

void Distribution::InitializeBeam(PartBunch &beam) {

    /*
     * Set emission time, the current beam bunch number and
     * set the beam energy bin structure, if it has one.
     */
    beam.setTEmission(tEmission_m);
    beam.setNumBunch(1);
    if (numberOfEnergyBins_m > 0)
        beam.SetEnergyBins(numberOfEnergyBins_m);
    beam.LastSection = 0;
}

void Distribution::InjectBeam(PartBunch &beam) {

    WriteOutFileInjection();

    std::vector<double> id1 = Attributes::getRealArray(itsAttr[AttributesT::ID1]);
    std::vector<double> id2 = Attributes::getRealArray(itsAttr[AttributesT::ID2]);

    bool hasID1 = (id1.size() != 0);
    bool hasID2 = (id2.size() != 0);

    if (hasID1 || hasID2)
        *gmsg << "* Use special ID1 or ID2 particle in distribution" << endl;


    size_t numberOfParticles = tOrZDist_m.size();
    beam.create(numberOfParticles);
    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

        beam.R[partIndex] = Vector_t(xDist_m.at(partIndex),
                                     yDist_m.at(partIndex),
                                     tOrZDist_m.at(partIndex));

        beam.P[partIndex] = Vector_t(pxDist_m.at(partIndex),
                                     pyDist_m.at(partIndex),
                                     pzDist_m.at(partIndex));

        beam.Q[partIndex] = beam.getChargePerParticle();
        beam.Ef[partIndex] = Vector_t(0.0);
        beam.Bf[partIndex] = Vector_t(0.0);
        beam.PType[partIndex] = ParticleType::REGULAR;
        beam.TriID[partIndex] = 0;

        if (numberOfEnergyBins_m > 0) {
            size_t binNumber = FindEBin(tOrZDist_m.at(partIndex));
            beam.Bin[partIndex] = binNumber;
            beam.iterateEmittedBin(binNumber);
        } else
            beam.Bin[partIndex] = 0;

        if (hasID1) {
            if (beam.ID[partIndex] == 1) {
                beam.R[partIndex] = Vector_t(id1[0],id1[1],id1[2]);
                beam.P[partIndex] = Vector_t(id1[3],id1[4],id1[5]);
            }
        }

        if (hasID2) {
            if (beam.ID[partIndex] == 2) {
                beam.R[partIndex] = Vector_t(id2[0],id2[1],id2[2]);
                beam.P[partIndex] = Vector_t(id2[3],id2[4],id2[5]);
            }
        }
    }

    xDist_m.clear();
    pxDist_m.clear();
    yDist_m.clear();
    pyDist_m.clear();
    tOrZDist_m.clear();
    pzDist_m.clear();

    beam.boundp();

    if (distrTypeT_m != DistrTypeT::FROMFILE)
        beam.correctEnergy(avrgpz_m);
    beam.calcEMean();
}

bool Distribution::GetIfDistEmitting() {
    return emitting_m;
}

int Distribution::GetLastEmittedEnergyBin() {
    return currentEnergyBin_m;
}

double Distribution::GetMaxTOrZ() {

    double maxTOrZ = 0.0;
    std::vector<double>::iterator longIt;
    for (longIt = tOrZDist_m.begin(); longIt != tOrZDist_m.end(); longIt++) {
        if (longIt == tOrZDist_m.begin())
            maxTOrZ = *longIt;
        else
            if (maxTOrZ < *longIt)
                maxTOrZ = *longIt;
    }

    std::vector<double> maxTOrZs;
    for (int nodeIndex = 0; nodeIndex < Ippl::getNodes(); nodeIndex++) {
        if (nodeIndex == Ippl::myNode())
            maxTOrZs.push_back(maxTOrZ);
        else
            maxTOrZs.push_back(0.0);

        reduce(maxTOrZs.at(nodeIndex), maxTOrZs.at(nodeIndex), OpAddAssign());
    }

    std::vector<double>::iterator maxIt;
    for (maxIt = maxTOrZs.begin(); maxIt != maxTOrZs.end(); maxIt++) {
        if (maxIt == maxTOrZs.begin())
            maxTOrZ = *maxIt;
        else if (maxTOrZ < *maxIt)
            maxTOrZ = *maxIt;
    }

    return maxTOrZ;
}

double Distribution::GetMinTOrZ() {

    double minTOrZ = 0.0;
    std::vector<double>::iterator longIt;
    for (longIt = tOrZDist_m.begin(); longIt != tOrZDist_m.end(); longIt++) {
        if (longIt == tOrZDist_m.begin())
            minTOrZ = *longIt;
        else
            if (minTOrZ > *longIt)
                minTOrZ = *longIt;
    }

    std::vector<double> minTOrZs;
    for (int nodeIndex = 0; nodeIndex < Ippl::getNodes(); nodeIndex++) {
        if (nodeIndex == Ippl::myNode())
            minTOrZs.push_back(minTOrZ);
        else
            minTOrZs.push_back(0.0);

        reduce(minTOrZs.at(nodeIndex), minTOrZs.at(nodeIndex), OpAddAssign());
    }

    std::vector<double>::iterator minIt;
    for (minIt = minTOrZs.begin(); minIt != minTOrZs.end(); minIt++) {
        if (minIt == minTOrZs.begin())
            minTOrZ = *minIt;
        else if (minTOrZ > *minIt)
            minTOrZ = *minIt;
    }

    return minTOrZ;
}

size_t Distribution::GetNumberOfEmissionSteps() {
    return static_cast<size_t> (numberOfEnergyBins_m * numberOfSampleBins_m);
}

int Distribution::GetNumberOfEnergyBins() {
    return numberOfEnergyBins_m;
}

double Distribution::GetEmissionDeltaT() {
    return tBin_m / numberOfSampleBins_m;
}

double Distribution::GetEnergyBinDeltaT() {
    return tBin_m;
}

double Distribution::GetWeight() {
    return std::abs(Attributes::getReal(itsAttr[AttributesT::WEIGHT]));
}

std::vector<double>& Distribution::GetXDist() {
    return xDist_m;
}

std::vector<double>& Distribution::GetBGxDist() {
    return pxDist_m;
}

std::vector<double>& Distribution::GetYDist() {
    return yDist_m;
}

std::vector<double>& Distribution::GetBGyDist() {
    return pyDist_m;
}

std::vector<double>& Distribution::GetTOrZDist() {
    return tOrZDist_m;
}

std::vector<double>& Distribution::GetBGzDist() {
    return pzDist_m;
}

void Distribution::PrintDist(Inform &os, size_t numberOfParticles) const {

    if (numberOfParticles > 0) {
        os << "* Number of particles: "
           << numberOfParticles * (Options::cZero && !(distrTypeT_m == DistrTypeT::FROMFILE)? 2: 1)
           << endl
           << "* " << endl;
    }

    switch (distrTypeT_m) {

    case DistrTypeT::FROMFILE:
        PrintDistFromFile(os);
        break;
    case DistrTypeT::GAUSS:
        PrintDistGauss(os);
        break;
    case DistrTypeT::BINOMIAL:
        PrintDistBinomial(os);
        break;
    case DistrTypeT::FLATTOP:
        PrintDistFlattop(os);
        break;
    case DistrTypeT::SURFACEEMISSION:
        PrintDistSurfEmission(os);
        break;
    case DistrTypeT::SURFACERANDCREATE:
        PrintDistSurfAndCreate(os);
        break;
    case DistrTypeT::GUNGAUSSFLATTOPTH:
        PrintDistFlattop(os);
        break;
    case DistrTypeT::ASTRAFLATTOPTH:
        PrintDistFlattop(os);
        break;
    case DistrTypeT::MATCHEDGAUSS:
        PrintDistMatchedGauss(os);
        break;
    default:
        INFOMSG("Distribution unknown." << endl;);
        break;
    }

}

void Distribution::PrintDistBinomial(Inform &os) const {

    os << "* Distribution type: BINOMIAL" << endl;
    os << "* " << endl;
    os << "* SIGMAX   = " << sigmaR_m[0] << " [m]" << endl;
    os << "* SIGMAY   = " << sigmaR_m[1] << " [m]" << endl;
    if (emitting_m)
        os << "* SIGMAT   = " << sigmaR_m[2] << " [sec]" << endl;
    else
        os << "* SIGMAZ   = " << sigmaR_m[2] << " [m]" << endl;
    os << "* SIGMAPX  = " << sigmaP_m[0] << " [Beta Gamma]" << endl;
    os << "* SIGMAPY  = " << sigmaP_m[1] << " [Beta Gamma]" << endl;
    os << "* SIGMAPZ  = " << sigmaP_m[2] << " [Beta Gamma]" << endl;
    os << "* MX       = " << mBinomial_m[0] << endl;
    os << "* MY       = " << mBinomial_m[1] << endl;
    if (emitting_m)
        os << "* MT       = " << mBinomial_m[2] << endl;
    else
        os << "* MZ       = " << mBinomial_m[2] << endl;
    os << "* CORRX    = " << correlationMatrix_m(1, 0) << endl;
    os << "* CORRY    = " << correlationMatrix_m(3, 2) << endl;
    os << "* CORRZ    = " << correlationMatrix_m(5, 4) << endl;
    os << "* R61      = " << correlationMatrix_m(5, 0) << endl;
    os << "* R62      = " << correlationMatrix_m(5, 1) << endl;
    os << "* R51      = " << correlationMatrix_m(4, 0) << endl;
    os << "* R52      = " << correlationMatrix_m(4, 1) << endl;
}

void Distribution::PrintDistFlattop(Inform &os) const {

    switch (distrTypeT_m) {

    case DistrTypeT::ASTRAFLATTOPTH:
        os << "* Distribution type: ASTRAFLATTOPTH" << endl;
        break;

    case DistrTypeT::GUNGAUSSFLATTOPTH:
        os << "* Distribution type: GUNGAUSSFLATTOPTH" << endl;
        break;

    default:
        os << "* Distribution type: FLATTOP" << endl;
        break;

    }
    os << "* " << endl;

    if (laserProfile_m != NULL) {

        os << "* Transverse profile determined by laser image: " << endl
           << endl
           << "* Laser profile: " << laserProfileFileName_m << endl
           << "* Laser image:   " << laserImageName_m << endl
           << "* Intensity cut: " << laserIntensityCut_m << endl;

    } else {

        os << "* SIGMAX                        = " << sigmaR_m[0] << " [m]" << endl;
        os << "* SIGMAY                        = " << sigmaR_m[1] << " [m]" << endl;

    }

    if (emitting_m) {

        if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH) {

            os << "* Time Rise                     = " << tRise_m
               << " [sec]" << endl;
            os << "* TPULSEFWHM                    = " << tPulseLengthFWHM_m
               << " [sec]" << endl;

        } else {
            os << "* Sigma Time Rise               = " << sigmaTRise_m
               << " [sec]" << endl;
            os << "* TPULSEFWHM                    = " << tPulseLengthFWHM_m
               << " [sec]" << endl;
            os << "* Sigma Time Fall               = " << sigmaTFall_m
               << " [sec]" << endl;
            os << "* Longitudinal cutoff           = " << cutoffR_m[2]
               << " [units of Sigma Time]" << endl;
            os << "* Flat top modulation amplitude = "
               << Attributes::getReal(itsAttr[AttributesT::FTOSCAMPLITUDE])
               << " [Percent of distribution amplitude]" << endl;
            os << "* Flat top modulation periods   = "
               << std::abs(Attributes::getReal(itsAttr[AttributesT::FTOSCPERIODS]))
               << endl;
        }

    } else
        os << "* SIGMAZ                        = " << sigmaR_m[2]
           << " [m]" << endl;
}

void Distribution::PrintDistFromFile(Inform &os) const {
    os << "* Distribution type: FROMFILE" << endl;
    os << "* " << endl;
    os << "* Input file:        "
       << Attributes::getString(itsAttr[AttributesT::FNAME]) << endl;
}


void Distribution::PrintDistMatchedGauss(Inform &os) const {
    os << "* Distribution type: MATCHEDGAUSS" << endl;
    os << "* SIGMAX   = " << sigmaR_m[0] << " [m]" << endl;
    os << "* SIGMAY   = " << sigmaR_m[1] << " [m]" << endl;
    os << "* SIGMAZ   = " << sigmaR_m[2] << " [m]" << endl;
    os << "* SIGMAPX  = " << sigmaP_m[0] << " [Beta Gamma]" << endl;
    os << "* SIGMAPY  = " << sigmaP_m[1] << " [Beta Gamma]" << endl;
    os << "* SIGMAPZ  = " << sigmaP_m[2] << " [Beta Gamma]" << endl;
    os << "* AVRGPZ   = " << avrgpz_m <<    " [Beta Gamma]" << endl;

    os << "* CORRX    = " << correlationMatrix_m(1, 0) << endl;
    os << "* CORRY    = " << correlationMatrix_m(3, 2) << endl;
    os << "* CORRZ    = " << correlationMatrix_m(5, 4) << endl;
    os << "* R61      = " << correlationMatrix_m(5, 0) << endl;
    os << "* R62      = " << correlationMatrix_m(5, 1) << endl;
    os << "* R51      = " << correlationMatrix_m(4, 0) << endl;
    os << "* R52      = " << correlationMatrix_m(4, 1) << endl;
    os << "* CUTOFFX  = " << cutoffR_m[0] << " [units of SIGMAX]" << endl;
    os << "* CUTOFFY  = " << cutoffR_m[1] << " [units of SIGMAY]" << endl;
    os << "* CUTOFFZ  = " << cutoffR_m[2] << " [units of SIGMAZ]" << endl;
    os << "* CUTOFFPX = " << cutoffP_m[0] << " [units of SIGMAPX]" << endl;
    os << "* CUTOFFPY = " << cutoffP_m[1] << " [units of SIGMAPY]" << endl;
    os << "* CUTOFFPZ = " << cutoffP_m[2] << " [units of SIGMAPY]" << endl;
}

void Distribution::PrintDistGauss(Inform &os) const {
    os << "* Distribution type: GAUSS" << endl;
    os << "* " << endl;
    if (emitting_m) {
        os << "* Sigma Time Rise               = " << sigmaTRise_m
           << " [sec]" << endl;
        os << "* TPULSEFWHM                    = " << tPulseLengthFWHM_m
           << " [sec]" << endl;
        os << "* Sigma Time Fall               = " << sigmaTFall_m
           << " [sec]" << endl;
        os << "* Longitudinal cutoff           = " << cutoffR_m[2]
           << " [units of Sigma Time]" << endl;
        os << "* Flat top modulation amplitude = "
           << Attributes::getReal(itsAttr[AttributesT::FTOSCAMPLITUDE])
           << " [Percent of distribution amplitude]" << endl;
        os << "* Flat top modulation periods   = "
           << std::abs(Attributes::getReal(itsAttr[AttributesT::FTOSCPERIODS]))
           << endl;
        os << "* SIGMAX                        = " << sigmaR_m[0] << " [m]" << endl;
        os << "* SIGMAY                        = " << sigmaR_m[1] << " [m]" << endl;
        os << "* SIGMAPX                       = " << sigmaP_m[0]
           << " [Beta Gamma]" << endl;
        os << "* SIGMAPY                       = " << sigmaP_m[1]
           << " [Beta Gamma]" << endl;
        os << "* CORRX                         = " << correlationMatrix_m(1, 0) << endl;
        os << "* CORRY                         = " << correlationMatrix_m(3, 2) << endl;
        os << "* CUTOFFX                       = " << cutoffR_m[0]
           << " [units of SIGMAX]" << endl;
        os << "* CUTOFFY                       = " << cutoffR_m[1]
           << " [units of SIGMAY]" << endl;
        os << "* CUTOFFPX                      = " << cutoffP_m[0]
           << " [units of SIGMAPX]" << endl;
        os << "* CUTOFFPY                      = " << cutoffP_m[1]
           << " [units of SIGMAPY]" << endl;
    } else {
        os << "* SIGMAX   = " << sigmaR_m[0] << " [m]" << endl;
        os << "* SIGMAY   = " << sigmaR_m[1] << " [m]" << endl;
        os << "* SIGMAZ   = " << sigmaR_m[2] << " [m]" << endl;
        os << "* SIGMAPX  = " << sigmaP_m[0] << " [Beta Gamma]" << endl;
        os << "* SIGMAPY  = " << sigmaP_m[1] << " [Beta Gamma]" << endl;
        os << "* SIGMAPZ  = " << sigmaP_m[2] << " [Beta Gamma]" << endl;
        os << "* AVRGPZ   = " << avrgpz_m <<    " [Beta Gamma]" << endl;

        os << "* CORRX    = " << correlationMatrix_m(1, 0) << endl;
        os << "* CORRY    = " << correlationMatrix_m(3, 2) << endl;
        os << "* CORRZ    = " << correlationMatrix_m(5, 4) << endl;
        os << "* R61      = " << correlationMatrix_m(5, 0) << endl;
        os << "* R62      = " << correlationMatrix_m(5, 1) << endl;
        os << "* R51      = " << correlationMatrix_m(4, 0) << endl;
        os << "* R52      = " << correlationMatrix_m(4, 1) << endl;
        os << "* CUTOFFX  = " << cutoffR_m[0] << " [units of SIGMAX]" << endl;
        os << "* CUTOFFY  = " << cutoffR_m[1] << " [units of SIGMAY]" << endl;
        os << "* CUTOFFZ  = " << cutoffR_m[2] << " [units of SIGMAZ]" << endl;
        os << "* CUTOFFPX = " << cutoffP_m[0] << " [units of SIGMAPX]" << endl;
        os << "* CUTOFFPY = " << cutoffP_m[1] << " [units of SIGMAPY]" << endl;
        os << "* CUTOFFPZ = " << cutoffP_m[2] << " [units of SIGMAPY]" << endl;
    }
}

void Distribution::PrintDistSurfEmission(Inform &os) const {

    os << "* Distribution type: SURFACEEMISION" << endl;
    os << "* " << endl;
    os << "* * Number of electrons for surface emission  "
       << Attributes::getReal(itsAttr[AttributesT::NPDARKCUR]) << endl;
    os << "* * Initialized electrons inward margin for surface emission  "
       << Attributes::getReal(itsAttr[AttributesT::INWARDMARGIN]) << endl;
    os << "* * E field threshold (MV), only in position r with E(r)>EINITHR that "
       << "particles will be initialized   "
       << Attributes::getReal(itsAttr[AttributesT::EINITHR]) << endl;
    os << "* * Field enhancement for surface emission  "
       << Attributes::getReal(itsAttr[AttributesT::FNBETA]) << endl;
    os << "* * Maximum number of electrons emitted from a single triangle for "
       << "Fowler-Nordheim emission  "
       << Attributes::getReal(itsAttr[AttributesT::FNMAXEMI]) << endl;
    os << "* * Field Threshold for Fowler-Nordheim emission (MV/m)  "
       <<  Attributes::getReal(itsAttr[AttributesT::FNFIELDTHR]) << endl;
    os << "* * Empirical constant A for Fowler-Nordheim emission model  "
       <<  Attributes::getReal(itsAttr[AttributesT::FNA]) << endl;
    os << "* * Empirical constant B for Fowler-Nordheim emission model  "
       <<  Attributes::getReal(itsAttr[AttributesT::FNB]) << endl;
    os << "* * Constant for image charge effect parameter y(E) in Fowler-Nordheim "
       << "emission model  "
       <<  Attributes::getReal(itsAttr[AttributesT::FNY]) << endl;
    os << "* * Zero order constant for image charge effect parameter v(y) in "
       << "Fowler-Nordheim emission model  "
       <<  Attributes::getReal(itsAttr[AttributesT::FNVYZERO]) << endl;
    os << "* * Second order constant for image charge effect parameter v(y) in "
       << "Fowler-Nordheim emission model  "
       <<  Attributes::getReal(itsAttr[AttributesT::FNVYSECOND]) << endl;
    os << "* * Select secondary model type(0:no secondary emission; 1:Furman-Pivi; "
       << "2 or larger: Vaughan's model  "
       <<  Attributes::getReal(itsAttr[AttributesT::SECONDARYFLAG]) << endl;
    os << "* * Secondary emission mode type(true: emit n true secondaries; false: "
       << "emit one particle with n times charge  "
       << Attributes::getBool(itsAttr[AttributesT::NEMISSIONMODE]) << endl;
    os << "* * Sey_0 in Vaughan's model "
       <<  Attributes::getReal(itsAttr[AttributesT::VSEYZERO]) << endl;
    os << "* * Energy related to sey_0 in Vaughan's model in eV  "
       <<  Attributes::getReal(itsAttr[AttributesT::VEZERO]) << endl;
    os << "* * Sey max in Vaughan's model  "
       <<  Attributes::getReal(itsAttr[AttributesT::VSEYMAX]) << endl;
    os << "* * Emax in Vaughan's model in eV  "
       <<  Attributes::getReal(itsAttr[AttributesT::VEMAX]) << endl;
    os << "* * Fitting parameter denotes the roughness of surface for impact "
       << "energy in Vaughan's model  "
       <<  Attributes::getReal(itsAttr[AttributesT::VKENERGY]) << endl;
    os << "* * Fitting parameter denotes the roughness of surface for impact angle "
       << "in Vaughan's model  "
       <<  Attributes::getReal(itsAttr[AttributesT::VKTHETA]) << endl;
    os << "* * Thermal velocity of Maxwellian distribution of secondaries in "
       << "Vaughan's model  "
       <<  Attributes::getReal(itsAttr[AttributesT::VVTHERMAL]) << endl;
    os << "* * VW denote the velocity scalar for Parallel plate benchmark  "
       <<  Attributes::getReal(itsAttr[AttributesT::VW]) << endl;
    os << "* * Material type number of the cavity surface for Furman-Pivi's model, "
       << "0 for copper, 1 for stainless steel  "
       <<  Attributes::getReal(itsAttr[AttributesT::SURFMATERIAL]) << endl;
}

void Distribution::PrintDistSurfAndCreate(Inform &os) const {

    os << "* Distribution type: SURFACERANDCREATE" << endl;
    os << "* " << endl;
    os << "* * Number of electrons initialized on the surface as primaries  "
       << Attributes::getReal(itsAttr[AttributesT::NPDARKCUR]) << endl;
    os << "* * Initialized electrons inward margin for surface emission  "
       << Attributes::getReal(itsAttr[AttributesT::INWARDMARGIN]) << endl;
    os << "* * E field threshold (MV), only in position r with E(r)>EINITHR that "
       << "particles will be initialized   "
       << Attributes::getReal(itsAttr[AttributesT::EINITHR]) << endl;
}

void Distribution::PrintEmissionModel(Inform &os) const {

    os << "* ------------- THERMAL EMITTANCE MODEL --------------------------------------------" << endl;

    switch (emissionModel_m) {

    case EmissionModelT::NONE:
        PrintEmissionModelNone(os);
        break;
    case EmissionModelT::ASTRA:
        PrintEmissionModelAstra(os);
        break;
    case EmissionModelT::NONEQUIL:
        PrintEmissionModelNonEquil(os);
        break;
    default:
        break;
    }

    os << "* ----------------------------------------------------------------------------------" << endl;

}

void Distribution::PrintEmissionModelAstra(Inform &os) const {
    os << "*  THERMAL EMITTANCE in ASTRA MODE" << endl;
    os << "*  Kinetic energy (thermal emittance) = "
       << std::abs(Attributes::getReal(itsAttr[AttributesT::EKIN]))
       << " [eV]  " << endl;
}

void Distribution::PrintEmissionModelNone(Inform &os) const {
    os << "*  THERMAL EMITTANCE in NONE MODE" << endl;
    os << "*  Kinetic energy added to longitudinal momentum = "
       << std::abs(Attributes::getReal(itsAttr[AttributesT::EKIN]))
       << " [eV]  " << endl;
}

void Distribution::PrintEmissionModelNonEquil(Inform &os) const {
    os << "*  THERMAL EMITTANCE in NONEQUIL MODE" << endl;
    os << "*  Cathode work function     = " << cathodeWorkFunc_m << " [eV]  " << endl;
    os << "*  Cathode Fermi energy      = " << cathodeFermiEnergy_m << " [eV]  " << endl;
    os << "*  Cathode temperature       = " << cathodeTemp_m << " [K]  " << endl;
    os << "*  Photocathode laser energy = " << laserEnergy_m << " [eV]  " << endl;
}

void Distribution::PrintEnergyBins(Inform &os) const {

    os << "* " << endl;
    for (int binIndex = numberOfEnergyBins_m - 1; binIndex >=0; binIndex--) {
        size_t sum = 0;
        for (int sampleIndex = 0; sampleIndex < numberOfSampleBins_m; sampleIndex++)
            sum += gsl_histogram_get(energyBinHist_m,
                                     binIndex * numberOfSampleBins_m + sampleIndex);

        os << "* Energy Bin # " << numberOfEnergyBins_m - binIndex
           << " contains " << sum << " particles" << endl;
    }
    os << "* " << endl;

}

bool Distribution::Rebin() {
    /*
     * Allow a rebin (numberOfEnergyBins_m = 0) if all particles are emitted or
     * injected.
     */
    if (!emitting_m) {
        numberOfEnergyBins_m = 0;
        return true;
    } else {
        return false;
    }
}

void Distribution::ReflectDistribution(size_t &numberOfParticles) {

    size_t currentNumPart = tOrZDist_m.size();
    for (size_t partIndex = 0; partIndex < currentNumPart; partIndex++) {
        xDist_m.push_back(-xDist_m.at(partIndex));
        pxDist_m.push_back(-pxDist_m.at(partIndex));
        yDist_m.push_back(-yDist_m.at(partIndex));
        pyDist_m.push_back(-pyDist_m.at(partIndex));
        tOrZDist_m.push_back(tOrZDist_m.at(partIndex));
        pzDist_m.push_back(pzDist_m.at(partIndex));
    }
    numberOfParticles = tOrZDist_m.size();
    reduce(numberOfParticles, numberOfParticles, OpAddAssign());
}

void Distribution::ScaleDistCoordinates() {

    for (size_t particleIndex = 0; particleIndex < tOrZDist_m.size(); particleIndex++) {
        xDist_m.at(particleIndex) *= Attributes::getReal(itsAttr[AttributesT::XMULT]);
        pxDist_m.at(particleIndex) *= Attributes::getReal(itsAttr[AttributesT::PXMULT]);
        yDist_m.at(particleIndex) *= Attributes::getReal(itsAttr[AttributesT::YMULT]);
        pyDist_m.at(particleIndex) *= Attributes::getReal(itsAttr[AttributesT::PYMULT]);

        if (emitting_m)
            tOrZDist_m.at(particleIndex) *= Attributes::getReal(itsAttr[AttributesT::TMULT]);
        else
            tOrZDist_m.at(particleIndex) *= Attributes::getReal(itsAttr[AttributesT::ZMULT]);

        pzDist_m.at(particleIndex) *= Attributes::getReal(itsAttr[AttributesT::PZMULT]);
    }
}

void Distribution::SetAttributes() {
    itsAttr[AttributesT::DISTRIBUTION]
        = Attributes::makeString("DISTRIBUTION","Distribution type: FROMFILE, "
                                 "GAUSS, "
                                 "BINOMIAL, "
                                 "FLATTOP, "
                                 "SURFACEEMISSION, "
                                 "SURFACERANDCREATE, "
                                 "GUNGAUSSFLATTOPTH, "
                                 "ASTRAFLATTOPTH, "
                                 "GAUSS, "
                                 "GAUSSMATCHED");

    itsAttr[AttributesT::LINE]
        = Attributes::makeString("LINE", "Beamline that contains a cyclotron or ring ", "");
    itsAttr[AttributesT::FMAPFN]
        = Attributes::makeString("FMAPFN", "File for reading fieldmap used to create matched distibution ", "");
    itsAttr[AttributesT::FMTYPE]
        = Attributes::makeString("FMTYPE", "File format for reading fieldmap used to create matched distribution ", "");
    itsAttr[AttributesT::EX]
        = Attributes::makeReal("EX", "Projected normalized emittance EX (m-rad), used to create matched distibution ", 1E-6);
    itsAttr[AttributesT::EY]
        = Attributes::makeReal("EY", "Projected normalized emittance EY (m-rad) used to create matched distibution ", 1E-6);
    itsAttr[AttributesT::ET]
        = Attributes::makeReal("ET", "Projected normalized emittance EY (m-rad) used to create matched distibution ", 1E-6);
    itsAttr[AttributesT::E2]
        = Attributes::makeReal("E2", "If E2<Eb, we compute the tunes from the beams energy Eb to E2 with dE=0.25 MeV ", 0.0);
    itsAttr[AttributesT::RESIDUUM]
        = Attributes::makeReal("RESIDUUM", "Residuum for the closed orbit finder and sigma matrix generator ", 1e-8);
    itsAttr[AttributesT::MAXSTEPSCO]
        = Attributes::makeReal("MAXSTEPSCO", "Maximum steps used to find closed orbit ", 100);
    itsAttr[AttributesT::MAXSTEPSSI]
        = Attributes::makeReal("MAXSTEPSSI", "Maximum steps used to find matched distribution ",2000);
    itsAttr[AttributesT::ORDERMAPS]
      = Attributes::makeReal("ORDERMAPS", "Order used in the field expansion ", 7);
    itsAttr[AttributesT::MAGSYM]
      = Attributes::makeReal("MAGSYM", "Number of sector magnets ", 0);
    itsAttr[AttributesT::RGUESS]
      = Attributes::makeReal("RGUESS", "Guess value of radius (m) for closed orbit finder ", -1);

    itsAttr[AttributesT::FNAME]
      = Attributes::makeString("FNAME", "File for reading in 6D particle "
			       "coordinates.", "");


    itsAttr[AttributesT::WRITETOFILE]
      = Attributes::makeBool("WRITETOFILE", "Write initial distribution to file.",
			     false);

    itsAttr[AttributesT::WEIGHT]
      = Attributes::makeReal("WEIGHT", "Weight of distribution when used in a "
			     "distribution list.", 1.0);

    itsAttr[AttributesT::INPUTMOUNITS]
        = Attributes::makeString("INPUTMOUNITS", "Tell OPAL what input units are for momentum."
                                 " Currently \"NONE\" or \"EV\".", "");

    // Attributes for beam emission.
    itsAttr[AttributesT::EMITTED]
        = Attributes::makeBool("EMITTED", "Emitted beam, from cathode, as opposed to "
                               "an injected beam.", false);
    itsAttr[AttributesT::EMISSIONSTEPS]
        = Attributes::makeReal("EMISSIONSTEPS", "Number of time steps to use during emission.",
                               1);
    itsAttr[AttributesT::EMISSIONMODEL]
        = Attributes::makeString("EMISSIONMODEL", "Model used to emit electrons from a "
                                 "photocathode.", "None");
    itsAttr[AttributesT::EKIN]
        = Attributes::makeReal("EKIN", "Kinetic energy used in ASTRA thermal emittance "
                               "model (eV). (Thermal energy added in with random "
                               "number generator.)", 1.0);
    itsAttr[AttributesT::ELASER]
        = Attributes::makeReal("ELASER", "Laser energy (eV) for photocathode "
                               "emission. (Default 255 nm light.)", 4.86);
    itsAttr[AttributesT::W]
        = Attributes::makeReal("W", "Workfunction of photocathode material (eV)."
                               " (Default atomically clean copper.)", 4.31);
    itsAttr[AttributesT::FE]
        = Attributes::makeReal("FE", "Fermi energy of photocathode material (eV)."
                               " (Default atomically clean copper.)", 7.0);
    itsAttr[AttributesT::CATHTEMP]
        = Attributes::makeReal("CATHTEMP", "Temperature of photocathode (K)."
                               " (Default 300 K.)", 300.0);
    itsAttr[AttributesT::NBIN]
        = Attributes::makeReal("NBIN", "Number of energy bins to use when doing a space "
                               "charge solve.", 0.0);

    // Parameters for shifting or scaling initial distribution.
    itsAttr[AttributesT::XMULT] = Attributes::makeReal("XMULT", "Multiplier for x dimension.", 1.0);
    itsAttr[AttributesT::YMULT] = Attributes::makeReal("YMULT", "Multiplier for y dimension.", 1.0);
    itsAttr[AttributesT::ZMULT] = Attributes::makeReal("TMULT", "Multiplier for z dimension.", 1.0);
    itsAttr[AttributesT::TMULT] = Attributes::makeReal("TMULT", "Multiplier for t dimension.", 1.0);

    itsAttr[AttributesT::PXMULT] = Attributes::makeReal("PXMULT", "Multiplier for px dimension.", 1.0);
    itsAttr[AttributesT::PYMULT] = Attributes::makeReal("PYMULT", "Multiplier for py dimension.", 1.0);
    itsAttr[AttributesT::PZMULT] = Attributes::makeReal("PZMULT", "Multiplier for pz dimension.", 1.0);

    itsAttr[AttributesT::OFFSETX]
        = Attributes::makeReal("OFFSETX", "Average x offset of distribution.", 0.0);
    itsAttr[AttributesT::OFFSETY]
        = Attributes::makeReal("OFFSETY", "Average y offset of distribution.", 0.0);
    itsAttr[AttributesT::OFFSETZ]
        = Attributes::makeReal("OFFSETZ", "Average z offset of distribution.", 0.0);
    itsAttr[AttributesT::OFFSETT]
        = Attributes::makeReal("OFFSETT", "Average t offset of distribution.", 0.0);

    itsAttr[AttributesT::OFFSETPX]
        = Attributes::makeReal("OFFSETPX", "Average px offset of distribution.", 0.0);
    itsAttr[AttributesT::OFFSETPY]
        = Attributes::makeReal("OFFSETPY", "Average py offset of distribution.", 0.0);
    itsAttr[AttributesT::OFFSETPZ]
        = Attributes::makeReal("OFFSETPZ", "Average pz offset of distribution.", 0.0);

    // Parameters for defining an initial distribution.
    itsAttr[AttributesT::SIGMAX] = Attributes::makeReal("SIGMAX", "SIGMAx (m)", 0.0);
    itsAttr[AttributesT::SIGMAY] = Attributes::makeReal("SIGMAY", "SIGMAy (m)", 0.0);
    itsAttr[AttributesT::SIGMAR] = Attributes::makeReal("SIGMAR", "SIGMAr (m)", 0.0);
    itsAttr[AttributesT::SIGMAZ] = Attributes::makeReal("SIGMAZ", "SIGMAz (m)", 0.0);
    itsAttr[AttributesT::SIGMAT] = Attributes::makeReal("SIGMAT", "SIGMAt (m)", 0.0);
    itsAttr[AttributesT::TPULSEFWHM]
        = Attributes::makeReal("TPULSEFWHM", "Pulse FWHM for emitted distribution.", 0.0);
    itsAttr[AttributesT::TRISE]
        = Attributes::makeReal("TRISE", "Rise time for emitted distribution.", 0.0);
    itsAttr[AttributesT::TFALL]
        = Attributes::makeReal("TFALL", "Fall time for emitted distribution.", 0.0);
    itsAttr[AttributesT::SIGMAPX] = Attributes::makeReal("SIGMAPX", "SIGMApx", 0.0);
    itsAttr[AttributesT::SIGMAPY] = Attributes::makeReal("SIGMAPY", "SIGMApy", 0.0);
    itsAttr[AttributesT::SIGMAPZ] = Attributes::makeReal("SIGMAPZ", "SIGMApz", 0.0);

    itsAttr[AttributesT::MX]
        = Attributes::makeReal("MX", "Defines the binomial distribution in x, "
                               "0.0 ... infinity.", 0.0);
    itsAttr[AttributesT::MY]
        = Attributes::makeReal("MY", "Defines the binomial distribution in y, "
                               "0.0 ... infinity.", 0.0);
    itsAttr[AttributesT::MZ]
        = Attributes::makeReal("MZ", "Defines the binomial distribution in z, "
                               "0.0 ... infinity.", 0.0);
    itsAttr[AttributesT::MT]
        = Attributes::makeReal("MT", "Defines the binomial distribution in t, "
                               "0.0 ... infinity.", 1.0);

    itsAttr[AttributesT::CUTOFFX] = Attributes::makeReal("CUTOFFX", "Distribution cutoff x "
                                                         "direction in units of sigma.", 3.0);
    itsAttr[AttributesT::CUTOFFY] = Attributes::makeReal("CUTOFFY", "Distribution cutoff r "
                                                         "direction in units of sigma.", 3.0);
    itsAttr[AttributesT::CUTOFFR] = Attributes::makeReal("CUTOFFR", "Distribution cutoff radial "
                                                         "direction in units of sigma.", 3.0);
    itsAttr[AttributesT::CUTOFFLONG]
        = Attributes::makeReal("CUTOFFLONG", "Distribution cutoff z or t direction in "
                               "units of sigma.", 3.0);
    itsAttr[AttributesT::CUTOFFPX] = Attributes::makeReal("CUTOFFPX", "Distribution cutoff px "
                                                          "dimension in units of sigma.", 3.0);
    itsAttr[AttributesT::CUTOFFPY] = Attributes::makeReal("CUTOFFPY", "Distribution cutoff py "
                                                          "dimension in units of sigma.", 3.0);
    itsAttr[AttributesT::CUTOFFPZ] = Attributes::makeReal("CUTOFFPZ", "Distribution cutoff pz "
                                                          "dimension in units of sigma.", 3.0);

    itsAttr[AttributesT::FTOSCAMPLITUDE]
        = Attributes::makeReal("FTOSCAMPLITUDE", "Amplitude of oscillations superimposed "
                               "on flat top portion of emitted GAUSS "
                               "distribtuion (in percent of flat top "
                               "amplitude)",0.0);
    itsAttr[AttributesT::FTOSCPERIODS]
        = Attributes::makeReal("FTOSCPERIODS", "Number of oscillations superimposed on "
                               "flat top portion of emitted GAUSS "
                               "distribution", 0.0);

    /*
     * TODO: Find out what these correlations really mean and write
     * good descriptions for them.
     */
    itsAttr[AttributesT::CORRX]
        = Attributes::makeReal("CORRX", "x/px correlation, (R12 in transport "
                               "notation).", 0.0);
    itsAttr[AttributesT::CORRY]
        = Attributes::makeReal("CORRY", "y/py correlation, (R34 in transport "
                               "notation).", 0.0);
    itsAttr[AttributesT::CORRZ]
        = Attributes::makeReal("CORRZ", "z/pz correlation, (R56 in transport "
                               "notation).", 0.0);
    itsAttr[AttributesT::CORRT]
        = Attributes::makeReal("CORRT", "t/pt correlation, (R56 in transport "
                               "notation).", 0.0);

    itsAttr[AttributesT::R51]
        = Attributes::makeReal("R51", "x/z correlation, (R51 in transport "
                               "notation).", 0.0);
    itsAttr[AttributesT::R52]
        = Attributes::makeReal("R52", "xp/z correlation, (R52 in transport "
                               "notation).", 0.0);

    itsAttr[AttributesT::R61]
        = Attributes::makeReal("R61", "x/pz correlation, (R61 in transport "
                               "notation).", 0.0);
    itsAttr[AttributesT::R62]
        = Attributes::makeReal("R62", "xp/pz correlation, (R62 in transport "
                               "notation).", 0.0);

    itsAttr[AttributesT::R]
        = Attributes::makeRealArray("R", "r correlation");

    // Parameters for using laser profile to generate a distribution.
    itsAttr[AttributesT::LASERPROFFN]
        = Attributes::makeString("LASERPROFFN", "File containing a measured laser image "
                                 "profile (x,y).", "");
    itsAttr[AttributesT::IMAGENAME]
        = Attributes::makeString("IMAGENAME", "Name of the laser image.", "");
    itsAttr[AttributesT::INTENSITYCUT]
        = Attributes::makeReal("INTENSITYCUT", "For background subtraction of laser "
                               "image profile, in % of max intensity.",
                               0.0);
    itsAttr[AttributesT::FLIPX]
        = Attributes::makeBool("FLIPX", "Flip laser profile horizontally", false);
    itsAttr[AttributesT::FLIPY]
        = Attributes::makeBool("FLIPY", "Flip laser profile vertically", false);
    itsAttr[AttributesT::ROTATE90]
        = Attributes::makeBool("ROTATE90", "Rotate laser profile 90 degrees counter clockwise", false);
    itsAttr[AttributesT::ROTATE180]
        = Attributes::makeBool("ROTATE180", "Rotate laser profile 180 degrees counter clockwise", false);
    itsAttr[AttributesT::ROTATE270]
        = Attributes::makeBool("ROTATE270", "Rotate laser profile 270 degrees counter clockwise", false);

    // Dark current and field emission model parameters.
    itsAttr[AttributesT::NPDARKCUR]
        = Attributes::makeReal("NPDARKCUR", "Number of dark current particles.", 1000.0);
    itsAttr[AttributesT::INWARDMARGIN]
        = Attributes::makeReal("INWARDMARGIN", "Inward margin of initialized dark "
                               "current particle positions.", 0.001);
    itsAttr[AttributesT::EINITHR]
        = Attributes::makeReal("EINITHR", "E field threshold (MV/m), only in position r "
                               "with E(r)>EINITHR that particles will be "
                               "initialized.", 0.0);
    itsAttr[AttributesT::FNA]
        = Attributes::makeReal("FNA", "Empirical constant A for Fowler-Nordheim "
                               "emission model.", 1.54e-6);
    itsAttr[AttributesT::FNB]
        = Attributes::makeReal("FNB", "Empirical constant B for Fowler-Nordheim "
                               "emission model.", 6.83e9);
    itsAttr[AttributesT::FNY]
        = Attributes::makeReal("FNY", "Constant for image charge effect parameter y(E) "
                               "in Fowler-Nordheim emission model.", 3.795e-5);
    itsAttr[AttributesT::FNVYZERO]
        = Attributes::makeReal("FNVYZERO", "Zero order constant for v(y) function in "
                               "Fowler-Nordheim emission model.", 0.9632);
    itsAttr[AttributesT::FNVYSECOND]
        = Attributes::makeReal("FNVYSECOND", "Second order constant for v(y) function "
                               "in Fowler-Nordheim emission model.", 1.065);
    itsAttr[AttributesT::FNPHIW]
        = Attributes::makeReal("FNPHIW", "Work function of gun surface material (eV).",
                               4.65);
    itsAttr[AttributesT::FNBETA]
        = Attributes::makeReal("FNBETA", "Field enhancement factor for Fowler-Nordheim "
                               "emission.", 50.0);
    itsAttr[AttributesT::FNFIELDTHR]
        = Attributes::makeReal("FNFIELDTHR", "Field threshold for Fowler-Nordheim "
                               "emission (MV/m).", 30.0);
    itsAttr[AttributesT::FNMAXEMI]
        = Attributes::makeReal("FNMAXEMI", "Maximum number of electrons emitted from a "
                               "single triangle for Fowler-Nordheim "
                               "emission.", 20.0);
    itsAttr[AttributesT::SECONDARYFLAG]
        = Attributes::makeReal("SECONDARYFLAG", "Select the secondary model type 0:no "
                               "secondary emission; 1:Furman-Pivi; "
                               "2 or larger: Vaughan's model.", 0.0);
    itsAttr[AttributesT::NEMISSIONMODE]
        = Attributes::makeBool("NEMISSIONMODE", "Secondary emission mode type true: "
                               "emit n true secondaries; false: emit "
                               "one particle with n times charge.", true);
    itsAttr[AttributesT::VSEYZERO]
        = Attributes::makeReal("VSEYZERO", "Sey_0 in Vaughan's model.", 0.5);
    itsAttr[AttributesT::VEZERO]
        = Attributes::makeReal("VEZERO", "Energy related to sey_0 in Vaughan's model "
                               "in eV.", 12.5);
    itsAttr[AttributesT::VSEYMAX]
        = Attributes::makeReal("VSEYMAX", "Sey max in Vaughan's model.", 2.22);
    itsAttr[AttributesT::VEMAX]
        = Attributes::makeReal("VEMAX", "Emax in Vaughan's model in eV.", 165.0);
    itsAttr[AttributesT::VKENERGY]
        = Attributes::makeReal("VKENERGY", "Fitting parameter denotes the roughness of "
                               "surface for impact energy in Vaughan's "
                               "model.", 1.0);
    itsAttr[AttributesT::VKTHETA]
        = Attributes::makeReal("VKTHETA", "Fitting parameter denotes the roughness of "
                               "surface for impact angle in Vaughan's "
                               "model.", 1.0);
    itsAttr[AttributesT::VVTHERMAL]
        = Attributes::makeReal("VVTHERMAL", "Thermal velocity of Maxwellian distribution "
                               "of secondaries in Vaughan's model.", 7.268929821 * 1e5);
    itsAttr[AttributesT::VW]
        = Attributes::makeReal("VW", "VW denote the velocity scalar for parallel plate "
                               "benchmark.", 1.0);
    itsAttr[AttributesT::SURFMATERIAL]
        = Attributes::makeReal("SURFMATERIAL", "Material type number of the cavity "
                               "surface for Furman-Pivi's model, 0 "
                               "for copper, 1 for stainless steel.", 0.0);

    /*
     *   Feature request Issue #14
     */

    itsAttr[AttributesT::ID1]
        = Attributes::makeRealArray("ID1", "User defined particle with ID=1: x,y,z,px,py,pz");
    itsAttr[AttributesT::ID2]
        = Attributes::makeRealArray("ID2", "User defined particle with ID=2: x,y,z,px,py,pz");

    /*
     * Legacy attributes (or ones that need to be implemented.)
     */

    // Parameters for emitting a distribution.
    itsAttr[ LegacyAttributesT::DEBIN]
        = Attributes::makeReal("DEBIN", "Energy band for a bin in keV that defines "
                               "when to combine bins. That is, when the energy "
                               "spread of all bins is below this number "
                               "combine bins into a single bin.", 1000000.0);
    itsAttr[ LegacyAttributesT::SBIN]
        = Attributes::makeReal("SBIN", "Number of sample bins to use per energy bin "
                               "when emitting a distribution.", 100.0);
    itsAttr[ LegacyAttributesT::TEMISSION]
        = Attributes::makeReal("TEMISSION", "Time it takes to emit distribution.", 0.0);

    itsAttr[ LegacyAttributesT::SIGLASER]
        = Attributes::makeReal("SIGLASER","Sigma of (uniform) laser spot size.", 0.0);
    itsAttr[ LegacyAttributesT::AG]
        = Attributes::makeReal("AG", "Acceleration gradient at photocathode surface "
                               "(MV/m).", 0.0);


    /*
     * Specific to type GAUSS and GUNGAUSSFLATTOPTH and ASTRAFLATTOPTH. These
     * last two distribution will eventually just become a subset of GAUSS.
     */
    itsAttr[ LegacyAttributesT::SIGMAPT] = Attributes::makeReal("SIGMAPT", "SIGMApt", 0.0);

    itsAttr[ LegacyAttributesT::TRANSVCUTOFF]
        = Attributes::makeReal("TRANSVCUTOFF", "Transverse cut-off in units of sigma.",
                               -3.0);

    itsAttr[ LegacyAttributesT::CUTOFF]
        = Attributes::makeReal("CUTOFF", "Longitudinal cutoff for Gaussian in units "
                               "of sigma.", 3.0);


    // Mixed use attributes (used by more than one distribution type).
    itsAttr[ LegacyAttributesT::Z]
        = Attributes::makeReal("Z", "Average longitudinal position of distribution "
                               "in z coordinates (m).", -99.0);
    itsAttr[ LegacyAttributesT::T]
        = Attributes::makeReal("T", "Average longitudinal position of distribution "
                               "in t coordinates (s).", 0.0);
    itsAttr[ LegacyAttributesT::PT]
        = Attributes::makeReal("PT", "Average longitudinal momentum of distribution.",
                               0.0);


    // Attributes that are not yet implemented.
    itsAttr[ LegacyAttributesT::ALPHAX]
        = Attributes::makeReal("ALPHAX", "Courant Snyder parameter.", 0.0);
    itsAttr[ LegacyAttributesT::ALPHAY]
        = Attributes::makeReal("ALPHAY", "Courant Snyder parameter.", 0.0);
    itsAttr[ LegacyAttributesT::BETAX]
        = Attributes::makeReal("BETAX", "Courant Snyder parameter.", 1.0);
    itsAttr[ LegacyAttributesT::BETAY]
        = Attributes::makeReal("BETAY", "Courant Snyder parameter.", 1.0);

    itsAttr[ LegacyAttributesT::DX]
        = Attributes::makeReal("DX", "Dispersion in x (R16 in Transport notation).", 0.0);
    itsAttr[ LegacyAttributesT::DDX]
        = Attributes::makeReal("DDX", "First derivative of DX.", 0.0);

    itsAttr[ LegacyAttributesT::DY]
        = Attributes::makeReal("DY", "Dispersion in y (R36 in Transport notation).", 0.0);
    itsAttr[ LegacyAttributesT::DDY]
        = Attributes::makeReal("DDY", "First derivative of DY.", 0.0);

    registerOwnership(AttributeHandler::STATEMENT);
}

void Distribution::SetFieldEmissionParameters() {

    darkCurrentParts_m = static_cast<size_t> (Attributes::getReal(itsAttr[AttributesT::NPDARKCUR]));
    darkInwardMargin_m = Attributes::getReal(itsAttr[AttributesT::INWARDMARGIN]);
    eInitThreshold_m = Attributes::getReal(itsAttr[AttributesT::EINITHR]);
    workFunction_m = Attributes::getReal(itsAttr[AttributesT::FNPHIW]);
    fieldEnhancement_m = Attributes::getReal(itsAttr[AttributesT::FNBETA]);
    maxFN_m = static_cast<size_t> (Attributes::getReal(itsAttr[AttributesT::FNMAXEMI]));
    fieldThrFN_m = Attributes::getReal(itsAttr[AttributesT::FNFIELDTHR]);
    paraFNA_m = Attributes::getReal(itsAttr[AttributesT::FNA]);
    paraFNB_m = Attributes::getReal(itsAttr[AttributesT::FNB]);
    paraFNY_m = Attributes::getReal(itsAttr[AttributesT::FNY]);
    paraFNVYZe_m = Attributes::getReal(itsAttr[AttributesT::FNVYZERO]);
    paraFNVYSe_m = Attributes::getReal(itsAttr[AttributesT::FNVYSECOND]);
    secondaryFlag_m = Attributes::getReal(itsAttr[AttributesT::SECONDARYFLAG]);
    ppVw_m = Attributes::getReal(itsAttr[AttributesT::VW]);
    vVThermal_m = Attributes::getReal(itsAttr[AttributesT::VVTHERMAL]);

}

void Distribution::SetDistToEmitted(bool emitted) {
    emitting_m = emitted;
}

void Distribution::SetDistType() {

    distT_m = Attributes::getString(itsAttr[AttributesT::DISTRIBUTION]);
    if (distT_m == "FROMFILE")
        distrTypeT_m = DistrTypeT::FROMFILE;
    else if(distT_m == "GAUSS")
        distrTypeT_m = DistrTypeT::GAUSS;
    else if(distT_m == "BINOMIAL")
        distrTypeT_m = DistrTypeT::BINOMIAL;
    else if(distT_m == "FLATTOP")
        distrTypeT_m = DistrTypeT::FLATTOP;
    else if(distT_m == "GUNGAUSSFLATTOPTH")
        distrTypeT_m = DistrTypeT::GUNGAUSSFLATTOPTH;
    else if(distT_m == "ASTRAFLATTOPTH")
        distrTypeT_m = DistrTypeT::ASTRAFLATTOPTH;
    else if(distT_m == "SURFACEEMISSION")
        distrTypeT_m = DistrTypeT::SURFACEEMISSION;
    else if(distT_m == "SURFACERANDCREATE")
        distrTypeT_m = DistrTypeT::SURFACERANDCREATE;
    else if(distT_m == "GAUSSMATCHED")
        distrTypeT_m = DistrTypeT::MATCHEDGAUSS;

}

void Distribution::SetEmissionTime(double &maxT, double &minT) {

    if (addedDistributions_m.size() == 0) {

        switch (distrTypeT_m) {

        case DistrTypeT::FLATTOP:
            tEmission_m = tPulseLengthFWHM_m + (cutoffR_m[2] - sqrt(2.0 * log(2.0)))
                * (sigmaTRise_m + sigmaTFall_m);
            break;

        case DistrTypeT::GAUSS:
            tEmission_m = tPulseLengthFWHM_m + (cutoffR_m[2] - sqrt(2.0 * log(2.0)))
                * (sigmaTRise_m + sigmaTFall_m);
            break;

        case DistrTypeT::ASTRAFLATTOPTH:
            /*
             * Don't do anything. Emission time is set during the distribution
             * creation. Only this distribution type does it this way. This is
             * a legacy feature.
             */
            break;

        case DistrTypeT::GUNGAUSSFLATTOPTH:
            tEmission_m = tPulseLengthFWHM_m + (cutoffR_m[2] - sqrt(2.0 * log(2.0)))
                * (sigmaTRise_m + sigmaTFall_m);
            break;

        default:
            /*
             * Increase maxT and decrease minT by 0.05% of total time
             * to ensure that no particles fall on the leading edge of
             * the first emission time step or the trailing edge of the last
             * emission time step.
             */
            double deltaT = maxT - minT;
            maxT += deltaT * 0.0005;
            minT -= deltaT * 0.0005;
            tEmission_m = (maxT - minT);
            break;
        }

    } else {
        /*
         * Increase maxT and decrease minT by 0.05% of total time
         * to ensure that no particles fall on the leading edge of
         * the first emission time step or the trailing edge of the last
         * emission time step.
         */
        double deltaT = maxT - minT;
        maxT += deltaT * 0.0005;
        minT -= deltaT * 0.0005;
        tEmission_m = (maxT - minT);
    }
    tBin_m = tEmission_m / numberOfEnergyBins_m;
}

void Distribution::SetDistParametersBinomial(double massIneV) {

    /*
     * Set Distribution parameters. Do all the necessary checks depending
     * on the input attributes.
     */
    correlationMatrix_m(1, 0) = Attributes::getReal(itsAttr[AttributesT::CORRX]);
    correlationMatrix_m(3, 2) = Attributes::getReal(itsAttr[AttributesT::CORRY]);
    correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[AttributesT::CORRT]);
    correlationMatrix_m(4, 0) = Attributes::getReal(itsAttr[AttributesT::R51]);
    correlationMatrix_m(4, 1) = Attributes::getReal(itsAttr[AttributesT::R52]);
    correlationMatrix_m(5, 0) = Attributes::getReal(itsAttr[AttributesT::R61]);
    correlationMatrix_m(5, 1) = Attributes::getReal(itsAttr[AttributesT::R62]);

    //CORRZ overrides CORRT. We initially use CORRT for legacy compatibility.
    if (Attributes::getReal(itsAttr[AttributesT::CORRZ]) != 0.0)
        correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[AttributesT::CORRZ]);


    sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAX])),
                        std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAY])),
                        std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAT])));

    // SIGMAZ overrides SIGMAT. We initially use SIGMAT for legacy compatibility.
    if (Attributes::getReal(itsAttr[AttributesT::SIGMAZ]) != 0.0)
        sigmaR_m[2] = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAZ]));

    sigmaP_m = Vector_t(std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPX])),
                        std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPY])),
                        std::abs(Attributes::getReal(itsAttr[LegacyAttributesT::SIGMAPT])));

    // SIGMAPZ overrides SIGMAPT. We initially use SIGMAPT for legacy compatibility.
    if (Attributes::getReal(itsAttr[AttributesT::SIGMAPZ]) != 0.0)
        sigmaP_m[2] = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPZ]));

    // Check what input units we are using for momentum.
    switch (inputMoUnits_m) {

    case InputMomentumUnitsT::EV:
        sigmaP_m[0] = ConverteVToBetaGamma(sigmaP_m[0], massIneV);
        sigmaP_m[1] = ConverteVToBetaGamma(sigmaP_m[1], massIneV);
        sigmaP_m[2] = ConverteVToBetaGamma(sigmaP_m[2], massIneV);
        break;

    default:
        break;
    }

    mBinomial_m = Vector_t(std::abs(Attributes::getReal(itsAttr[AttributesT::MX])),
                           std::abs(Attributes::getReal(itsAttr[AttributesT::MY])),
                           std::abs(Attributes::getReal(itsAttr[AttributesT::MT])));

    cutoffR_m = Vector_t(Attributes::getReal(itsAttr[AttributesT::CUTOFFX]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFY]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFLONG]));

    cutoffP_m = Vector_t(Attributes::getReal(itsAttr[AttributesT::CUTOFFPX]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFPY]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFPZ]));

    if (mBinomial_m[2] == 0.0
        || Attributes::getReal(itsAttr[AttributesT::MZ]) != 0.0)
        mBinomial_m[2] = std::abs(Attributes::getReal(itsAttr[AttributesT::MZ]));

    if (emitting_m) {
        sigmaR_m[2] = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAT]));
        mBinomial_m[2] = std::abs(Attributes::getReal(itsAttr[AttributesT::MT]));
        correlationMatrix_m(5, 4) = std::abs(Attributes::getReal(itsAttr[AttributesT::CORRT]));
    }
}

void Distribution::SetDistParametersFlattop(double massIneV) {

    /*
     * Set distribution parameters. Do all the necessary checks depending
     * on the input attributes.
     */
    sigmaP_m = Vector_t(std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPX])),
                        std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPY])),
                        std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPZ])));

    // Check what input units we are using for momentum.
    switch (inputMoUnits_m) {

    case InputMomentumUnitsT::EV:
        sigmaP_m[0] = ConverteVToBetaGamma(sigmaP_m[0], massIneV);
        sigmaP_m[1] = ConverteVToBetaGamma(sigmaP_m[1], massIneV);
        sigmaP_m[2] = ConverteVToBetaGamma(sigmaP_m[2], massIneV);
        break;

    default:
        break;
    }

    cutoffR_m = Vector_t(Attributes::getReal(itsAttr[AttributesT::CUTOFFX]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFY]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFLONG]));

    correlationMatrix_m(1, 0) = Attributes::getReal(itsAttr[AttributesT::CORRX]);
    correlationMatrix_m(3, 2) = Attributes::getReal(itsAttr[AttributesT::CORRY]);
    correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[AttributesT::CORRT]);

    // CORRZ overrides CORRT.
    if (Attributes::getReal(itsAttr[AttributesT::CORRZ]) != 0.0)
        correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[AttributesT::CORRZ]);


    if (emitting_m) {
        INFOMSG("emitting"<<endl);
        sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAX])),
                            std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAY])),
                            0.0);

        sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAT]));
        sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAT]));

        tPulseLengthFWHM_m = std::abs(Attributes::getReal(itsAttr[AttributesT::TPULSEFWHM]));

        /*
         * If TRISE and TFALL are defined > 0.0 then these attributes
         * override SIGMAT.
         */
        if (std::abs(Attributes::getReal(itsAttr[AttributesT::TRISE])) > 0.0
            || std::abs(Attributes::getReal(itsAttr[AttributesT::TFALL])) > 0.0) {

            double timeRatio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));

            if (std::abs(Attributes::getReal(itsAttr[AttributesT::TRISE])) > 0.0)
                sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[AttributesT::TRISE]))
                    / timeRatio;

            if (std::abs(Attributes::getReal(itsAttr[AttributesT::TFALL])) > 0.0)
                sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[AttributesT::TFALL]))
                    / timeRatio;

        }

        // For an emitted beam, the longitudinal cutoff >= 0.
        cutoffR_m[2] = std::abs(cutoffR_m[2]);

    } else {

        sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAX])),
                            std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAY])),
                            std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAZ])));

    }

    if (std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAR])) > 0.0) {
        sigmaR_m[0] = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAR]));
        sigmaR_m[1] = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAR]));
    }

    // Set laser profile/
    laserProfileFileName_m = Attributes::getString(itsAttr[AttributesT::LASERPROFFN]);
    if (!(laserProfileFileName_m == std::string(""))) {
        laserImageName_m = Attributes::getString(itsAttr[AttributesT::IMAGENAME]);
        laserIntensityCut_m = std::abs(Attributes::getReal(itsAttr[AttributesT::INTENSITYCUT]));
        short flags = 0;
        if (Attributes::getBool(itsAttr[AttributesT::FLIPX])) flags |= LaserProfile::FLIPX;
        if (Attributes::getBool(itsAttr[AttributesT::FLIPY])) flags |= LaserProfile::FLIPY;
        if (Attributes::getBool(itsAttr[AttributesT::ROTATE90])) flags |= LaserProfile::ROTATE90;
        if (Attributes::getBool(itsAttr[AttributesT::ROTATE180])) flags |= LaserProfile::ROTATE180;
        if (Attributes::getBool(itsAttr[AttributesT::ROTATE270])) flags |= LaserProfile::ROTATE270;

        laserProfile_m = new LaserProfile(laserProfileFileName_m,
                                          laserImageName_m,
                                          laserIntensityCut_m,
                                          flags);
    }

    // Legacy for ASTRAFLATTOPTH.
    if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH)
        tRise_m = std::abs(Attributes::getReal(itsAttr[AttributesT::TRISE]));

}

void Distribution::SetDistParametersGauss(double massIneV) {

    /*
     * Set distribution parameters. Do all the necessary checks depending
     * on the input attributes.
     * In case of DistrTypeT::MATCHEDGAUSS we only need to set the cutoff parameters
     */


    cutoffP_m = Vector_t(Attributes::getReal(itsAttr[AttributesT::CUTOFFPX]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFPY]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFPZ]));


    cutoffR_m = Vector_t(Attributes::getReal(itsAttr[AttributesT::CUTOFFX]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFY]),
                         Attributes::getReal(itsAttr[AttributesT::CUTOFFLONG]));

    if  (distrTypeT_m != DistrTypeT::MATCHEDGAUSS) {
        sigmaP_m = Vector_t(std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPX])),
                            std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPY])),
                            std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPZ])));

        // SIGMAPZ overrides SIGMAPT. We initially use SIGMAPT for legacy compatibility.
        if (Attributes::getReal(itsAttr[AttributesT::SIGMAPZ]) != 0.0)
            sigmaP_m[2] = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAPZ]));

        // Check what input units we are using for momentum.
        switch (inputMoUnits_m) {

        case InputMomentumUnitsT::EV:
            sigmaP_m[0] = ConverteVToBetaGamma(sigmaP_m[0], massIneV);
            sigmaP_m[1] = ConverteVToBetaGamma(sigmaP_m[1], massIneV);
            sigmaP_m[2] = ConverteVToBetaGamma(sigmaP_m[2], massIneV);
            break;

        default:
            break;
        }

        std::vector<double> cr = Attributes::getRealArray(itsAttr[AttributesT::R]);

        if(cr.size()>0) {
            if(cr.size() == 15) {
                *gmsg << "* Use r to specify correlations" << endl;
                unsigned int k = 0;
                for (unsigned int i = 0; i < 5; ++ i) {
                    for (unsigned int j = i + 1; j < 6; ++ j, ++ k) {
                        correlationMatrix_m(j, i) = cr.at(k);
                    }
                }

            }
            else {
                throw OpalException("Distribution::SetDistParametersGauss",
                                    "Inconsitent set of correlations specified, check manual");
            }
        }
        else {
            correlationMatrix_m(1, 0) = Attributes::getReal(itsAttr[AttributesT::CORRX]);
            correlationMatrix_m(3, 2) = Attributes::getReal(itsAttr[AttributesT::CORRY]);
            correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[AttributesT::CORRT]);
            correlationMatrix_m(4, 0) = Attributes::getReal(itsAttr[AttributesT::R51]);
            correlationMatrix_m(4, 1) = Attributes::getReal(itsAttr[AttributesT::R52]);
            correlationMatrix_m(5, 0) = Attributes::getReal(itsAttr[AttributesT::R61]);
            correlationMatrix_m(5, 1) = Attributes::getReal(itsAttr[AttributesT::R62]);

            // CORRZ overrides CORRT.
            if (Attributes::getReal(itsAttr[AttributesT::CORRZ]) != 0.0)
                correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[AttributesT::CORRZ]);
        }
    }

    if (emitting_m) {

        sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAX])),
                            std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAY])),
                            0.0);

        sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAT]));
        sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAT]));

        tPulseLengthFWHM_m = std::abs(Attributes::getReal(itsAttr[AttributesT::TPULSEFWHM]));

        /*
         * If TRISE and TFALL are defined then these attributes
         * override SIGMAT.
         */
        if (std::abs(Attributes::getReal(itsAttr[AttributesT::TRISE])) > 0.0
            || std::abs(Attributes::getReal(itsAttr[AttributesT::TFALL])) > 0.0) {

            double timeRatio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));

            if (std::abs(Attributes::getReal(itsAttr[AttributesT::TRISE])) > 0.0)
                sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[AttributesT::TRISE]))
                    / timeRatio;

            if (std::abs(Attributes::getReal(itsAttr[AttributesT::TFALL])) > 0.0)
                sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[AttributesT::TFALL]))
                    / timeRatio;

        }

        // For and emitted beam, the longitudinal cutoff >= 0.
        cutoffR_m[2] = std::abs(cutoffR_m[2]);

    } else {
        if  (distrTypeT_m != DistrTypeT::MATCHEDGAUSS) {
            sigmaR_m = Vector_t(std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAX])),
                                std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAY])),
                                std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAT])));

            // SIGMAZ overrides SIGMAT. We initially use SIGMAT for legacy compatibility.
            if (Attributes::getReal(itsAttr[AttributesT::SIGMAZ]) != 0.0)
                sigmaR_m[2] = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAZ]));

        }
    }

    if (std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAR])) > 0.0) {
        sigmaR_m[0] = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAR]));
        sigmaR_m[1] = std::abs(Attributes::getReal(itsAttr[AttributesT::SIGMAR]));
        cutoffR_m[0] = Attributes::getReal(itsAttr[AttributesT::CUTOFFR]);
        cutoffR_m[1] = Attributes::getReal(itsAttr[AttributesT::CUTOFFR]);
    }
}

void Distribution::SetupEmissionModel(PartBunch &beam) {

    std::string model = Attributes::getString(itsAttr[AttributesT::EMISSIONMODEL]);
    if (model == "ASTRA")
        emissionModel_m = EmissionModelT::ASTRA;
    else if (model == "NONEQUIL")
        emissionModel_m = EmissionModelT::NONEQUIL;
    else
        emissionModel_m = EmissionModelT::NONE;

    /*
     * The ASTRAFLATTOPTH  of GUNGAUSSFLATTOPTH distributions always uses the
     * ASTRA emission model.
     */
    if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH
        || distrTypeT_m == DistrTypeT::GUNGAUSSFLATTOPTH)
        emissionModel_m = EmissionModelT::ASTRA;

    switch (emissionModel_m) {

    case EmissionModelT::NONE:
        SetupEmissionModelNone(beam);
        break;
    case EmissionModelT::ASTRA:
        SetupEmissionModelAstra(beam);
        break;
    case EmissionModelT::NONEQUIL:
        SetupEmissionModelNonEquil();
        break;
    default:
        break;
    }

}

void Distribution::SetupEmissionModelAstra(PartBunch &beam) {

    double wThermal = std::abs(Attributes::getReal(itsAttr[AttributesT::EKIN]));
    pTotThermal_m = ConverteVToBetaGamma(wThermal, beam.getM());
    pmean_m = Vector_t(0.0, 0.0, 0.5 * pTotThermal_m);

    gsl_rng_env_setup();
    randGenEmit_m = gsl_rng_alloc(gsl_rng_default);
}

void Distribution::SetupEmissionModelNone(PartBunch &beam) {

    double wThermal = std::abs(Attributes::getReal(itsAttr[AttributesT::EKIN]));
    pTotThermal_m = ConverteVToBetaGamma(wThermal, beam.getM());
    double avgPz = std::accumulate(pzDist_m.begin(), pzDist_m.end(), 0.0);
    size_t numParticles = pzDist_m.size();
    reduce(avgPz, avgPz, OpAddAssign());
    reduce(numParticles, numParticles, OpAddAssign());

    pmean_m = Vector_t(0.0, 0.0, pTotThermal_m + avgPz / numParticles);
}

void Distribution::SetupEmissionModelNonEquil() {

    cathodeWorkFunc_m = std::abs(Attributes::getReal(itsAttr[AttributesT::W]));
    laserEnergy_m = std::abs(Attributes::getReal(itsAttr[AttributesT::ELASER]));
    cathodeFermiEnergy_m = std::abs(Attributes::getReal(itsAttr[AttributesT::FE]));
    cathodeTemp_m = std::abs(Attributes::getReal(itsAttr[AttributesT::CATHTEMP]));

    /*
     * Upper limit on energy distribution theoretically goes to infinity.
     * Practically we limit to a probability of 1 part in a billion.
     */
    emitEnergyUpperLimit_m = cathodeFermiEnergy_m
        + Physics::kB * cathodeTemp_m * log(1.0e9 - 1.0);

    gsl_rng_env_setup();
    randGenEmit_m = gsl_rng_alloc(gsl_rng_default);

    // TODO: get better estimate of pmean
    pmean_m = 0.5 * (cathodeWorkFunc_m + emitEnergyUpperLimit_m);
}

void Distribution::SetupEnergyBins(double maxTOrZ, double minTOrZ) {

    energyBinHist_m = gsl_histogram_alloc(numberOfSampleBins_m * numberOfEnergyBins_m);

    if (emitting_m)
        gsl_histogram_set_ranges_uniform(energyBinHist_m, -tEmission_m, 0.0);
    else
        gsl_histogram_set_ranges_uniform(energyBinHist_m, minTOrZ, maxTOrZ);

}

void Distribution::SetupParticleBins(double massIneV, PartBunch &beam) {

    numberOfEnergyBins_m
        = static_cast<int>(std::abs(Attributes::getReal(itsAttr[AttributesT::NBIN])));

    if (numberOfEnergyBins_m > 0) {
        if (energyBins_m)
            delete energyBins_m;

        int sampleBins = static_cast<int>(std::abs(Attributes::getReal(itsAttr[LegacyAttributesT::SBIN])));
        energyBins_m = new PartBins(numberOfEnergyBins_m, sampleBins);

        double dEBins = Attributes::getReal(itsAttr[ LegacyAttributesT::DEBIN]);
        energyBins_m->setRebinEnergy(dEBins);

        if (Attributes::getReal(itsAttr[ LegacyAttributesT::PT])!=0.0)
            WARNMSG("PT & PZ are obsolet and will be ignored. The moments of the beam is defined with PC or use OFFSETPZ" << endl);

        // we get gamma from PC of the beam
        double gamma = 0.0;
        const double pz    = beam.getP()/beam.getM();
        gamma = sqrt(pow(pz, 2.0) + 1.0);
        energyBins_m->setGamma(gamma);

    } else {
        energyBins_m = NULL;
    }
}

void Distribution::ShiftBeam(double &maxTOrZ, double &minTOrZ) {

    std::vector<double>::iterator tOrZIt;
    if (emitting_m) {

        if (addedDistributions_m.size() == 0) {

            if (distrTypeT_m == DistrTypeT::ASTRAFLATTOPTH) {
                for (tOrZIt = tOrZDist_m.begin(); tOrZIt != tOrZDist_m.end(); tOrZIt++)
                    *tOrZIt -= tEmission_m / 2.0;

                minTOrZ -= tEmission_m / 2.0;
                maxTOrZ -= tEmission_m / 2.0;
            } else if (distrTypeT_m == DistrTypeT::GAUSS
                       || distrTypeT_m == DistrTypeT::FLATTOP
                       || distrTypeT_m == DistrTypeT::GUNGAUSSFLATTOPTH) {
                for (tOrZIt = tOrZDist_m.begin(); tOrZIt != tOrZDist_m.end(); tOrZIt++)
                    *tOrZIt -= tEmission_m;

                minTOrZ -= tEmission_m;
                maxTOrZ -= tEmission_m;
            } else {
                for (tOrZIt = tOrZDist_m.begin(); tOrZIt != tOrZDist_m.end(); tOrZIt++)
                    *tOrZIt -= maxTOrZ;

                minTOrZ -= maxTOrZ;
                maxTOrZ -= maxTOrZ;
            }

        } else {
            for (tOrZIt = tOrZDist_m.begin(); tOrZIt != tOrZDist_m.end(); tOrZIt++)
                *tOrZIt -= maxTOrZ;

            minTOrZ -= maxTOrZ;
            maxTOrZ -= maxTOrZ;
        }

    } else {
        if (minTOrZ < 0.0) {
            for (tOrZIt = tOrZDist_m.begin(); tOrZIt != tOrZDist_m.end(); tOrZIt++)
                *tOrZIt -= minTOrZ;
            maxTOrZ -= minTOrZ;
            minTOrZ -= minTOrZ;
        }
    }

}

void Distribution::ShiftDistCoordinates(double massIneV) {

    double deltaX = Attributes::getReal(itsAttr[AttributesT::OFFSETX]);
    double deltaY = Attributes::getReal(itsAttr[AttributesT::OFFSETY]);

    /*
     * OFFSETZ overrides T if it is nonzero. We initially use T
     * for legacy compatiblity. OFFSETT always overrides T, even
     * when zero, for an emitted beam.
     */
    double deltaTOrZ = Attributes::getReal(itsAttr[ LegacyAttributesT::T]);
    if (emitting_m)
        deltaTOrZ = Attributes::getReal(itsAttr[AttributesT::OFFSETT]);
    else {
        if (Attributes::getReal(itsAttr[AttributesT::OFFSETZ]) != 0.0)
            deltaTOrZ = Attributes::getReal(itsAttr[AttributesT::OFFSETZ]);
    }

    double deltaPx = Attributes::getReal(itsAttr[AttributesT::OFFSETPX]);
    double deltaPy = Attributes::getReal(itsAttr[AttributesT::OFFSETPY]);
    double deltaPz = Attributes::getReal(itsAttr[AttributesT::OFFSETPZ]);

    if (Attributes::getReal(itsAttr[LegacyAttributesT::PT])!=0.0)
        WARNMSG("PT & PZ are obsolet and will be ignored. The moments of the beam is defined with PC" << endl);

    // Check input momentum units.
    switch (inputMoUnits_m) {
    case InputMomentumUnitsT::EV:
        deltaPx = ConverteVToBetaGamma(deltaPx, massIneV);
        deltaPy = ConverteVToBetaGamma(deltaPy, massIneV);
        deltaPz = ConverteVToBetaGamma(deltaPz, massIneV);
        break;
    default:
        break;
    }

    for (size_t particleIndex = 0; particleIndex < tOrZDist_m.size(); particleIndex++) {
        xDist_m.at(particleIndex) += deltaX;
        pxDist_m.at(particleIndex) += deltaPx;
        yDist_m.at(particleIndex) += deltaY;
        pyDist_m.at(particleIndex) += deltaPy;
        tOrZDist_m.at(particleIndex) += deltaTOrZ;
        pzDist_m.at(particleIndex) += deltaPz;
    }
}

void Distribution::WriteOutFileHeader() {

    if (Attributes::getBool(itsAttr[AttributesT::WRITETOFILE])) {

        if (Ippl::myNode() == 0) {
            *gmsg << "* **********************************************************************************" << endl;
            *gmsg << "* Write initial distribution to file \"data/distribution.data\"" << endl;
            *gmsg << "* **********************************************************************************" << endl;
            std::ofstream outputFile("data/distribution.data");
            if (outputFile.bad()) {
                *gmsg << "Unable to open output file \"data/distribution.data\"" << endl;
            } else {
                outputFile.setf(std::ios::left);
                outputFile << "# ";
                if (emitting_m) {
                    outputFile.width(17);
                    outputFile << "x [m]";
                    outputFile.width(17);
                    outputFile << "px [betax gamma]";
                    outputFile.width(17);
                    outputFile << "y [m]";
                    outputFile.width(17);
                    outputFile << "py [betay gamma]";
                    outputFile.width(17);
                    outputFile << "t [s]";
                    outputFile.width(17);
                    outputFile << "pz [betaz gamma]" ;
                    outputFile.width(17);
                    outputFile << "Bin Number" << std::endl;
                } else {
                    outputFile.width(17);
                    outputFile << "x [m]";
                    outputFile.width(17);
                    outputFile << "px [betax gamma]";
                    outputFile.width(17);
                    outputFile << "y [m]";
                    outputFile.width(17);
                    outputFile << "py [betay gamma]";
                    outputFile.width(17);
                    outputFile << "t [s]";
                    outputFile.width(17);
                    outputFile << "pz [betaz gamma]";
                    if (numberOfEnergyBins_m > 0) {
                        outputFile.width(17);
                        outputFile << "Bin Number";
                    }
                    outputFile << std::endl;
                }
            }
            outputFile.close();
        }
    }
}

void Distribution::WriteOutFileEmission() {

    if (Attributes::getBool(itsAttr[AttributesT::WRITETOFILE])) {

        // Gather particles to be written onto node 0.
        int currentNode = 1;
        for (int nodeIndex = 1; nodeIndex < Ippl::getNodes(); nodeIndex++) {

            size_t numberOfParticles = 0;
            if (Ippl::myNode() == currentNode) {
                if (!xWrite_m.empty())
                    numberOfParticles = xWrite_m.size();
            }
            reduce(numberOfParticles, numberOfParticles, OpAddAssign());

            for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

                double x = 0.0;
                double px = 0.0;
                double y = 0.0;
                double py = 0.0;
                double tOrZ = 0.0;
                double pz = 0.0;
                size_t bin = 0;
                if (Ippl::myNode() == currentNode) {
                    x = xWrite_m.at(partIndex);
                    px = pxWrite_m.at(partIndex);
                    y = yWrite_m.at(partIndex);
                    py = pyWrite_m.at(partIndex);
                    tOrZ = tOrZWrite_m.at(partIndex);
                    pz = pzWrite_m.at(partIndex);
                    bin = binWrite_m.at(partIndex);
                }
                reduce(x, x, OpAddAssign());
                reduce(px, px, OpAddAssign());
                reduce(y, y, OpAddAssign());
                reduce(py, py, OpAddAssign());
                reduce(tOrZ, tOrZ, OpAddAssign());
                reduce(pz, pz, OpAddAssign());
                reduce(bin, bin, OpAddAssign());

                if (Ippl::myNode() == 0) {
                    xWrite_m.push_back(x);
                    pxWrite_m.push_back(px);
                    yWrite_m.push_back(y);
                    pyWrite_m.push_back(py);
                    tOrZWrite_m.push_back(tOrZ);
                    pzWrite_m.push_back(pz);
                    binWrite_m.push_back(bin);
                }

            }
            currentNode++;
        }

        // Write to file only on node 0.
        if (Ippl::myNode() == 0 && !xWrite_m.empty() > 0) {
            std::ofstream outputFile("data/distribution.data", std::fstream::app);
            if (outputFile.bad()) {
                *gmsg << "Unable to write to file \"data/distribution.data\"" << endl;
            } else {

                outputFile.precision(9);
                outputFile.setf(std::ios::scientific);
                outputFile.setf(std::ios::right);
                for (size_t partIndex = 0; partIndex < xWrite_m.size(); partIndex++) {

                    outputFile.width(17);
                    outputFile << xWrite_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << pxWrite_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << yWrite_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << pyWrite_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << tOrZWrite_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << pzWrite_m.at(partIndex);
                    outputFile.width(17);
                    outputFile << binWrite_m.at(partIndex) << std::endl;

                }

            }
            outputFile.close();
        }
    }

    // Clear write vectors.
    xWrite_m.clear();
    pxWrite_m.clear();
    yWrite_m.clear();
    pyWrite_m.clear();
    tOrZWrite_m.clear();
    pzWrite_m.clear();
    binWrite_m.clear();

}

void Distribution::WriteOutFileInjection() {

    if (Attributes::getBool(itsAttr[AttributesT::WRITETOFILE])) {

        // Nodes take turn writing particles to file.
        for (int nodeIndex = 0; nodeIndex < Ippl::getNodes(); nodeIndex++) {

            // Write to file if its our turn.
            size_t numberOfParticles = 0;
            if (Ippl::myNode() == nodeIndex) {

                std::ofstream outputFile("data/distribution.data", std::fstream::app);
                if (outputFile.bad()) {
		     *gmsg << "Node " << Ippl::myNode()
			   << " unable to write to file \"data/distribution.data\""
			   << endl;
                } else {

                    outputFile.precision(9);
                    outputFile.setf(std::ios::scientific);
                    outputFile.setf(std::ios::right);

                    numberOfParticles = tOrZDist_m.size();
                    for (size_t partIndex = 0; partIndex < numberOfParticles; partIndex++) {

                        outputFile.width(17);
                        outputFile << xDist_m.at(partIndex);
                        outputFile.width(17);
                        outputFile << pxDist_m.at(partIndex);
                        outputFile.width(17);
                        outputFile << yDist_m.at(partIndex);
                        outputFile.width(17);
                        outputFile << pyDist_m.at(partIndex);
                        outputFile.width(17);
                        outputFile << tOrZDist_m.at(partIndex);
                        outputFile.width(17);
                        outputFile << pzDist_m.at(partIndex);
                        if (numberOfEnergyBins_m > 0) {
                            size_t binNumber = FindEBin(tOrZDist_m.at(partIndex));
                            outputFile.width(17);
                            outputFile << binNumber;
                        }
                        outputFile << std::endl;

                    }
                }
                outputFile.close();
            }

            // Wait for writing node before moving on.
            reduce(numberOfParticles, numberOfParticles, OpAddAssign());
        }
    }
}
// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
