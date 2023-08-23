//
// Class Distribution
//   This class defines the initial beam that is injected or emitted into the simulation.
//
// Copyright (c) 2008 - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Distribution/Distribution.h"

#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/PartBunch.h"
#include "BasicActions/Option.h"
#include "DataSource/DataConnect.h"
#include "Distribution/LaserProfile.h"
#include "Elements/OpalBeamline.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/H5PartWrapper.h"
#include "Utilities/EarlyLeaveException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Utility/IpplTimings.h"

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_erf.h>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include <sys/time.h>

#include <cmath>
#include <cfloat>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>

extern Inform *gmsg;

constexpr double SMALLESTCUTOFF = 1e-12;

namespace {
    matrix_t getUnit6x6() {
        matrix_t unit6x6(6, 6, 0.0); // Initialize a 6x6 matrix with all elements as 0.0
        for (unsigned int i = 0; i < 6u; ++i) {
            unit6x6(i, i) = 1.0; // Set diagonal elements to 1.0
        }
        return unit6x6;
    }
}

Distribution::Distribution():
    Definition(DISTRIBUTION::SIZE, "DISTRIBUTION",
               "The DISTRIBUTION statement defines data for the 6D particle distribution."),
    distrTypeT_m(DistributionType::NODIST)
{
    itsAttr[DISTRIBUTION::TYPE] = Attributes::makePredefinedString("TYPE","Distribution type.", {"GAUSS","FROMFILE"});
   
    itsAttr[DISTRIBUTION::FNAME] = Attributes::makeString("FNAME", "File for reading in 6D particle coordinates.", "");

    // Parameters for defining an initial distribution.
    itsAttr[DISTRIBUTION::SIGMAX] = Attributes::makeReal("SIGMAX", "SIGMAx (m)", 0.0);
    itsAttr[DISTRIBUTION::SIGMAY] = Attributes::makeReal("SIGMAY", "SIGMAy (m)", 0.0);
    itsAttr[DISTRIBUTION::SIGMAZ] = Attributes::makeReal("SIGMAZ", "SIGMAz (m)", 0.0);
    itsAttr[DISTRIBUTION::SIGMAPX] = Attributes::makeReal("SIGMAPX", "SIGMApx", 0.0);
    itsAttr[DISTRIBUTION::SIGMAPY] = Attributes::makeReal("SIGMAPY", "SIGMApy", 0.0);
    itsAttr[DISTRIBUTION::SIGMAPZ] = Attributes::makeReal("SIGMAPZ", "SIGMApz", 0.0);

    registerOwnership(AttributeHandler::STATEMENT);
}


Distribution::Distribution(const std::string &name, Distribution *parent):
    Definition(name, parent)
{
    gsl_rng_env_setup();
    randGen_m = gsl_rng_alloc(gsl_rng_default);
}

Distribution::~Distribution() {

    gsl_rng_free(randGen_m);
}


/**
 * Calculate the local number of particles evenly and adjust node 0
 * such that n is matched exactly.
 * @param n total number of particles
 * @return n / #cores
 * @param
 */

size_t Distribution::getNumOfLocalParticlesToCreate(size_t n) {

    size_t locNumber = n / Ippl::getNodes();

    // make sure the total number is exact
    size_t remainder  = n % Ippl::getNodes();
    size_t myNode = Ippl::myNode();
    if (myNode < remainder)
        ++ locNumber;

    return locNumber;
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

void Distribution::createOpalT(PartBunch<double, 3> *beam,
                               std::vector<Distribution *> addedDistributions,
                               size_t &numberOfParticles) {


}

void Distribution::create(size_t &numberOfParticles, double massIneV, double charge) {

    /*
     * If Options::cZero is true than we reflect generated distribution
     * to ensure that the transverse averages are 0.0.
     *
     * For now we just cut the number of generated particles in half.
     */
    size_t numberOfLocalParticles = numberOfParticles;
    numberOfLocalParticles = (numberOfParticles + 1) / 2;
    
    size_t mySeed = Options::seed;

    if (Options::seed == -1) {
        struct timeval tv;
        gettimeofday(&tv,0);
        mySeed = tv.tv_sec + tv.tv_usec;
    }

    *gmsg << level2 << "* Generation of distribution with seed = " << mySeed << "\n"
          << "* isn't scalable with number of particles and cores." << endl;

    gsl_rng_set(randGen_m, mySeed);

    switch (distrTypeT_m) {

    case DistributionType::GAUSS:
        createDistributionGauss(numberOfLocalParticles, massIneV);
        break;
    default:
        throw OpalException("Distribution::create",
                            "Unknown \"TYPE\" of \"DISTRIBUTION\"");
    }

    if (Options::seed != -1)
        Options::seed = gsl_rng_uniform_int(randGen_m, gsl_rng_max(randGen_m));

}

Distribution *Distribution::find(const std::string &name) {
    Distribution *dist = dynamic_cast<Distribution *>(OpalData::getInstance()->find(name));

    if (dist == 0) {
        throw OpalException("Distribution::find()", "Distribution \"" + name + "\" not found.");
    }

    return dist;
}

Inform &Distribution::printInfo(Inform &os) const {

    os << "\n"
       << "* ************* D I S T R I B U T I O N ********************************************" << endl;
    os << "* " << endl;
    if (OpalData::getInstance()->inRestartRun()) {
        os << "* In restart. Distribution read in from .h5 file." << endl;
    } else {
        printDist(os, 0);

        os << "* " << endl;
        os << "* Distribution is injected." << endl;
        
    }
    os << "* " << endl;
    os << "* **********************************************************************************" << endl;

    return os;
}


void Distribution::setDistParametersGauss(double massIneV) {

    /*
     * Set distribution parameters. Do all the necessary checks depending
     * on the input attributes.
     * In case of DistributionType::MATCHEDGAUSS we only need to set the cutoff parameters
     */

    /*
    cutoffP_m = Vector_t<double, 3>(Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFPX]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFPY]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFPZ]));


    cutoffR_m = Vector_t<double, 3>(Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFX]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFY]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFLONG]));

    if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAR])) > 0.0) {
        cutoffR_m[0] = Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFR]);
        cutoffR_m[1] = Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFR]);
    }
    */
    
    setSigmaR_m();
    
}
void Distribution::createDistributionGauss(size_t numberOfParticles, double massIneV) {

    setDistParametersGauss(massIneV);

}




void Distribution::printDist(Inform &os, size_t numberOfParticles) const {

    

}


void Distribution::printDistGauss(Inform &os) const {
    os << "* Distribution type: GAUSS" << endl;
    os << "* " << endl;
    os << "* SIGMAX     = " << sigmaR_m[0] << " [m]" << endl;
    os << "* SIGMAY     = " << sigmaR_m[1] << " [m]" << endl;
    os << "* SIGMAZ     = " << sigmaR_m[2] << " [m]" << endl;
    os << "* SIGMAPX    = " << sigmaP_m[0] << " [Beta Gamma]" << endl;
    os << "* SIGMAPY    = " << sigmaP_m[1] << " [Beta Gamma]" << endl;
    os << "* SIGMAPZ    = " << sigmaP_m[2] << " [Beta Gamma]" << endl;
}

gsl_qrng* Distribution::selectRandomGenerator(std::string,unsigned int dimension)
{
    gsl_qrng *quasiRandGen = nullptr;
    if (Options::rngtype != std::string("RANDOM")) {
        INFOMSG("RNGTYPE= " << Options::rngtype << endl);
        if (Options::rngtype == std::string("HALTON")) {
            quasiRandGen = gsl_qrng_alloc(gsl_qrng_halton, dimension);
        } else if (Options::rngtype == std::string("SOBOL")) {
            quasiRandGen = gsl_qrng_alloc(gsl_qrng_sobol, dimension);
        } else if (Options::rngtype == std::string("NIEDERREITER")) {
            quasiRandGen = gsl_qrng_alloc(gsl_qrng_niederreiter_2, dimension);
        } else {
            INFOMSG("RNGTYPE= " << Options::rngtype << " not known, using HALTON" << endl);
            quasiRandGen = gsl_qrng_alloc(gsl_qrng_halton, dimension);
        }
    }
    return quasiRandGen;
}

void Distribution::setAttributes() {
}

void Distribution::setDistType() {

    static const std::map<std::string, DistributionType> typeStringToDistType_s = {
        {"NODIST",            DistributionType::NODIST},
        {"GAUSS",             DistributionType::GAUSS}
    };

    distT_m = Attributes::getString(itsAttr[DISTRIBUTION::TYPE]);

    if (distT_m.empty()) {
        throw OpalException("Distribution::setDistType",
                            "The attribute \"TYPE\" isn't set for the \"DISTRIBUTION\"!");
    } else {
        distrTypeT_m = typeStringToDistType_s.at(distT_m);
    }
}

void Distribution::setSigmaR_m() {
    sigmaR_m = Vector_t<double, 3>(std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAX])),
                        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAY])),
                        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAZ])));
 }

void Distribution::setSigmaP_m() {
    sigmaP_m = Vector_t<double, 3>(std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAPX])),
                        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAPY])),
                        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAPZ])));
}
