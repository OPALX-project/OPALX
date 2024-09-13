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
#include "BasicActions/Option.h"
#include "Distribution/LaserProfile.h"
#include "Elements/OpalBeamline.h"
#include "OPALTypes.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/H5PartWrapper.h"
#include "Utilities/EarlyLeaveException.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Utility/IpplTimings.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_erf.h>

#include <boost/filesystem.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/regex.hpp>

#include <sys/time.h>

#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>

extern Inform* gmsg;

using GeneratorPool = typename Kokkos::Random_XorShift64_Pool<>;

using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, Dim>, 1>::view_type;

namespace DISTRIBUTION {
    enum { TYPE, FNAME, SIGMAX, SIGMAY, SIGMAZ, SIGMAPX, SIGMAPY, SIGMAPZ, CORR,
           CUTOFFPX, CUTOFFPY, CUTOFFPZ, CUTOFFX, CUTOFFY, CUTOFFLONG, CORRX, CORRY,
           CORRZ, CORRT, SIGMAT, TPULSEFWHM, TRISE, TFALL, FTOSCAMPLITUDE, FTOSCPERIODS, SIZE };
}

/*
namespace {
    matrix_t getUnit6x6() {
        matrix_t unit6x6(6, 6, 0.0);  // Initialize a 6x6 matrix with all elements as 0.0
        for (unsigned int i = 0; i < 6u; ++i) {
            unit6x6(i, i) = 1.0;  // Set diagonal elements to 1.0
        }
        return unit6x6;
    }
}
*/

Distribution::Distribution()
    : Definition(
        DISTRIBUTION::SIZE, "DISTRIBUTION",
        "The DISTRIBUTION statement defines data for the 6D particle distribution."),
      distrTypeT_m(DistributionType::NODIST) {
    itsAttr[DISTRIBUTION::TYPE] =
        Attributes::makePredefinedString("TYPE", "Distribution type.", {"GAUSS", "MULTIVARIATEGAUSS", "FLATTOP", "FROMFILE"});

    itsAttr[DISTRIBUTION::FNAME] =
        Attributes::makeString("FNAME", "File for reading in 6D particle coordinates.", "");

    // Parameters for defining an initial distribution.
    itsAttr[DISTRIBUTION::SIGMAX]  = Attributes::makeReal("SIGMAX", "SIGMAx (m)", 0.0);
    itsAttr[DISTRIBUTION::SIGMAY]  = Attributes::makeReal("SIGMAY", "SIGMAy (m)", 0.0);
    itsAttr[DISTRIBUTION::SIGMAZ]  = Attributes::makeReal("SIGMAZ", "SIGMAz (m)", 0.0);
    itsAttr[DISTRIBUTION::SIGMAPX] = Attributes::makeReal("SIGMAPX", "SIGMApx", 0.0);
    itsAttr[DISTRIBUTION::SIGMAPY] = Attributes::makeReal("SIGMAPY", "SIGMApy", 0.0);
    itsAttr[DISTRIBUTION::SIGMAPZ] = Attributes::makeReal("SIGMAPZ", "SIGMApz", 0.0);

    itsAttr[DISTRIBUTION::CORR] = Attributes::makeRealArray("CORR", "r correlation");

    itsAttr[DISTRIBUTION::CUTOFFPX] = Attributes::makeReal("CUTOFFPX", "Distribution cutoff px dimension in units of sigma.", 3.0);
    itsAttr[DISTRIBUTION::CUTOFFPY] = Attributes::makeReal("CUTOFFPY", "Distribution cutoff py dimension in units of sigma.", 3.0);
    itsAttr[DISTRIBUTION::CUTOFFPZ] = Attributes::makeReal("CUTOFFPZ", "Distribution cutoff pz dimension in units of sigma.", 3.0);

    itsAttr[DISTRIBUTION::CUTOFFX] = Attributes::makeReal("CUTOFFX", "Distribution cutoff x direction in units of sigma.", 3.0);
    itsAttr[DISTRIBUTION::CUTOFFY] = Attributes::makeReal("CUTOFFY", "Distribution cutoff r direction in units of sigma.", 3.0);
    itsAttr[DISTRIBUTION::CUTOFFLONG] = Attributes::makeReal("CUTOFFLONG", "Distribution cutoff z or t direction in units of sigma.", 3.0);

    itsAttr[DISTRIBUTION::CORRX] = Attributes::makeReal("CORRX", "x/px correlation, (R12 in transport notation).", 0.0);
    itsAttr[DISTRIBUTION::CORRY] = Attributes::makeReal("CORRY", "y/py correlation, (R34 in transport notation).", 0.0);
    itsAttr[DISTRIBUTION::CORRZ] = Attributes::makeReal("CORRZ", "z/pz correlation, (R56 in transport notation).", 0.0);
    itsAttr[DISTRIBUTION::CORRT] = Attributes::makeReal("CORRT", "t/pt correlation, (R56 in transport notation).", 0.0);

    itsAttr[DISTRIBUTION::SIGMAT] = Attributes::makeReal("SIGMAT", "SIGMAt (m)", 0.0);
    itsAttr[DISTRIBUTION::TPULSEFWHM] = Attributes::makeReal("TPULSEFWHM", "Pulse FWHM for emitted distribution.", 0.0);
    itsAttr[DISTRIBUTION::TRISE] = Attributes::makeReal("TRISE", "Rise time for emitted distribution.", 0.0);
    itsAttr[DISTRIBUTION::TFALL] = Attributes::makeReal("TFALL", "Fall time for emitted distribution.", 0.0);

    itsAttr[DISTRIBUTION::FTOSCAMPLITUDE]
        = Attributes::makeReal("FTOSCAMPLITUDE", "Amplitude of oscillations superimposed "
                               "on flat top portion of emitted GAUSS "
                               "distribtuion (in percent of flat top "
                               "amplitude)",0.0);

    itsAttr[DISTRIBUTION::FTOSCPERIODS]
        = Attributes::makeReal("FTOSCPERIODS", "Number of oscillations superimposed on "
                               "flat top portion of emitted GAUSS "
                               "distribution", 0.0);

    registerOwnership(AttributeHandler::STATEMENT);
}

Distribution::Distribution(const std::string& name, Distribution* parent)
    : Definition(name, parent) {
}

Distribution::~Distribution() {
}

/**
 * Calculate the local number of particles evenly and adjust node 0
 * such that n is matched exactly.
 * @param n total number of particles
 * @return n / #cores
 * @param
 */

size_t Distribution::getNumOfLocalParticlesToCreate(size_t n) {
    size_t locNumber = n / ippl::Comm->size();

    // make sure the total number is exact
    size_t remainder = n % ippl::Comm->size();
    size_t myNode    = ippl::Comm->rank();
    if (myNode < remainder)
        ++locNumber;

    return locNumber;
}

/// Distribution can only be replaced by another distribution.
bool Distribution::canReplaceBy(Object* object) {
    return dynamic_cast<Distribution*>(object) != 0;
}

Distribution* Distribution::clone(const std::string& name) {
    return new Distribution(name, this);
}

void Distribution::execute() {
    setAttributes();
    update();
}

void Distribution::update() {
}

void Distribution::create(size_t& numberOfParticles, double massIneV, double charge, ippl::ParticleAttrib<ippl::Vector<double, 3>>& R, ippl::ParticleAttrib<ippl::Vector<double, 3>>& P, std::shared_ptr<ParticleContainer_t> &pc, std::shared_ptr<FieldContainer_t> &fc, Vector_t<double, 3> nr) {
// moved to SamplingBase
}

Distribution* Distribution::find(const std::string& name) {
    Distribution* dist = dynamic_cast<Distribution*>(OpalData::getInstance()->find(name));

    if (dist == 0) {
        throw OpalException("Distribution::find()", "Distribution \"" + name + "\" not found.");
    }

    return dist;
}

Inform& Distribution::printInfo(Inform& os) const {
    os << "\n"
       << "* ************* D I S T R I B U T I O N ********************************************"
       << endl;
    os << "* " << endl;
    if (OpalData::getInstance()->inRestartRun()) {
        os << "* In restart. Distribution read in from .h5 file." << endl;
    } else {
        switch (distrTypeT_m) {
            case DistributionType::GAUSS:
                printDistGauss(os);
                break;
            case DistributionType::MULTIVARIATEGAUSS:
                printDistMultiVariateGauss(os);
                break;
            case DistributionType::FLATTOP:
                printDistFlatTop(os);
                break;
            default:
                throw OpalException("Distribution Param", "Unknown \"TYPE\" of \"DISTRIBUTION\"");
         }
        os << "* " << endl;
        os << "* Distribution is injected." << endl;
    }
    os << "* " << endl;
    os << "* **********************************************************************************"
       << endl;

    return os;
}

void Distribution::setAvrgPz(double avrgpz){
    avrgpz_m = avrgpz;
}

void Distribution::setDistParametersGauss() {
    /*
     * Set distribution parameters. Do all the necessary checks depending
     * on the input attributes.
     * In case of DistributionType::MATCHEDGAUSS we only need to set the cutoff parameters
     */

    cutoffR_m = 3.;
    cutoffP_m = 3.;
    /*
    cutoffP_m = ippl::Vector<double, 3>(Attributes::getReal(itsAttr[DISTRIBUTION::CUTOFFPX]),
                         Attributes::getReal(itsAttr[DISTRIBUTION::CUTOFFPY]),
                         Attributes::getReal(itsAttr[DISTRIBUTION::CUTOFFPZ]));


    cutoffR_m = ippl::Vector<double, 3>(Attributes::getReal(itsAttr[DISTRIBUTION::CUTOFFX]),
                         Attributes::getReal(itsAttr[DISTRIBUTION::CUTOFFY]),
                         Attributes::getReal(itsAttr[DISTRIBUTION::CUTOFFLONG]));
    */

    //if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAR])) > 0.0) {
    //    cutoffR_m[0] = Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFR]);
    //    cutoffR_m[1] = Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFR]);
    //}

    setSigmaR_m();
    setSigmaP_m();

    avrgpz_m = 0.0;
}

void Distribution::setDistParametersMultiVariateGauss() {

    cutoffR_m = 3.;
    cutoffP_m = 3.;

    // initialize the covariance matrix to identity
    for (unsigned int i = 0; i < 6; ++ i) {
        for (unsigned int j = 0; j < 6; ++ j) {
            if (i==j)
               correlationMatrix_m[i][j] = 1.0;
            else
               correlationMatrix_m[i][j] = 0.0;
        }
    }

    // set diagonal elements first
    setSigmaR_m();
    setSigmaP_m();

    for (unsigned int i = 0; i < 3; ++ i){
        correlationMatrix_m[2*i  ][2*i  ] = sigmaR_m[i]*sigmaR_m[i];
        correlationMatrix_m[2*i+1][2*i+1] = sigmaP_m[i]*sigmaP_m[i];
    }

    std::vector<double> cr = Attributes::getRealArray(itsAttr[DISTRIBUTION::CORR]);

    if (!cr.empty()) {
            // read off-diagonal correlation matrix from input file
            if (cr.size() == 15) {
                *gmsg << "* Use r to specify correlations" << endl;
                unsigned int k = 0;
                for (unsigned int i = 0; i < 5; ++ i) {
                    for (unsigned int j = i + 1; j < 6; ++ j, ++ k) {
                        correlationMatrix_m[j][i] = cr.at(k)*cr.at(k);
                        correlationMatrix_m[i][j] = correlationMatrix_m[j][i]; // impose symmetry
                    }
                }
            }
            else {
                throw OpalException("Distribution::SetDistParametersGauss",
                                    "Inconsistent set of correlations specified, check manual");
            }
    }

    avrgpz_m = 0.0;
}

void Distribution::setDistParametersFlatTop() {

    cutoffR_m = 3.;
    cutoffP_m = 3.;

    // set diagonal elements first
    setSigmaR_m();
    setSigmaP_m();

    // initialize the covariance matrix to identity
    for (unsigned int i = 0; i < 6; ++ i) {
        for (unsigned int j = 0; j < 6; ++ j) {
            if (i==j)
               correlationMatrix_m[i][j] = 1.0;
            else
               correlationMatrix_m[i][j] = 0.0;
        }
    }

    cutoffR_m = ippl::Vector<double, 3>(
        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::CUTOFFX])),
        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::CUTOFFY])),
        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::CUTOFFLONG])));

    correlationMatrix_m[1][0] = Attributes::getReal(itsAttr[DISTRIBUTION::CORRX]);
    correlationMatrix_m[3][2] = Attributes::getReal(itsAttr[DISTRIBUTION::CORRY]);
    correlationMatrix_m[5][4] = Attributes::getReal(itsAttr[DISTRIBUTION::CORRT]);

    if (Attributes::getReal(itsAttr[DISTRIBUTION::CORRZ]) != 0.0)
        correlationMatrix_m[5][4] = Attributes::getReal(itsAttr[DISTRIBUTION::CORRZ]);

    if (emitting_m) {
        sigmaR_m[2] = 0.0;

        sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAT]));
        sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAT]));
        tPulseLengthFWHM_m = std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::TPULSEFWHM]));

        FTOSCAmplitude_m = std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::FTOSCAMPLITUDE]));
        FTOSCPeriods_m = std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::FTOSCPERIODS]));

        // If TRISE and TFALL are defined > 0.0 then these attributes
        // override SIGMAT.
        //
        if (std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::TRISE])) > 0.0
            || std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::TFALL])) > 0.0) {

            double timeRatio = std::sqrt(2.0 * std::log(10.0)) - std::sqrt(2.0 * std::log(10.0 / 9.0));

            if (std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::TRISE])) > 0.0)
                sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::TRISE]))
                    / timeRatio;

            if (std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::TFALL])) > 0.0)
                sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::TFALL]))
                    / timeRatio;
        }

        // For an emitted beam, the longitudinal cutoff >= 0.
        cutoffR_m[2] = std::abs(cutoffR_m[2]);

    }
    /*
    cutoffR_m = Vector_t(Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFX]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFY]),
                         Attributes::getReal(itsAttr[Attrib::Distribution::CUTOFFLONG]));

    correlationMatrix_m(1, 0) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRX]);
    correlationMatrix_m(3, 2) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRY]);
    correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRT]);

    // CORRZ overrides CORRT.
    if (Attributes::getReal(itsAttr[Attrib::Distribution::CORRZ]) != 0.0)
        correlationMatrix_m(5, 4) = Attributes::getReal(itsAttr[Attrib::Distribution::CORRZ]);

    setSigmaR_m();
    if (emitting_m) {
        sigmaR_m[2] = 0.0;

        sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAT]));
        sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::SIGMAT]));

        tPulseLengthFWHM_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TPULSEFWHM]));

        //
         // If TRISE and TFALL are defined > 0.0 then these attributes
         // override SIGMAT.
         //
        if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE])) > 0.0
            || std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TFALL])) > 0.0) {

            double timeRatio = std::sqrt(2.0 * std::log(10.0)) - std::sqrt(2.0 * std::log(10.0 / 9.0));

            if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE])) > 0.0)
                sigmaTRise_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE]))
                    / timeRatio;

            if (std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TFALL])) > 0.0)
                sigmaTFall_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TFALL]))
                    / timeRatio;

        }

        // For an emitted beam, the longitudinal cutoff >= 0.
        cutoffR_m[2] = std::abs(cutoffR_m[2]);
    }

    // Set laser profile/
    laserProfileFileName_m = Attributes::getString(itsAttr[Attrib::Distribution::LASERPROFFN]);
    if (!(laserProfileFileName_m == std::string(""))) {
        laserImageName_m = Attributes::getString(itsAttr[Attrib::Distribution::IMAGENAME]);
        laserIntensityCut_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::INTENSITYCUT]));
        short flags = 0;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::FLIPX])) flags |= LaserProfile::FLIPX;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::FLIPY])) flags |= LaserProfile::FLIPY;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::ROTATE90])) flags |= LaserProfile::ROTATE90;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::ROTATE180])) flags |= LaserProfile::ROTATE180;
        if (Attributes::getBool(itsAttr[Attrib::Distribution::ROTATE270])) flags |= LaserProfile::ROTATE270;

        laserProfile_m = new LaserProfile(laserProfileFileName_m,
                                          laserImageName_m,
                                          laserIntensityCut_m,
                                          flags);
    }

    // Legacy for ASTRAFLATTOPTH.
    if (distrTypeT_m == DistributionType::ASTRAFLATTOPTH)
        tRise_m = std::abs(Attributes::getReal(itsAttr[Attrib::Distribution::TRISE]));
*/

}

void Distribution::createDistributionGauss(size_t numberOfParticles, double massIneV, ippl::ParticleAttrib<ippl::Vector<double, 3>>& R, ippl::ParticleAttrib<ippl::Vector<double, 3>>& P, std::shared_ptr<ParticleContainer_t> &pc, std::shared_ptr<FieldContainer_t> &fc, Vector_t<double, 3> nr) {
    // moved to Gaussian.hpp
}

void Distribution::printDist(Inform& os, size_t numberOfParticles) const {
}

void Distribution::printDistGauss(Inform& os) const {
    os << "* Distribution type: GAUSS" << endl;
    os << "* " << endl;
    os << "* SIGMAX     = " << sigmaR_m[0] << " [m]" << endl;
    os << "* SIGMAY     = " << sigmaR_m[1] << " [m]" << endl;
    os << "* SIGMAZ     = " << sigmaR_m[2] << " [m]" << endl;
    os << "* SIGMAPX    = " << sigmaP_m[0] << " [Beta Gamma]" << endl;
    os << "* SIGMAPY    = " << sigmaP_m[1] << " [Beta Gamma]" << endl;
    os << "* SIGMAPZ    = " << sigmaP_m[2] << " [Beta Gamma]" << endl;
}

void Distribution::printDistMultiVariateGauss(Inform& os)  const {
    os << "* Distribution type: MULTIVARIATEGAUSS" << endl;
    os << "* " << endl;
    os << "* SIGMAX     = " << sigmaR_m[0] << " [m]" << endl;
    os << "* SIGMAY     = " << sigmaR_m[1] << " [m]" << endl;
    os << "* SIGMAZ     = " << sigmaR_m[2] << " [m]" << endl;
    os << "* SIGMAPX    = " << sigmaP_m[0] << " [Beta Gamma]" << endl;
    os << "* SIGMAPY    = " << sigmaP_m[1] << " [Beta Gamma]" << endl;
    os << "* SIGMAPZ    = " << sigmaP_m[2] << " [Beta Gamma]" << endl;

    os << "* input cov matrix = ";
    for (unsigned int i = 0; i < 6; ++ i) {
        for (unsigned int j = 0; j < 6; ++ j) {
            os << correlationMatrix_m[i][j] << " ";
        }
        os << endl << "                     ";
    }
}

void Distribution::printDistFlatTop(Inform& os)  const {
    os << "* Distribution type: FLATTOP" << endl;
    os << "* " << endl;
    os << "* SIGMAX     = " << sigmaR_m[0] << " [m]" << endl;
    os << "* SIGMAY     = " << sigmaR_m[1] << " [m]" << endl;

    if (emitting_m) {
            os << "* Sigma Time Rise               = " << sigmaTRise_m
               << " [sec]" << endl;
            os << "* TPULSEFWHM                    = " << tPulseLengthFWHM_m
               << " [sec]" << endl;
            os << "* Sigma Time Fall               = " << sigmaTFall_m
               << " [sec]" << endl;
            os << "* Longitudinal cutoff           = " << cutoffR_m[2]
               << " [units of Sigma Time]" << endl;
            //os << "* Flat top modulation amplitude = "
            //   << Attributes::getReal(itsAttr[DISTRIBUTION::FTOSCAMPLITUDE])
            //   << " [Percent of distribution amplitude]" << endl;
            //os << "* Flat top modulation periods   = "
            //   << std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::FTOSCPERIODS]))
            //   << endl;
    }
    else{
        os << "* SIGMAZ                        = " << sigmaR_m[2] << " [m]" << endl;
    }
}

void Distribution::setAttributes() {
    setDist();
}

void Distribution::setDist() {
    // set distribution type
    setDistType();
    // set distribution parameters
    switch (distrTypeT_m) {
        case DistributionType::GAUSS:
            setDistParametersGauss();
            break;
        case DistributionType::MULTIVARIATEGAUSS:
            setDistParametersMultiVariateGauss();
            break;
        case DistributionType::FLATTOP:
            setDistParametersFlatTop();
            break;
        default:
            throw OpalException("Distribution Param", "Unknown \"TYPE\" of \"DISTRIBUTION\"");
    }
}

void Distribution::setDistType() {
    static const std::map<std::string, DistributionType> typeStringToDistType_s = {
        {"NODIST", DistributionType::NODIST},
        {"GAUSS", DistributionType::GAUSS},
        {"MULTIVARIATEGAUSS", DistributionType::MULTIVARIATEGAUSS},
        {"FLATTOP", DistributionType::FLATTOP}
    };

    distT_m = Attributes::getString(itsAttr[DISTRIBUTION::TYPE]);

    if (distT_m.empty()) {
        throw OpalException(
            "Distribution::setDistType",
            "The attribute \"TYPE\" isn't set for the \"DISTRIBUTION\"!");
    } else {
        distrTypeT_m = typeStringToDistType_s.at(distT_m);
    }
}

void Distribution::setSigmaR_m() {
    sigmaR_m = ippl::Vector<double, 3>(
        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAX])),
        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAY])),
        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAZ])));
}

void Distribution::setSigmaP_m() {
    sigmaP_m = ippl::Vector<double, 3>(
        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAPX])),
        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAPY])),
        std::abs(Attributes::getReal(itsAttr[DISTRIBUTION::SIGMAPZ])));
}

