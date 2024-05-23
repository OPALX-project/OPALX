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
#ifndef OPAL_Distribution_HH
#define OPAL_Distribution_HH

#include "Ippl.h"

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Attributes/Attributes.h"

#include <memory>

#include "Algorithms/BoostMatrix.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Attributes/Attributes.h"
#include "Distribution/Distribution.h"
#include "Manager/BaseManager.h"
#include "Manager/PicManager.h"
#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/FieldSolver.hpp"
#include "PartBunch/LoadBalancer.hpp"
#include "PartBunch/ParticleContainer.hpp"

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

#include <fstream>
#include <string>
#include <vector>

class Beam;
class Beamline;
class H5PartWrapper;

enum class DistributionType : short { NODIST = -1, GAUSS };

using ParticleContainer_t = ParticleContainer<double, 3>;
using FieldContainer_t = FieldContainer<double, 3>;

class Distribution : public Definition {
public:
    Distribution();

    virtual ~Distribution();

    virtual bool canReplaceBy(Object* object);
    virtual Distribution* clone(const std::string& name);
    virtual void execute();
    virtual void update();

    void create(size_t &numberOfParticles, double massIneV, double charge, ippl::ParticleAttrib<ippl::Vector<double, 3>>& R, ippl::ParticleAttrib<ippl::Vector<double, 3>>& P, std::shared_ptr<ParticleContainer<double, 3>> &pc, std::shared_ptr<FieldContainer_t> &fc, Vector_t<double, 3> nr);

    size_t getNumOfLocalParticlesToCreate(size_t n);

    // void createOpalT(
    //     PartBunch_t* beam, std::vector<Distribution*> addedDistributions,
    //     size_t& numberOfParticles);

    // void createOpalT(PartBunch_t* beam, size_t& numberOfParticles);

    static Distribution* find(const std::string& name);

    std::string getTypeofDistribution();
    DistributionType getType() const;

    Inform& printInfo(Inform& os) const;

    ippl::Vector<double, 3> get_pmean() const;
    ippl::Vector<double, 3> get_xmean() const;

    ippl::Vector<double, 3> getSigmaR() const;
    ippl::Vector<double, 3> getSigmaP() const;

    void setDistType();

    void setDist();

    void setAvrgPz(double avrgpz);

    double getAvrgpz() const;

    ippl::Vector<double, 3> getCutoffR() const;
    ippl::Vector<double, 3> getCutoffP() const;

private:
    enum class EmissionModel : unsigned short { NONE, ASTRA, NONEQUIL };

    enum class InputMomentumUnits : unsigned short { NONE, EVOVERC };

#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(GaussTest, FullSigmaTest1);
    FRIEND_TEST(GaussTest, FullSigmaTest2);
#endif

    Distribution(const std::string& name, Distribution* parent);

    // Not implemented.
    Distribution(const Distribution&) = delete;

    void operator=(const Distribution&) = delete;

    void addDistributions();

    /*!
     * Create the particle distribution.
     * @param numberOfParticles to create
     * @param massIneV particle charge in eV
     * @param charge of the particle type in elementary charge
     */
    //void create(size_t& numberOfParticles, double massIneV, double charge, auto R);
    void calcPartPerDist(size_t numberOfParticles);
    void checkParticleNumber(size_t& numberOfParticles);

    void chooseInputMomentumUnits(InputMomentumUnits inputMoUnits);
    size_t getNumberOfParticlesInFile(std::ifstream& inputFile);

    class BinomialBehaviorSplitter {
    public:
        virtual ~BinomialBehaviorSplitter() {
        }

        virtual double get(double rand) = 0;
    };

    class MDependentBehavior : public BinomialBehaviorSplitter {
    public:
        MDependentBehavior(const MDependentBehavior& rhs) : ami_m(rhs.ami_m) {
        }

        MDependentBehavior(double a) {
            ami_m = 1.0 / a;
        }

        virtual double get(double rand);

    private:
        double ami_m;
    };

    class GaussianLikeBehavior : public BinomialBehaviorSplitter {
    public:
        virtual double get(double rand);
    };

    void createDistributionGauss(size_t numberOfParticles, double massIneV, ippl::ParticleAttrib<ippl::Vector<double, 3>>& R, ippl::ParticleAttrib<ippl::Vector<double, 3>>& P, std::shared_ptr<ParticleContainer<double, 3>> &pc, std::shared_ptr<FieldContainer_t> &fc, Vector_t<double, 3> nr);

    //
    //
    // void initializeBeam(PartBunch_t* beam);
    void printDist(Inform& os, size_t numberOfParticles) const;
    void printDistGauss(Inform& os) const;

    void setAttributes();

    void setDistParametersGauss();

    void setSigmaR_m();

    void setSigmaP_m();

    /*
      private member of Distribution
    */

    std::string distT_m;  /// Distribution type strings.

    PartData particleRefData_m;  /// Reference data for particle type (charge,
                                 /// mass etc.)

    size_t totalNumberParticles_m;

    ippl::Vector<double, 3> pmean_m, xmean_m, sigmaR_m, sigmaP_m, cutoffR_m, cutoffP_m;

    DistributionType distrTypeT_m;

    double avrgpz_m;
};

inline Inform& operator<<(Inform& os, const Distribution& d) {
    return d.printInfo(os);
}

inline ippl::Vector<double, 3> Distribution::getSigmaR() const {
    return sigmaR_m;
}

inline ippl::Vector<double, 3> Distribution::getSigmaP() const {
    return sigmaP_m;
}

inline ippl::Vector<double, 3> Distribution::get_pmean() const {
    return pmean_m;
}

inline ippl::Vector<double, 3> Distribution::get_xmean() const {
    return xmean_m;
}

inline ippl::Vector<double, 3> Distribution::getCutoffR() const {
    return cutoffR_m;
}

inline ippl::Vector<double, 3> Distribution::getCutoffP() const {
    return cutoffP_m;
}

inline DistributionType Distribution::getType() const {
    return distrTypeT_m;
}

inline double Distribution::getAvrgpz() const {
    return avrgpz_m;
}

inline std::string Distribution::getTypeofDistribution() {
    return distT_m;
}

#endif  
/*
// OPAL_Distribution_HH
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
#ifndef OPAL_Distribution_HH
#define OPAL_Distribution_HH

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"

#include "Attributes/Attributes.h"

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

#include <fstream>
#include <string>
#include <vector>

class Beam;
class Beamline;

template <class T, unsigned Dim>
class PartBunch;

class H5PartWrapper;

enum class DistributionType : short { NODIST = -1, GAUSS };

namespace DISTRIBUTION {
    enum { TYPE, FNAME, SIGMAX, SIGMAY, SIGMAZ, SIGMAPX, SIGMAPY, SIGMAPZ, SIZE, CUTOFFPX, CUTOFFPY, CUTOFFPZ, CUTOFFR };
}

class Distribution : public Definition {
public:
    Distribution();

    virtual ~Distribution();

    virtual bool canReplaceBy(Object* object);
    virtual Distribution* clone(const std::string& name);
    virtual void execute();
    virtual void update();

    // void create(size_t &numberOfParticles, double massIneV, double charge);

    size_t getNumOfLocalParticlesToCreate(size_t n);

    void createOpalT(
        PartBunch_t* beam, std::vector<Distribution*> addedDistributions,
        size_t& numberOfParticles);
    void createOpalT(PartBunch_t* beam, size_t& numberOfParticles);

    static Distribution* find(const std::string& name);

    std::string getTypeofDistribution();
    DistributionType getType() const;

    Inform& printInfo(Inform& os) const;

    ippl::Vector<double, 3> get_pmean() const;
    ippl::Vector<double, 3> get_xmean() const;

    void setDistType();

    void setAvrgPz(double avrgpz);

private:
    enum class EmissionModel : unsigned short { NONE, ASTRA, NONEQUIL };

    enum class InputMomentumUnits : unsigned short { NONE, EVOVERC };

#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(GaussTest, FullSigmaTest1);
    FRIEND_TEST(GaussTest, FullSigmaTest2);
#endif

    Distribution(const std::string& name, Distribution* parent);

    // Not implemented.
    Distribution(const Distribution&) = delete;

    void operator=(const Distribution&) = delete;

    void addDistributions();

    ///!
     // Create the particle distribution.
     // @param numberOfParticles to create
     // @param massIneV particle charge in eV
     // @param charge of the particle type in elementary charge
     ///
    void create(size_t& numberOfParticles, double massIneV, double charge);
    void calcPartPerDist(size_t numberOfParticles);
    void checkParticleNumber(size_t& numberOfParticles);

    void chooseInputMomentumUnits(InputMomentumUnits inputMoUnits);
    size_t getNumberOfParticlesInFile(std::ifstream& inputFile);

    class BinomialBehaviorSplitter {
    public:
        virtual ~BinomialBehaviorSplitter() {
        }

        virtual double get(double rand) = 0;
    };

    class MDependentBehavior : public BinomialBehaviorSplitter {
    public:
        MDependentBehavior(const MDependentBehavior& rhs) : ami_m(rhs.ami_m) {
        }

        MDependentBehavior(double a) {
            ami_m = 1.0 / a;
        }

        virtual double get(double rand);

    private:
        double ami_m;
    };

    class GaussianLikeBehavior : public BinomialBehaviorSplitter {
    public:
        virtual double get(double rand);
    };

    void createDistributionGauss(size_t numberOfParticles, double massIneV);

    void initializeBeam(PartBunch_t* beam);
    void printDist(Inform& os, size_t numberOfParticles) const;
    void printDistGauss(Inform& os) const;

    gsl_qrng* selectRandomGenerator(std::string, unsigned int dimension);

    void setAttributes();

    void setDistParametersGauss(double massIneV);

    void setSigmaR_m();

    void setSigmaP_m();


    
      //private member of Distribution
    

    std::string distT_m;  /// Distribution type strings.

    PartData particleRefData_m;  /// Reference data for particle type (charge,
                                 /// mass etc.)

    gsl_rng* randGen_m;  /// Random number generator

    size_t totalNumberParticles_m;

    ippl::Vector<double, 3> pmean_m, xmean_m, sigmaR_m, sigmaP_m;

    DistributionType distrTypeT_m;

    //double avrgpz_m;
};

inline Inform& operator<<(Inform& os, const Distribution& d) {
    return d.printInfo(os);
}

inline ippl::Vector<double, 3> Distribution::get_pmean() const {
    return pmean_m;
}

inline ippl::Vector<double, 3> Distribution::get_xmean() const {
    return xmean_m;
}


inline DistributionType Distribution::getType() const {
    return distrTypeT_m;
}

inline std::string Distribution::getTypeofDistribution() {
    return distT_m;
}

#endif  // OPAL_Distribution_HH
*/
