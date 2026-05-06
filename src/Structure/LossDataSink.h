//
// Class LossDataSink
//   This class writes file attributes to describe phase space of loss files
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef LOSSDATASINK_H_
#define LOSSDATASINK_H_

//////////////////////////////////////////////////////////////
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/OpalParticle.h"
 
#include "H5hut.h"

#include <fstream>
#include <optional>
#include <set>
#include <string>
#include <vector>

struct SetStatistics {
    SetStatistics();

    std::string outputName_m;
    double spos_m;
    double refTime_m;  // ns
    double tmean_m;    // ns
    double trms_m;     // ns
    double totalCharge_m;
    double totalMass_m;
    double meanKineticEnergy_m;
    double stdKineticEnergy_m;
    double meanGamma_m;
    unsigned long nTotal_m;
    Vector_t<double, 3> RefPartR_m;
    Vector_t<double, 3> RefPartP_m;
    Vector_t<double, 3> rmin_m;
    Vector_t<double, 3> rmax_m;
    Vector_t<double, 3> rmean_m;
    Vector_t<double, 3> pmean_m;
    Vector_t<double, 3> rrms_m;
    Vector_t<double, 3> prms_m;
    Vector_t<double, 3> rprms_m;
    Vector_t<double, 3> normEmit_m;
    Vector_t<double, 3> geomEmit_m;
    Vector_t<double, 3> maxR_m;
    Vector_t<double, 3> rsqsum_m;
    Vector_t<double, 3> psqsum_m;
    Vector_t<double, 3> rpsum_m;
    Vector_t<double, 3> eps2_m;
    Vector_t<double, 3> eps_norm_m;
    Vector_t<double, 3> fac_m;

    Vector_t<double, 3> sixtyEightPercentile_m;
    Vector_t<double, 3> ninetyFivePercentile_m;
    Vector_t<double, 3> ninetyNinePercentile_m;
    Vector_t<double, 3> ninetyNine_NinetyNinePercentile_m;

    Vector_t<double, 3> normalizedEps68Percentile_m;
    Vector_t<double, 3> normalizedEps95Percentile_m;
    Vector_t<double, 3> normalizedEps99Percentile_m;
    Vector_t<double, 3> normalizedEps99_99Percentile_m;
};

namespace std {
    template <>
    struct less<SetStatistics> {
        bool operator()(const SetStatistics& x, const SetStatistics& y) const {
            return x.spos_m < y.spos_m;
        }
    };
}  // namespace std
enum class CollectionType : unsigned short { SPATIAL = 0, TEMPORAL };

/*
  - In the destructor we do ALL the file handling
  - h5hut_mode_m defines h5hut or ASCII
 */
class LossDataSink {
public:
    LossDataSink() = default;

    LossDataSink(std::string outfn, bool hdf5Save, CollectionType = CollectionType::TEMPORAL);

    LossDataSink(const LossDataSink& rsh);
    ~LossDataSink() noexcept(false);

    bool inH5Mode() const {
        return h5hut_mode_m;
    }

    void save(
        unsigned int numSets = 1, OpalData::OpenMode openMode = OpalData::OpenMode::UNDEFINED);

    void addReferenceParticle(
        const Vector_t<double, 3>& x, const Vector_t<double, 3>& p, double time, double spos,
        long long globalTrackStep);

    void addParticle(
        const OpalParticle&,
        const std::optional<std::pair<int, short int>>& turnBunchNumPair = std::nullopt);

    size_t size() const;

    std::set<SetStatistics> computeStatistics(unsigned int numSets);

private:
    void writeFileAttribString(const char* attribute, const char* value);
    void writeStepAttribFloat64(const char* attribute, const h5_float64_t* value, h5_int64_t size);
    void writeStepAttribInt64(const char* attribute, const h5_int64_t* value, h5_int64_t size);
    void writeDataFloat64(const char* name, const h5_float64_t* value);
    void writeDataInt64(const char* name, const h5_int64_t* value);

    void setStep();
    void getNumSteps();
    void setNumParticles(h5_int64_t num);

    void openFile(const char* fname, h5_int32_t mode, h5_prop_t props);
    void closeFile();

    void openASCII() {
        if (ippl::Comm->rank() == 0) {
            os_m.open(fileName_m.c_str(), std::ios::out);
        }
    }
    void openH5(h5_int32_t mode = H5_O_WRONLY);

    void appendASCII() {
        if (ippl::Comm->rank() == 0) {
            os_m.open(fileName_m.c_str(), std::ios::app);
        }
    }

    void writeHeaderASCII();
    void writeHeaderH5();

    void saveASCII();
    void saveH5(unsigned int setIdx);

    void closeASCII() {
        if (ippl::Comm->rank() == 0) {
            os_m.close();
        }
    }

    bool hasNoParticlesToDump() const;
    bool hasTurnInformations() const;

    // void reportOnError(h5_int64_t rc, const char* file, int line);

    void splitSets(unsigned int numSets);
    SetStatistics computeSetStatistics(unsigned int setIdx);

    // filename without extension
    std::string fileName_m;

    // write either in ASCII or H5hut format
    bool h5hut_mode_m = false;

    // used to write out data in ASCII mode
    std::ofstream os_m;

    /// used to write out data in H5hut mode
    h5_file_t H5file_m = 0;

    std::string outputName_m;

    /// Current record, or time step, of H5 file.
    h5_int64_t H5call_m = 0;

    std::vector<OpalParticle> particles_m;
    std::vector<size_t> bunchNumber_m;
    std::vector<size_t> turnNumber_m;

    std::vector<Vector_t<double, 3>> RefPartR_m;
    std::vector<Vector_t<double, 3>> RefPartP_m;
    std::vector<h5_int64_t> globalTrackStep_m;
    std::vector<double> refTime_m;
    std::vector<double> spos_m;

    std::vector<unsigned long> startSet_m;

    CollectionType collectionType_m = CollectionType::TEMPORAL;
};

inline size_t LossDataSink::size() const {
    return particles_m.size();
}

inline std::set<SetStatistics> LossDataSink::computeStatistics(unsigned int numStatistics) {
    std::set<SetStatistics> stats;

    splitSets(numStatistics);

    for (unsigned int i = 0; i < numStatistics; ++i) {
        auto setStats = computeSetStatistics(i);
        if (setStats.nTotal_m > 0) {
            stats.insert(setStats);
        }
    }

    return stats;
}

#endif
