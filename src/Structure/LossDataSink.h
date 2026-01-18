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
#include "Ippl.h"
 

#include <boost/optional.hpp>

#include <fstream>
#include <functional>
#include <set>
#include <string>
#include <vector>
#include <cstdint>

#include "H5hut.h"
#include <Kokkos_Core.hpp>

struct SetStatistics {
    SetStatistics();

    std::string outputName_m;
    double spos_m;
    double refTime_m;  // ns
    double tmean_m;    // ns
    double trms_m;     // ns
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
    Vector_t<double, 3> rsqsum_m;
    Vector_t<double, 3> psqsum_m;
    Vector_t<double, 3> rpsum_m;
    Vector_t<double, 3> eps2_m;
    Vector_t<double, 3> eps_norm_m;
    Vector_t<double, 3> fac_m;
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
    struct DeviceData {
        using RView_t = ippl::ParticleAttrib<Vector_t<double, 3>>::view_type;
        using ScalarView_t = ippl::ParticleAttrib<double>::view_type;
        using IdView_t = Kokkos::View<int64_t*>;
        using IndexView_t = Kokkos::View<size_t*>;
        using TimeView_t = Kokkos::View<double*>;
        using TurnView_t = Kokkos::View<int*>;
        using BunchView_t = Kokkos::View<short int*>;

        RView_t R;
        RView_t P;
        ScalarView_t Q;
        ScalarView_t M;
        TimeView_t Time;
        IdView_t ID;

        IndexView_t lossIndex;
        size_t localCount = 0;
        size_t lossCount = 0;

        TurnView_t turn;
        BunchView_t bunch;

        bool hasLossIndex() const {
            return lossCount > 0 || lossIndex.extent(0) > 0;
        }
        bool hasTurnBunch() const {
            return turn.extent(0) > 0 && bunch.extent(0) > 0;
        }
        bool hasTime() const {
            return Time.extent(0) > 0;
        }
    };

    LossDataSink() = default;

    LossDataSink(std::string outfn, bool hdf5Save, CollectionType = CollectionType::TEMPORAL);

    LossDataSink(const LossDataSink& rsh);
    ~LossDataSink() noexcept(false);

    bool inH5Mode() {
        return h5hut_mode_m;
    }

    void save(
        unsigned int numSets = 1, OpalData::OpenMode openMode = OpalData::OpenMode::UNDEFINED);

    void addReferenceParticle(
        const Vector_t<double, 3>& x, const Vector_t<double, 3>& p, double time, double spos,
        long long globalTrackStep);

    void addParticle(
        const OpalParticle&,
        const boost::optional<std::pair<int, short int>>& turnBunchNumPair = boost::none);

    void bindDeviceData(const DeviceData& data);
    void clearDeviceData();
    bool hasDeviceData() const {
        return deviceData_m.is_initialized();
    }

    size_t size() const;

    std::set<SetStatistics> computeStatistics(unsigned int numSets);

private:
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

    void reportOnError(h5_int64_t rc, const char* file, int line);

    void splitSets(unsigned int numSets);
    SetStatistics computeSetStatistics(unsigned int setIdx);
    std::set<SetStatistics> computeStatisticsFromHost(unsigned int numSets);

    // filename without extension
    std::string fileName_m;

    // write either in ASCII or H5hut format
    bool h5hut_mode_m;

    // used to write out data in ASCII mode
    std::ofstream os_m;

    /// used to write out data in H5hut mode
    h5_file_t H5file_m;

    std::string outputName_m;

    /// Current record, or time step, of H5 file.
    h5_int64_t H5call_m;

    std::vector<OpalParticle> particles_m;
    std::vector<size_t> bunchNumber_m;
    std::vector<size_t> turnNumber_m;

    std::vector<Vector_t<double, 3>> RefPartR_m;
    std::vector<Vector_t<double, 3>> RefPartP_m;
    std::vector<h5_int64_t> globalTrackStep_m;
    std::vector<double> refTime_m;
    std::vector<double> spos_m;

    std::vector<unsigned long> startSet_m;

    CollectionType collectionType_m;

    boost::optional<DeviceData> deviceData_m;
};

inline size_t LossDataSink::size() const {
    if (!particles_m.empty()) {
    return particles_m.size();
}
    if (deviceData_m) {
        if (deviceData_m->lossCount > 0) {
            return deviceData_m->lossCount;
        }
        if (deviceData_m->localCount > 0) {
            return deviceData_m->localCount;
        }
        return deviceData_m->R.extent(0);
        }
    return particles_m.size();
}

#endif
