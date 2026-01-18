//
// Class Monitor
//   Defines the abstract interface for a beam position monitor.
//
// Copyright (c) 2000 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#include "AbsBeamline/Monitor.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "AbstractObjects/OpalData.h"
#include "Fields/Fieldmap.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Structure/MonitorStatisticsWriter.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>

#include <fstream>
#include <memory>

std::map<double, SetStatistics> Monitor::statFileEntries_sm;
const double Monitor::halfLength_s = 0.005;

Monitor::Monitor() : Monitor("") {
}

Monitor::Monitor(const Monitor& right)
    : Component(right),
      filename_m(right.filename_m),
      plane_m(right.plane_m),
      type_m(right.type_m),
      numPassages_m(0) {
}

Monitor::Monitor(const std::string& name)
    : Component(name),
      filename_m(""),
      plane_m(OFF),
      type_m(CollectionType::SPATIAL),
      numPassages_m(0) {
}

Monitor::~Monitor() {
}

void Monitor::accept(BeamlineVisitor& visitor) const {
    visitor.visitMonitor(*this);
}

bool Monitor::apply(
    const size_t& i, const double& t, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {
    // const Vector_t<double, 3>& R         = RefPartBunch_m->R(i);
    // const Vector_t<double, 3>& P         = RefPartBunch_m->P(i);
    // const double& dt                     = RefPartBunch_m->dt(i);
    // const Vector_t<double, 3> singleStep = Physics::c * dt * Util::getBeta(P);
    // if (online_m && type_m == CollectionType::SPATIAL) {
    //     if (dt * R(2) < 0.0 && dt * (R(2) + singleStep(2)) > 0.0) {
    //         // if R(2) is negative then frac should be positive and vice versa
    //         double frac = -R(2) / singleStep(2);

    //         lossDs_m->addParticle(OpalParticle(
    //             RefPartBunch_m->ID[i], R + frac * singleStep, P, t + frac * dt,
    //             RefPartBunch_m->Q[i], RefPartBunch_m->M[i]));
    //     }
    // }

    return false;
}

bool Monitor::apply() {
    if (!online_m || type_m != CollectionType::SPATIAL || !lossDs_m || !RefPartBunch_m) {
        return false;
    }

    auto pc = RefPartBunch_m->getParticleContainer();
    const size_t localNum = RefPartBunch_m->getLocalNum();
    if (localNum == 0) {
        deviceData_m.lossCount = 0;
        lossDs_m->bindDeviceData(deviceData_m);
        return false;
    }

    if (lossCapacity_m != localNum) {
        lossCapacity_m = localNum;
        lossIndex_m    = LossDataSink::DeviceData::IndexView_t("MonitorLossIndex", lossCapacity_m);
        lossTime_m     = LossDataSink::DeviceData::TimeView_t("MonitorLossTime", lossCapacity_m);
        lossCount_m    = Kokkos::View<size_t>("MonitorLossCount");
        Kokkos::deep_copy(lossCount_m, size_t(0));

        deviceData_m.lossIndex  = lossIndex_m;
        deviceData_m.Time       = lossTime_m;
        deviceData_m.localCount = localNum;
        deviceData_m.lossCount  = 0;
        lossDs_m->bindDeviceData(deviceData_m);
    }

    auto Rview  = pc->R.getView();
    auto Pview  = pc->P.getView();
    auto dtview = pc->dt.getView();
    auto lossIndex = lossIndex_m;
    auto lossTime  = lossTime_m;
    auto lossCount = lossCount_m;
    const double t0 = RefPartBunch_m->getT();
    const double c  = Physics::c;

    Kokkos::deep_copy(lossCount_m, size_t(0));

    Kokkos::parallel_for(
        "Monitor::applyLossIndex", ippl::getRangePolicy(Rview), KOKKOS_LAMBDA(const int idx) {
            const double Rz = Rview(idx)(2);
            const double dt = dtview(idx);
            if (dt == 0.0) {
                return;
            }

            const double px = Pview(idx)(0);
            const double py = Pview(idx)(1);
            const double pz = Pview(idx)(2);
            const double gamma = Kokkos::sqrt(1.0 + px * px + py * py + pz * pz);
            const double betaZ = pz / gamma;
            const double singleStepZ = c * dt * betaZ;
            if (singleStepZ == 0.0) {
                return;
            }

            if (dt * Rz < 0.0 && dt * (Rz + singleStepZ) > 0.0) {
                const double frac = -Rz / singleStepZ;
                const double time = t0 + frac * dt;
                const size_t slot = Kokkos::atomic_fetch_add(&lossCount(), size_t(1));
                if (slot < lossIndex.extent(0)) {
                    lossIndex(slot) = static_cast<size_t>(idx);
                    lossTime(idx)   = time;
                }
            }
        });
    Kokkos::fence();

    size_t hostCount = 0;
    Kokkos::deep_copy(hostCount, lossCount_m);
    deviceData_m.lossCount = hostCount;
    lossDs_m->bindDeviceData(deviceData_m);

    return false;
}

bool Monitor::apply(
    const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
    Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    return false;
}

void Monitor::driftToCorrectPositionAndSave(
    const Vector_t<double, 3>& refR, const Vector_t<double, 3>& refP) {
    // const double cdt                           = Physics::c * RefPartBunch_m->getdT();
    // const Vector_t<double, 3> driftPerTimeStep = cdt * Util::getBeta(refP);
    // const double tau                           = -refR(2) / driftPerTimeStep(2);
    // const CoordinateSystemTrafo update(
    //     refR + tau * driftPerTimeStep, getQuaternion(refP, Vector_t<double, 3>(0, 0, 1)));
    // const CoordinateSystemTrafo refToLocalCSTrafo =
    //     update * (getCSTrafoGlobal2Local() * RefPartBunch_m->toLabTrafo_m);

    // for (OpalParticle particle : *RefPartBunch_m) {
    //     Vector_t<double, 3> beta = refToLocalCSTrafo.rotateTo(Util::getBeta(particle.getP()));
    //     Vector_t<double, 3> dS =
    //         (tau - 0.5) * cdt
    //         * beta;  // the particles are half a step ahead relative to the reference particle
    //     particle.setR(refToLocalCSTrafo.transformTo(particle.getR()) + dS);
    //     lossDs_m->addParticle(particle);
    // }
}

bool Monitor::applyToReferenceParticle(
    const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
    Vector_t<double, 3>&, Vector_t<double, 3>&) {
    // if (!OpalData::getInstance()->isInPrepState()) {
    //     const double dt                      = RefPartBunch_m->getdT();
    //     const double cdt                     = Physics::c * dt;
    //     const Vector_t<double, 3> singleStep = cdt * Util::getBeta(P);

    //     if (dt * R(2) < 0.0 && dt * (R(2) + singleStep(2)) > 0.0) {
    //         double frac            = -R(2) / singleStep(2);
    //         double time            = t + frac * dt;
    //         Vector_t<double, 3> dR = frac * singleStep;
    //         double ds              = euclidean_norm(dR + 0.5 * singleStep);
    //         lossDs_m->addReferenceParticle(
    //             csTrafoGlobal2Local_m.transformFrom(R + dR), csTrafoGlobal2Local_m.rotateFrom(P),
    //             time, RefPartBunch_m->get_sPos() + ds, RefPartBunch_m->getGlobalTrackStep());

    //         if (type_m == CollectionType::TEMPORAL) {
    //             driftToCorrectPositionAndSave(R, P);
    //             auto stats = lossDs_m->computeStatistics(1);
    //             if (!stats.empty()) {
    //                 statFileEntries_sm.insert(
    //                     std::make_pair(stats.begin()->spos_m, *stats.begin()));
    //                 OpalData::OpenMode openMode;
    //                 if (numPassages_m > 0) {
    //                     openMode = OpalData::OpenMode::APPEND;
    //                 } else {
    //                     openMode = OpalData::getInstance()->getOpenMode();
    //                 }
    //                 lossDs_m->save(1, openMode);
    //             }
    //         }

    //         ++numPassages_m;
    //     }
    // }
    return false;
}

void Monitor::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    RefPartBunch_m = bunch;
    endField       = startField + halfLength_s;
    startField -= halfLength_s;

    const size_t totalNum  = bunch->getTotalNum();
    double currentPosition = endField;
    if (totalNum > 0) {
        currentPosition = bunch->get_sPos();
    }

    filename_m = getOutputFN();

    // if (OpalData::getInstance()->getOpenMode() == OpalData::OpenMode::WRITE
    //     || currentPosition < startField) {
    //     namespace fs = boost::filesystem;

    //     fs::path lossFileName = fs::path(filename_m + ".h5");
    //     if (fs::exists(lossFileName)) {
    //         Ippl::Comm->barrier();
    //         if (ippl::Comm->rank() == 0) {
    //             fs::remove(lossFileName);
    //         }
    //         Ippl::Comm->barrier();
    //     }
    // }

    lossDs_m =
        std::unique_ptr<LossDataSink>(new LossDataSink(filename_m, !Options::asciidump, type_m));

    auto pc = RefPartBunch_m->getParticleContainer();
    deviceData_m.R = pc->R.getView();
    deviceData_m.P = pc->P.getView();
    deviceData_m.Q = pc->Q.getView();
    deviceData_m.M = pc->M.getView();
    deviceData_m.localCount = RefPartBunch_m->getLocalNum();
    deviceData_m.lossCount  = 0;

    lossCapacity_m = deviceData_m.localCount;
    lossIndex_m    = LossDataSink::DeviceData::IndexView_t("MonitorLossIndex", lossCapacity_m);
    lossTime_m     = LossDataSink::DeviceData::TimeView_t("MonitorLossTime", lossCapacity_m);
    lossCount_m    = Kokkos::View<size_t>("MonitorLossCount");
    Kokkos::deep_copy(lossCount_m, size_t(0));

    deviceData_m.lossIndex = lossIndex_m;
    deviceData_m.Time      = lossTime_m;
    lossDs_m->bindDeviceData(deviceData_m);
}

void Monitor::finalise() {
}

void Monitor::goOnline(const double&) {
    online_m = true;
}

void Monitor::goOffline() {
    auto stats = lossDs_m->computeStatistics(numPassages_m);
    for (auto& stat : stats) {
        statFileEntries_sm.insert(std::make_pair(stat.spos_m, stat));
    }

    if (type_m != CollectionType::TEMPORAL) {
        lossDs_m->save(numPassages_m);
    }
}

bool Monitor::bends() const {
    return false;
}

void Monitor::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = -halfLength_s;
    zEnd   = halfLength_s;
}

ElementType Monitor::getType() const {
    // return ElementType::MONITOR;
    throw std::runtime_error("Fix this function please");
}

void Monitor::writeStatistics() {
    if (statFileEntries_sm.size() == 0)
        return;

    std::string fileName =
        OpalData::getInstance()->getInputBasename() + std::string("_Monitors.stat");
    auto instance      = OpalData::getInstance();
    bool hasPriorTrack = instance->hasPriorTrack();
    bool inRestartRun  = instance->inRestartRun();

    auto it     = statFileEntries_sm.begin();
    double spos = it->first;
    Util::rewindLinesSDDS(fileName, spos, false);

    MonitorStatisticsWriter writer(fileName, hasPriorTrack || inRestartRun);

    for (const auto& entry : statFileEntries_sm) {
        writer.addRow(entry.second);
    }

    statFileEntries_sm.clear();
}
