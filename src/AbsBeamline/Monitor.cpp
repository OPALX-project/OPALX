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
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Structure/MonitorStatisticsWriter.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <Kokkos_Core.hpp>
#include <filesystem>
#include <limits>
#include <memory>

std::map<double, SetStatistics> Monitor::statFileEntries_sm;
const double Monitor::halfLength_s = 0.005;

Monitor::Monitor() : Monitor("") {}

Monitor::Monitor(const Monitor& right)
    : Component(right),
      filename_m(right.filename_m),
      plane_m(right.plane_m),
      type_m(right.type_m),
      numPassages_m(0) {}

Monitor::Monitor(const std::string& name)
    : Component(name),
      filename_m(""),
      plane_m(OFF),
      type_m(CollectionType::SPATIAL),
      numPassages_m(0) {}

Monitor::~Monitor() {}

void Monitor::accept(BeamlineVisitor& visitor) const { visitor.visitMonitor(*this); }

bool Monitor::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    if (!online_m || type_m != CollectionType::SPATIAL || !pc) {
        return false;
    }

    const size_t localNum = pc->getLocalNum();
    const double time     = RefPartBunch_m->getT() - 0.5 * RefPartBunch_m->getdT();
    Vector_t<double, 3> E(0.0);
    Vector_t<double, 3> B(0.0);

    for (size_t i = 0; i < localNum; ++i) {
        apply(i, time, E, B);
    }

    return false;
}

bool Monitor::apply(
        const size_t& i, const double& t, Vector_t<double, 3>& /*E*/,
        Vector_t<double, 3>& /*B*/) {
    const auto pc = RefPartBunch_m->getParticleContainer();
    if (!pc) {
        return false;
    }

    auto Rhost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
    auto Phost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
    auto dthost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->dt.getView());
    auto idHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());
    auto qHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getQView());
    auto mHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getMView());

    const Vector_t<double, 3> R = Rhost(i);
    const Vector_t<double, 3> P = Phost(i);
    const double dt             = dthost(i);
    const Vector_t<double, 3> ds = Physics::c * dt * Util::getBeta(P);
    const double q = (qHost.extent(0) == 1) ? qHost(0) : qHost(i);
    const double m = (mHost.extent(0) == 1) ? mHost(0) : mHost(i);
    const int64_t id = (idHost.extent(0) > i) ? idHost(i) : static_cast<int64_t>(i);

    if (online_m && type_m == CollectionType::SPATIAL
            && std::abs(ds(2)) > std::numeric_limits<double>::epsilon()) {
        const double z0 = R(2);
        const double z1 = R(2) + ds(2);
        if (z0 * z1 <= 0.0) {
            const double frac = -z0 / ds(2);
            if (frac >= 0.0 && frac <= 1.0) {
                lossDs_m->addParticle(OpalParticle(id, R + frac * ds, P, t + frac * dt, q, m));
            }
        }
    }

    return false;
}

bool Monitor::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {
    if (!online_m || type_m != CollectionType::SPATIAL) {
        return false;
    }

    const double dt              = RefPartBunch_m->getdT();
    const Vector_t<double, 3> ds = Physics::c * dt * Util::getBeta(P);
    if (std::abs(ds(2)) > std::numeric_limits<double>::epsilon()) {
        const double z0 = R(2);
        const double z1 = R(2) + ds(2);
        if (z0 * z1 <= 0.0) {
            const double frac = -z0 / ds(2);
            if (frac >= 0.0 && frac <= 1.0) {
                lossDs_m->addParticle(OpalParticle(-1, R + frac * ds, P, t + frac * dt, 0.0, 0.0));
            }
        }
    }
    return false;
}

void Monitor::driftToCorrectPositionAndSave(
        const Vector_t<double, 3>& refR, const Vector_t<double, 3>& refP) {
    const double cdt                       = Physics::c * RefPartBunch_m->getdT();
    const Vector_t<double, 3> driftPerStep = cdt * Util::getBeta(refP);
    if (std::abs(driftPerStep(2)) <= std::numeric_limits<double>::epsilon()) {
        return;
    }

    const double tau = -refR(2) / driftPerStep(2);
    const CoordinateSystemTrafo update(
            refR + tau * driftPerStep, getQuaternion(refP, Vector_t<double, 3>(0, 0, 1)));
    const CoordinateSystemTrafo refToLocal =
            update
            * (getCSTrafoGlobal2Local() * RefPartBunch_m->getParticleContainer()->getToLabTrafo());

    const auto pc = RefPartBunch_m->getParticleContainer();
    auto Rhost    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->R.getView());
    auto Phost    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->P.getView());
    auto idHost   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->ID.getView());
    auto qHost    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getQView());
    auto mHost    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->getMView());

    const size_t nLocal = pc->getLocalNum();
    for (size_t i = 0; i < nLocal; ++i) {
        const int64_t id = (idHost.extent(0) > i) ? idHost(i) : static_cast<int64_t>(i);
        const double q   = (qHost.extent(0) == 1) ? qHost(0) : qHost(i);
        const double m   = (mHost.extent(0) == 1) ? mHost(0) : mHost(i);
        OpalParticle particle(id, Rhost(i), Phost(i), 0.0, q, m);
        Vector_t<double, 3> beta = refToLocal.rotateTo(Util::getBeta(particle.getP()));
        Vector_t<double, 3> dS   = (tau - 0.5) * cdt * beta;
        particle.setR(refToLocal.transformTo(particle.getR()) + dS);
        lossDs_m->addParticle(particle);
    }
}

bool Monitor::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {
    if (!OpalData::getInstance()->isInPrepState()) {
        const double dt                       = RefPartBunch_m->getdT();
        const double cdt                      = Physics::c * dt;
        const Vector_t<double, 3> singleStep = cdt * Util::getBeta(P);

        if (std::abs(singleStep(2)) > std::numeric_limits<double>::epsilon()) {
            const double frac             = -R(2) / singleStep(2);
            const double z0               = R(2);
            const double z1               = R(2) + singleStep(2);
            if (z0 * z1 <= 0.0 && frac >= 0.0 && frac <= 1.0) {
                const double time             = t + frac * dt;
                const Vector_t<double, 3> dR  = frac * singleStep;
                const Vector_t<double, 3> dsR = dR + 0.5 * singleStep;
                const double ds               = euclidean_norm(dsR);
                lossDs_m->addReferenceParticle(
                        csTrafoGlobal2Local_m.transformFrom(R + dR), csTrafoGlobal2Local_m.rotateFrom(P),
                        time, RefPartBunch_m->getParticleContainer()->get_sPos() + ds,
                        RefPartBunch_m->getGlobalTrackStep());

                if (type_m == CollectionType::TEMPORAL) {
                    driftToCorrectPositionAndSave(R, P);
                    auto stats = lossDs_m->computeStatistics(1);
                    if (!stats.empty()) {
                        statFileEntries_sm.insert(std::make_pair(stats.begin()->spos_m, *stats.begin()));
                        OpalData::OpenMode openMode = numPassages_m > 0
                                                             ? OpalData::OpenMode::APPEND
                                                             : OpalData::getInstance()->getOpenMode();
                        lossDs_m->save(1, openMode);
                    }
                }

                ++numPassages_m;
            }
        }
    }
    return false;
}

void Monitor::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    RefPartBunch_m = bunch;
    endField       = startField + halfLength_s;
    startField -= halfLength_s;

    filename_m = getOutputFN();
    const auto pc         = bunch->getParticleContainer();
    const size_t totalNum = pc ? pc->getTotalNum() : 0;
    double currentPos     = endField;
    if (totalNum > 0) {
        currentPos = pc->get_sPos();
    }

    if (OpalData::getInstance()->getOpenMode() == OpalData::OpenMode::WRITE
            || currentPos < startField) {
        namespace fs = std::filesystem;
        fs::path lossFileName = fs::path(filename_m + ".h5");
        if (fs::exists(lossFileName)) {
            ippl::Comm->barrier();
            if (ippl::Comm->rank() == 0) {
                fs::remove(lossFileName);
            }
            ippl::Comm->barrier();
        }
    }

    lossDs_m = std::make_unique<LossDataSink>(filename_m, !Options::asciidump, type_m);
}

void Monitor::finalise() {}

void Monitor::goOnline(const double&) { online_m = true; }

void Monitor::goOffline() {
    auto stats = lossDs_m->computeStatistics(numPassages_m);
    for (auto& stat : stats) {
        statFileEntries_sm.insert(std::make_pair(stat.spos_m, stat));
    }

    if (type_m != CollectionType::TEMPORAL) {
        lossDs_m->save(numPassages_m);
    }
}

bool Monitor::bends() const { return false; }

void Monitor::getFieldExtend(double& zBegin, double& zEnd) const {
    zBegin = -halfLength_s;
    zEnd   = halfLength_s;
}

ElementType Monitor::getType() const {
    return ElementType::MONITOR;
}

void Monitor::writeStatistics() {
    if (statFileEntries_sm.size() == 0) return;

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
