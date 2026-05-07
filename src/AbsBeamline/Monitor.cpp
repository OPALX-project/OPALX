//
// Class Monitor
//   Defines the abstract interface for a beam position monitor.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#include <cmath>
#include <filesystem>
#include <fstream>
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
    if (!online_m || lossDs_m == nullptr || pc == nullptr || RefPartBunch_m == nullptr) {
        return false;
    }

    if (type_m != CollectionType::SPATIAL) {
        return false;
    }

    const long long globalStep = RefPartBunch_m->getGlobalTrackStep();
    const double sPos = pc->get_sPos();

    if (globalStep == 0 && std::abs(sPos) < 1.0e-14) {
        Inform msg("Monitor::apply(pc)");
        msg << level5
            << "Ignoring pre-tracking spatial particle crossing"
            << " globalStep=" << globalStep
            << " pc_spos=" << sPos
            << endl;
        return false;
    }

    const auto nLoc = pc->getLocalNum();
    if (nLoc == 0) {
        return false;
    }

    auto Rview  = pc->R.getView();
    auto Pview  = pc->P.getView();
    auto dtview = pc->dt.getView();
    auto IDview = pc->ID.getView();
    auto Qview  = pc->getQView();
    auto Mview  = pc->getMView();

    auto hR  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Rview);
    auto hP  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Pview);
    auto hdt = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dtview);
    auto hID = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), IDview);
    auto hQ  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Qview);
    auto hM  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Mview);

    const bool qmAreAttributes =
        (pc->getQMStorageMode() == ParticleContainer_t::QMStorageMode::Attributes);

    const double bunchTime = RefPartBunch_m->getT();

    for (size_t i = 0; i < nLoc; ++i) {
        const Vector_t<double, 3> R = hR(i);
        const Vector_t<double, 3> P = hP(i);
        const double dt = hdt(i);

        const Vector_t<double, 3> singleStep = Physics::c * dt * Util::getBeta(P);

        if (singleStep(2) == 0.0) {
            continue;
        }

        if (dt * R(2) < 0.0 && dt * (R(2) + singleStep(2)) > 0.0) {
            const double frac = -R(2) / singleStep(2);

            const Vector_t<double, 3> crossingR = R + frac * singleStep;
            const double crossingTime = 
                        bunchTime + 0.5 * RefPartBunch_m->getdT() + frac * dt;

            const std::size_t id = static_cast<std::size_t>(hID(i));
            const double q = qmAreAttributes ? hQ(i) : hQ(0);
            const double m = qmAreAttributes ? hM(i) : hM(0);

            lossDs_m->addParticle(OpalParticle(id, crossingR, P, crossingTime, q, m));
        }
    }

    return false;
}

bool Monitor::apply(
        const size_t& i, const double& t, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {

    if (!online_m || lossDs_m == nullptr || RefPartBunch_m == nullptr) {
        return false;
    }

    if (type_m != CollectionType::SPATIAL) {
        return false;
    }

    const std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    if (!pc) {
        return false;
    }

    auto Rview  = pc->R.getView();
    auto Pview  = pc->P.getView();
    auto dtview = pc->dt.getView();
    auto IDview = pc->ID.getView();
    auto Qview  = pc->getQView();
    auto Mview  = pc->getMView();

    const Vector_t<double, 3> R = Rview(i);
    const Vector_t<double, 3> P = Pview(i);
    const double dt = dtview(i);

    const Vector_t<double, 3> singleStep = Physics::c * dt * Util::getBeta(P);

    if (singleStep(2) == 0.0) {
        return false;
    }

    if (dt * R(2) < 0.0 && dt * (R(2) + singleStep(2)) > 0.0) {
        const double frac = -R(2) / singleStep(2);
        const Vector_t<double, 3> crossingR = R + frac * singleStep;
        const double crossingTime = t + frac * dt;

        const bool qmAreAttributes =
            (pc->getQMStorageMode() == ParticleContainer_t::QMStorageMode::Attributes);

        const std::size_t id = static_cast<std::size_t>(IDview(i));
        const double q = qmAreAttributes ? Qview(i) : Qview(0);
        const double m = qmAreAttributes ? Mview(i) : Mview(0);

        lossDs_m->addParticle(
            OpalParticle(id, crossingR, P, crossingTime, q, m));
    }

    return false;
}

bool Monitor::apply(
        const Vector_t<double, 3>& /*R*/, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {

    // Monitor is field-free.
    // This overload does not provide particle ID and may not provide the correct
    // per-particle dt/q/m information, so it cannot safely save monitor hits.
    return false;
}

void Monitor::driftToCorrectPositionAndSave(
        const Vector_t<double, 3>& refR, 
        const Vector_t<double, 3>& refP) {

    if (lossDs_m == nullptr || RefPartBunch_m == nullptr) {
        return;
    }

    const std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    if (!pc) {
        return;
    }

    const double dt  = RefPartBunch_m->getdT();
    const double cdt = Physics::c * dt;

    const Vector_t<double, 3> driftPerTimeStep = cdt * Util::getBeta(refP);

    if (driftPerTimeStep(2) == 0.0) {
        return;
    }

    const double tau = -refR(2) / driftPerTimeStep(2);

    const CoordinateSystemTrafo update(
        refR + tau * driftPerTimeStep,
        getQuaternion(refP, Vector_t<double, 3>(0.0, 0.0, 1.0)));

    const CoordinateSystemTrafo refToLocalCSTrafo =
        update * (getCSTrafoGlobal2Local() * pc->getToLabTrafo());

    auto Rview  = pc->R.getView();
    auto Pview  = pc->P.getView();
    auto IDview = pc->ID.getView();
    auto Qview  = pc->getQView();
    auto Mview  = pc->getMView();

    auto hR  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Rview);
    auto hP  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Pview);
    auto hID = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), IDview);
    auto hQ  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Qview);
    auto hM  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Mview);

    const auto nLoc = pc->getLocalNum();

    const bool qmAreAttributes =
        (pc->getQMStorageMode() == ParticleContainer_t::QMStorageMode::Attributes);

    // const double crossingTime =
    //     RefPartBunch_m->getT() + tau * RefPartBunch_m->getdT();
    const double particleTime = RefPartBunch_m->getT();

    for (size_t i = 0; i < nLoc; ++i) {
        const Vector_t<double, 3> globalR = hR(i);
        const Vector_t<double, 3> globalP = hP(i);

        const Vector_t<double, 3> beta =
            refToLocalCSTrafo.rotateTo(Util::getBeta(globalP));

        const Vector_t<double, 3> dS =
            (tau - 0.5) * cdt * beta;

        const Vector_t<double, 3> localR =
            refToLocalCSTrafo.transformTo(globalR) + dS;
        // const Vector_t<double, 3> localP =
        //     refToLocalCSTrafo.rotateTo(globalP);

        const std::size_t id = static_cast<std::size_t>(hID(i));

        const double q = qmAreAttributes ? hQ(i) : hQ(0);
        // const double m = qmAreAttributes ? hM(i) : hM(0);

        const double macroMassGeV = qmAreAttributes ? hM(i) : hM(0);
        const double macroWeight =
            std::abs(q) / std::abs(Physics::q_e);

        const double m =
            macroWeight > 0.0
                ? macroMassGeV * Units::GeV2MeV / macroWeight
                : macroMassGeV * Units::GeV2MeV;

        // lossDs_m->addParticle(
        //     OpalParticle(id, localR, localP, crossingTime, q, m));
        lossDs_m->addParticle(
            OpalParticle(id, localR, globalP, particleTime, q, m));
    }
}

bool Monitor::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {

    if (!online_m || lossDs_m == nullptr || RefPartBunch_m == nullptr
        || OpalData::getInstance()->isInPrepState()) {
        return false;
    }

    const std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    if (!pc) {
        return false;
    }

    const double dt = RefPartBunch_m->getdT();
    const double cdt = Physics::c * dt;
    const Vector_t<double, 3> singleStep = cdt * Util::getBeta(P);

    if (singleStep(2) == 0.0) {
        return false;
    }

    if (dt * R(2) < 0.0 && dt * (R(2) + singleStep(2)) > 0.0) {
        const long long globalStep = RefPartBunch_m->getGlobalTrackStep();
        const double sPos = pc->get_sPos();

        if (globalStep == 0 && std::abs(sPos) < 1.0e-14) {
            Inform msg("Monitor::applyToReferenceParticle");
            msg << level5
                << "Ignoring pre-tracking reference crossing"
                << " globalStep=" << globalStep
                << " pc_spos=" << sPos
                << endl;
            return false;
        }

        const double frac = -R(2) / singleStep(2);
        const double time = t + frac * dt;
        const Vector_t<double, 3> dR = frac * singleStep;

        const Vector_t<double, 3> dsVec = dR + 0.5 * singleStep;
        const double ds = euclidean_norm(dsVec);

        lossDs_m->addReferenceParticle(
            csTrafoGlobal2Local_m.transformFrom(R + dR),
            csTrafoGlobal2Local_m.rotateFrom(P),
            time,
            pc->get_sPos() + ds,
            RefPartBunch_m->getGlobalTrackStep());

        if (type_m == CollectionType::TEMPORAL) {
            driftToCorrectPositionAndSave(R, P);

            auto stats = lossDs_m->computeStatistics(1);
            if (!stats.empty()) {
                statFileEntries_sm.insert(std::make_pair(stats.begin()->spos_m, *stats.begin()));

                OpalData::OpenMode openMode =
                    (numPassages_m > 0)
                        ? OpalData::OpenMode::APPEND
                        : OpalData::getInstance()->getOpenMode();

                lossDs_m->save(1, openMode);
            }
        }

        ++numPassages_m;
    }

    return false;
}

void Monitor::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    RefPartBunch_m = bunch;
    
    endField   = startField + halfLength_s;
    startField -= halfLength_s;

    double currentPosition = endField;
    if (bunch != nullptr) {
        const std::shared_ptr<ParticleContainer_t> pc = bunch->getParticleContainer();
        if (pc) {
            currentPosition = pc->get_sPos();
        }
    }

    filename_m = getOutputFN();

    if (OpalData::getInstance()->getOpenMode() == OpalData::OpenMode::WRITE
        || currentPosition < startField) {
        namespace fs = std::filesystem;

        const fs::path lossFileName(filename_m + ".h5");
        if (fs::exists(lossFileName)) {
            ippl::Comm->barrier();
            if (ippl::Comm->rank() == 0) {
                fs::remove(lossFileName);
            }
            ippl::Comm->barrier();
        }
    }

    lossDs_m = std::make_unique<LossDataSink>(
        filename_m,
        !Options::asciidump,
        type_m);
}

void Monitor::finalise() {}

void Monitor::goOnline(const double&) { online_m = true; }

void Monitor::goOffline() {
    if (!lossDs_m) {
        return;
    }

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

// One small thing to keep in mind: if there are any switch statements 
// over ElementType elsewhere in OPALX, they may need a MONITOR case 
// later, even if the compiler does not complain yet.
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
