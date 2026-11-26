//
// Class Monitor
//   Defines a beam position monitor.
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

#include "Fields/Fieldmap.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/LossDataSink.h"
#include "Track/RunDescription.h"
#include "Utilities/Util.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>

const double Monitor::halfLength_s = 0.005;

Monitor::Monitor() : Monitor("") {}

Monitor::Monitor(const Monitor& right)
    : ElementBase(right),
      filename_m(right.filename_m),
      type_m(right.type_m),
      numPassages_m(0),
      asciiDump_m(right.asciiDump_m),
      computePercentiles_m(right.computePercentiles_m),
      restart_m(right.restart_m),
      appendOutput_m(right.appendOutput_m) {}

Monitor::Monitor(const std::string& name)
    : ElementBase(name),
      filename_m(""),
      type_m(CollectionType::SPATIAL),
      numPassages_m(0) {}

Monitor::~Monitor() {}

void Monitor::apply(const std::shared_ptr<ParticleContainer_t>& pc, const StepTime& time) {
    if (!online_m || lossDs_m == nullptr || pc == nullptr) {
        return;
    }

    if (type_m != CollectionType::SPATIAL) {
        return;
    }

    const long long globalStep = time.step;
    const double sPos          = pc->get_sPos();

    if (globalStep == 0 && std::abs(sPos) < 1.0e-14) {
        Inform msg("Monitor::apply(pc)");
        msg << level5 << "Ignoring pre-tracking spatial particle crossing"
            << " globalStep=" << globalStep << " pc_spos=" << sPos << endl;
        return;
    }

    const auto nLoc = pc->getLocalNum();
    if (nLoc == 0) {
        return;
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

    const bool hasSpin = pc->hasSpin();
    using SpinHostView = Kokkos::View<ippl::Vector<double, 3>*, Kokkos::HostSpace>;
    SpinHostView hPol;
    if (hasSpin) {
        auto Polview = pc->Pol.getView();
        hPol         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Polview);
    }

    const bool qmAreAttributes =
            (pc->getQMStorageMode() == ParticleContainer_t::QMStorageMode::Attributes);

    const double bunchTime = time.t;

    for (size_t i = 0; i < nLoc; ++i) {
        const Vector_t<double, 3> R = hR(i);
        const Vector_t<double, 3> P = hP(i);
        const double dt             = hdt(i);

        const Vector_t<double, 3> singleStep = Physics::c * dt * Util::getBeta(P);

        if (singleStep(2) == 0.0) {
            continue;
        }

        if (dt * R(2) < 0.0 && dt * (R(2) + singleStep(2)) > 0.0) {
            const double frac = -R(2) / singleStep(2);

            const Vector_t<double, 3> crossingR = R + frac * singleStep;
            const double crossingTime           = bunchTime + 0.5 * time.dt + frac * dt;

            const std::size_t id = static_cast<std::size_t>(hID(i));
            const double q       = qmAreAttributes ? hQ(i) : hQ(0);
            const double m       = qmAreAttributes ? hM(i) : hM(0);

            if (hasSpin) {
                const Vector_t<double, 3> pol(
                        static_cast<double>(hPol(i)[0]), static_cast<double>(hPol(i)[1]),
                        static_cast<double>(hPol(i)[2]));
                lossDs_m->addParticle(OpalParticle(id, crossingR, P, crossingTime, q, m, pol));
            } else {
                lossDs_m->addParticle(OpalParticle(id, crossingR, P, crossingTime, q, m));
            }
        }
    }
}

void Monitor::driftToCorrectPositionAndSave(
        const Vector_t<double, 3>& refR, const Vector_t<double, 3>& refP) {
    if (lossDs_m == nullptr || RefPartBunch_m == nullptr) {
        return;
    }

    const std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    if (!pc) {
        return;
    }
    pc->requireReferenceFrame("Monitor::driftToCorrectPositionAndSave");

    const double dt  = RefPartBunch_m->getdT();
    const double cdt = Physics::c * dt;

    const Vector_t<double, 3> driftPerTimeStep = cdt * Util::getBeta(refP);

    if (driftPerTimeStep(2) == 0.0) {
        return;
    }

    const double tau = -refR(2) / driftPerTimeStep(2);

    const CoordinateSystemTrafo update(
            refR + tau * driftPerTimeStep, getQuaternion(refP, Vector_t<double, 3>(0.0, 0.0, 1.0)));

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

    const bool hasSpin = pc->hasSpin();
    using SpinHostView = Kokkos::View<ippl::Vector<double, 3>*, Kokkos::HostSpace>;
    SpinHostView hPol;
    if (hasSpin) {
        auto Polview = pc->Pol.getView();
        hPol         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Polview);
    }

    const auto nLoc = pc->getLocalNum();

    const bool qmAreAttributes =
            (pc->getQMStorageMode() == ParticleContainer_t::QMStorageMode::Attributes);

    const double particleTime = RefPartBunch_m->getT();

    for (size_t i = 0; i < nLoc; ++i) {
        const Vector_t<double, 3> globalR = hR(i);
        const Vector_t<double, 3> globalP = hP(i);

        const Vector_t<double, 3> beta = refToLocalCSTrafo.rotateTo(Util::getBeta(globalP));

        const Vector_t<double, 3> dS = (tau - 0.5) * cdt * beta;

        const Vector_t<double, 3> localR = refToLocalCSTrafo.transformTo(globalR) + dS;

        const std::size_t id = static_cast<std::size_t>(hID(i));

        const double q = qmAreAttributes ? hQ(i) : hQ(0);

        const double macroMassGeV = qmAreAttributes ? hM(i) : hM(0);
        const double macroWeight  = std::abs(q) / std::abs(Physics::q_e);

        const double m = macroWeight > 0.0 ? macroMassGeV * Units::GeV2MeV / macroWeight
                                           : macroMassGeV * Units::GeV2MeV;

        if (hasSpin) {
            const Vector_t<double, 3> pol(
                    static_cast<double>(hPol(i)[0]), static_cast<double>(hPol(i)[1]),
                    static_cast<double>(hPol(i)[2]));
            lossDs_m->addParticle(OpalParticle(id, localR, globalP, particleTime, q, m, pol));
        } else {
            lossDs_m->addParticle(OpalParticle(id, localR, globalP, particleTime, q, m));
        }
    }
}

void Monitor::recordReference(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double t) {
    if (!online_m || lossDs_m == nullptr || RefPartBunch_m == nullptr || !recordReference_m) {
        return;
    }

    const std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    if (!pc) {
        return;
    }

    const double dt                      = RefPartBunch_m->getdT();
    const double cdt                     = Physics::c * dt;
    const Vector_t<double, 3> singleStep = cdt * Util::getBeta(P);

    if (singleStep(2) == 0.0) {
        return;
    }

    if (dt * R(2) < 0.0 && dt * (R(2) + singleStep(2)) > 0.0) {
        const long long globalStep = RefPartBunch_m->getGlobalTrackStep();
        const double sPos          = pc->get_sPos();

        if (globalStep == 0 && std::abs(sPos) < 1.0e-14) {
            Inform msg("Monitor::recordReference");
            msg << level5 << "Ignoring pre-tracking reference crossing"
                << " globalStep=" << globalStep << " pc_spos=" << sPos << endl;
            return;
        }

        const double frac            = -R(2) / singleStep(2);
        const double time            = t + frac * dt;
        const Vector_t<double, 3> dR = frac * singleStep;

        const Vector_t<double, 3> dsVec = dR + 0.5 * singleStep;
        const double ds                 = euclidean_norm(dsVec);

        lossDs_m->addReferenceParticle(
                csTrafoGlobal2Local_m.transformFrom(R + dR), csTrafoGlobal2Local_m.rotateFrom(P),
                time, pc->get_sPos() + ds, RefPartBunch_m->getGlobalTrackStep());

        if (type_m == CollectionType::TEMPORAL) {
            driftToCorrectPositionAndSave(R, P);

            auto stats = lossDs_m->computeStatistics(1);
            if (!stats.empty()) {
                lossDs_m->save(1, numPassages_m > 0 || appendOutput_m);
            }
        }

        ++numPassages_m;
    }
}

void Monitor::setRunValues(
        const opalx::run::OptionValues& options, const opalx::run::RunContext& context) {
    asciiDump_m          = options.asciidump;
    computePercentiles_m = options.computePercentiles;
    restart_m            = context.restart;
    appendOutput_m       = context.openMode == "APPEND";
}

void Monitor::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;

    // The monitor's path-length position (its ELEMEDGE, or 0 for a 6D-posed monitor); the
    // monitor plane straddles it by +-halfLength_s.
    const double elemEdge     = isElementPositionSet() ? getElementPosition() : 0.0;
    const double monitorStart = elemEdge - halfLength_s;

    double currentPosition = elemEdge + halfLength_s;
    if (bunch != nullptr) {
        const std::shared_ptr<ParticleContainer_t> pc = bunch->getParticleContainer();
        if (pc) {
            currentPosition = pc->get_sPos();
        }
    }

    filename_m = getOutputFN();

    namespace fs = std::filesystem;
    const fs::path lossFileName(filename_m + ".h5");
    if (restart_m) {
        if (!asciiDump_m && bunch != nullptr) {
            LossDataSink::rewindH5ToGlobalTrackStep(
                    lossFileName.string(), bunch->getGlobalTrackStep());
        }
    } else if (!appendOutput_m || currentPosition < monitorStart) {
        if (fs::exists(lossFileName)) {
            ippl::Comm->barrier();
            if (ippl::Comm->rank() == 0) {
                fs::remove(lossFileName);
            }
            ippl::Comm->barrier();
        }
    }

    lossDs_m =
            std::make_unique<LossDataSink>(filename_m, !asciiDump_m, type_m, computePercentiles_m);
    online_m = true;
}

void Monitor::goOffline() {
    if (!lossDs_m) {
        return;
    }

    if (type_m != CollectionType::TEMPORAL) {
        lossDs_m->save(numPassages_m, appendOutput_m);
    }
}

void Monitor::getFieldExtent(double& zBegin, double& zEnd) const {
    zBegin = -halfLength_s;
    zEnd   = halfLength_s;
}

BoundingBox Monitor::getLocalSupportBox(const double reach) const {
    // The base box covers the body [0, L), where an APERTURE deletes particles; add the
    // isInside() window |z| <= L/2 around the monitor plane.
    BoundingBox box           = ElementBase::getLocalSupportBox(reach);
    const auto [lower, upper] = box.getCorners();
    const double halfLength   = 0.5 * getGeometry().getElementLength();
    box.enlargeToContainPosition(Vector_t<double, 3>({lower(0), lower(1), -halfLength}));
    box.enlargeToContainPosition(Vector_t<double, 3>({upper(0), upper(1), halfLength}));
    return box;
}

/// \todo Check switch statements over ElementType and add a MONITOR case where needed.
ElementType Monitor::getType() const { return ElementType::MONITOR; }
