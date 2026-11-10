// Copyright (c) 2023, Chris Rogers
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

#include "AbsBeamline/OutputPlane.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "LossDataSink.h"
#include "PartBunch.h"

std::list<OutputPlane*> OutputPlane::outputPlanes_m;

OutputPlane::OutputPlane() : OutputPlane("") { geometry_m = Geometry::makeNull(); }

OutputPlane::OutputPlane(const std::string& name) : PluginElement(name) {}

OutputPlane::OutputPlane(const OutputPlane& right) : PluginElement(right) {
    normal_m       = right.normal_m;
    centre_m       = right.centre_m;
    radialExtent_m = right.radialExtent_m;
}

OutputPlane::~OutputPlane() { outputPlanes_m.remove(this); }

ElementBase* OutputPlane::clone() const {
    auto* result = new OutputPlane(*this);
    outputPlanes_m.push_back(result);
    return result;
}

void OutputPlane::accept(BeamlineVisitor& visitor) const { visitor.visitOutputPlane(*this); }

ElementType OutputPlane::getType() const { return ElementType::OUTPUTPLANE; }

void OutputPlane::doInitialise(PartBunch_t* /*bunch*/) {}

void OutputPlane::doGoOffline() {}

void OutputPlane::doFinalise() {}

bool OutputPlane::doCheck(
        PartBunch_t* /*bunch*/, int /*turnnumber*/, double /*t*/, double /*tstep*/) {
    return true;
}

void OutputPlane::getFieldExtent(double& zBegin, double& zEnd) const {
    zBegin = -0.005;
    zEnd   = 0.005;
}

void OutputPlane::observe(PartBunch_t* bunch) {
    for (auto* plane : outputPlanes_m) {
        if (!plane->online_m) {
            plane->initialise(bunch);
        }
        if (!plane->stash_m.empty()) {
            plane->detectCrossings(bunch, bunch->getT(), bunch->getdT());
            plane->writeParticles(1);
        }
        plane->stashParticles(bunch);
    }
}

void OutputPlane::store() {
    for (auto* plane : outputPlanes_m) {
        plane->goOffline();
    }
}

void OutputPlane::detectCrossings(
        const PartBunch_t* bunch, const double t, const double tStep) const {
    // Detect the crossings of particles with the output plane
    for (size_t i = 0; i < bunch->getParticleContainers().size(); ++i) {
        auto* pc = bunch->getParticleContainers()[i].get();
        if (i < stash_m.size()) {
            auto* stash             = stash_m[i].get();
            const auto r1           = stash->r();
            const auto r2           = pc->R.getView();
            const auto p            = pc->P.getView();
            const auto invalid      = pc->InvalidMask.getView();
            const auto crossing     = stash->r();
            const auto crossed      = stash->crossed();
            const auto crossingT    = stash->t();
            const auto crossingP    = stash->p();
            const auto centre       = centre_m;
            const auto normal       = normal_m;
            const auto radialExtent = radialExtent_m;
            const auto refR         = pc->getRefPartR();
            Kokkos::parallel_for(
                    "OutputPlane::detectCrossings", r1.extent(0), KOKKOS_LAMBDA(const size_t j) {
                        doDetectCrossing(
                                j, r1, r2, p, crossing, crossingP, crossed, invalid, crossingT, t,
                                tStep, centre, normal, radialExtent, refR);
                    });
        }
    }
    Kokkos::fence();
}

void OutputPlane::writeParticles(const int turnNumber) const {
    // Write out particles that crossed the output plane
    for (size_t i = 0; i < stash_m.size(); ++i) {
        auto* stash          = stash_m[i].get();
        const auto crossing  = Kokkos::create_mirror_view_and_copy(stash->r());
        const auto crossingP = Kokkos::create_mirror_view_and_copy(stash->p());
        const auto crossed   = Kokkos::create_mirror_view_and_copy(stash->crossed());
        const auto crossingT = Kokkos::create_mirror_view_and_copy(stash->t());
        for (size_t j = 0; j < crossed.extent(0); ++j) {
            if (crossed(j) && lossDs_m) {
                lossDs_m->addParticle(
                        OpalParticle(j, crossing(j), crossingP(j), crossingT(j), 0, 0),
                        std::make_pair(turnNumber, i));
            }
        }
    }
}

KOKKOS_FUNCTION void OutputPlane::doDetectCrossing(
        size_t i, VectorView_t r1, VectorView_t r2, VectorView_t p, VectorView_t crossing,
        VectorView_t crossingP, BooleanView_t crossed, BooleanView_t invalid,
        DoubleView_t crossingT, double t, double tStep, Vector3D_t centre, Vector3D_t normal,
        double radialExtent, Vector3D_t refR) {
    if (!invalid(i)) {
        Vector3D_t labR1 = r1(i);
        Vector3D_t labR2 = r2(i) + refR;
        auto d1          = normal.dot(labR1 - centre);
        auto d2          = normal.dot(labR2 - centre);
        crossed(i)       = d1 * d2 <= 0.0 && d1 != d2;
        if (crossed(i)) {
            double dt       = d1 / (d1 - d2);
            crossing(i)     = labR1 + dt * (labR2 - labR1);
            Vector3D_t diff = crossing(i) - centre;
            crossed(i)      = diff.Pnorm() <= radialExtent;
            if (crossed(i)) {
                crossingT(i) = t + dt * tStep;
                crossingP(i) = p(i);
            }
        }
    }
}

void OutputPlane::stashParticles(const PartBunch_t* bunch) {
    // Copy particle positions from the containers to the local store
    for (size_t i = 0; i < bunch->getParticleContainers().size(); ++i) {
        auto* pc = bunch->getParticleContainers()[i].get();
        if (i >= stash_m.size()) {
            stash_m.emplace_back(std::make_unique<StashedContainer>(pc));
        }
        auto* stash = stash_m[i].get();
        stash->stash(pc);
    }
}

OutputPlane::StashedContainer::StashedContainer(ParticleContainer_t* pc)
    : r_m("r", pc->R.getView().extent(0)),
      p_m("p", pc->P.getView().extent(0)),
      crossed_m("crossed", pc->R.getView().extent(0)),
      t_m("t", pc->R.getView().extent(0)) {}

void OutputPlane::StashedContainer::stash(ParticleContainer_t* pc) {
    Kokkos::deep_copy(r_m, pc->R.getView());
    Kokkos::deep_copy(p_m, pc->P.getView());
    Kokkos::deep_copy(t_m, pc->dt.getView());
    crossed_m = BooleanView_t("crossed", pc->R.getView().extent(0));
    Kokkos::deep_copy(crossed_m, false);
    const auto rView = r_m;
    const auto refR  = pc->getRefPartR();
    Kokkos::parallel_for(
            "OutputPlane::stash", rView.extent(0),
            KOKKOS_LAMBDA(const size_t j) { rView(j) += refR; });
}
