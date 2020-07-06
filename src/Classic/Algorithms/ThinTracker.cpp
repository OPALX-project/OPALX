// ------------------------------------------------------------------------
// $RCSfile: ThinTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThinTracker
//   The visitor class for tracking a bunch of particles through a beamline
//   using a thin-lens approximation for all elements.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/ThinTracker.h"
#include "AbsBeamline/BeamStripping.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"

#include "Algorithms/PartBunchBase.h"
#include "Algorithms/PartData.h"
#include "Algorithms/OpalParticle.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/Geometry.h"
#include "Beamlines/Beamline.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/ComplexErrorFun.h"
#include <cmath>
#include <complex>


// Class ThinTracker
// ------------------------------------------------------------------------

ThinTracker::ThinTracker(const Beamline &beamline, const PartData &reference,
                         bool revBeam, bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack)
{}


ThinTracker::ThinTracker(const Beamline &beamline,
                         PartBunchBase<double, 3> *bunch,
                         const PartData &reference,
                         bool revBeam, bool revTrack):
    Tracker(beamline, bunch, reference, revBeam, revTrack)
{}


ThinTracker::~ThinTracker()
{}


void ThinTracker::visitBeamStripping(const BeamStripping &/*bstp*/) {
//    applyDrift(flip_s * bstp.getElementLength());
}

void ThinTracker::visitCCollimator(const CCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}

void ThinTracker::visitDegrader(const Degrader &deg) {
    applyDrift(flip_s * deg.getElementLength());
}

void ThinTracker::visitCorrector(const Corrector &corr) {
    // Drift through first half of length.
    double length = flip_s * corr.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Apply kick.
    double scale = (flip_s * flip_B * corr.getElementLength() *
                    itsReference.getQ() * Physics::c) / itsReference.getP();
    const BDipoleField &field = corr.getField();

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        OpalParticle part = itsBunch_m->get_part(i);
        part.px() -= field.getBy() * scale;
        part.py() += field.getBx() * scale;
        itsBunch_m->set_part(part, i);
    }

    // Drift through second half of length.
    if(length) applyDrift(length / 2.0);
}


void ThinTracker::visitDrift(const Drift &drift) {
    applyDrift(flip_s * drift.getElementLength());
}

void ThinTracker::visitFlexibleCollimator(const FlexibleCollimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}

void ThinTracker::visitMarker(const Marker &) {
    // Do nothing.
}


void ThinTracker::visitMonitor(const Monitor &corr) {
    applyDrift(flip_s * corr.getElementLength());
}


void ThinTracker::visitMultipole(const Multipole &mult) {
    double length = flip_s * mult.getElementLength();
    const BMultipoleField &field = mult.getField();
    double scale = (flip_B * itsReference.getQ() * Physics::c) / itsReference.getP();

    if(length) {
        // Drift through first half of the length.
        applyDrift(length / 2.0);

        // Apply multipole kick, field is per unit length.
        scale *= length;
        applyThinMultipole(field, scale);

        // Drift through second half of the length, field is per unit length.
        applyDrift(length / 2.0);
    } else {
        // Thin multipole, field is integral(K*dl).
        scale *= flip_s;
        applyThinMultipole(field, scale);
    }
}


void ThinTracker::visitProbe(const Probe &) {
    // Do nothing.
}


void ThinTracker::visitRBend(const RBend &bend) {
    // Geometry.
    const RBendGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double angle = flip_s * geometry.getBendAngle();

    // Magnetic field.
    const BMultipoleField &field = bend.getField();

    // Drift to mid-plane.
    applyDrift(length / 2.0);

    // Apply multipole kick and linear approximation to geometric bend.
    double scale = (flip_B * itsReference.getQ() * Physics::c) / itsReference.getP();
    if(length) scale *= length;

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        OpalParticle part = itsBunch_m->get_part(i);
        int order = field.order();

        if(order > 0) {
            double x = part.x();
            double y = part.y();
            double kx = + field.normal(order);
            double ky = - field.skew(order);

            while(--order > 0) {
                double kxt = x * kx - y * ky;
                double kyt = x * ky + y * kx;
                kx = kxt + field.normal(order);
                ky = kyt - field.skew(order);
            }

            part.px() -= kx * scale - angle * (1.0 + part.pt());
            part.py() += ky * scale;
            part.t()  -= angle * x;
        }
        itsBunch_m->set_part(part, i);
    }

    // Drift to out-plane.
    applyDrift(length / 2.0);
}


void ThinTracker::visitRFCavity(const RFCavity &as) {
    // Drift through half length.
    double length = flip_s * as.getElementLength();
    if(length) applyDrift(length / 2.0);

    // Apply accelerating voltage.
    double freq = as.getFrequency();
    double peak = flip_s * as.getAmplitude() / itsReference.getP();
    double kin = itsReference.getM() / itsReference.getP();

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        OpalParticle part = itsBunch_m->get_part(i);
        double pt    = (part.pt() + 1.0);
        double speed = (Physics::c * pt) / std::sqrt(pt * pt + kin * kin);
        double phase = as.getPhase() + (freq * part.t()) / speed;
        part.pt() += peak * std::sin(phase) / pt;
        itsBunch_m->set_part(part, i);
    }

    if(length) applyDrift(length / 2.0);
}


void ThinTracker::visitSBend(const SBend &bend) {
    const PlanarArcGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double angle = flip_s * geometry.getBendAngle();

    // Magnetic field.
    const BMultipoleField &field = bend.getField();

    // Drift to mid-plane.
    applyDrift(length / 2.0);

    // Apply multipole kick and linear approximation to geometric bend.
    double scale = (flip_B * itsReference.getQ() * Physics::c) / itsReference.getP();
    if(length) scale *= length;

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        OpalParticle part = itsBunch_m->get_part(i);
        int order = field.order();

        if(order > 0) {
            double x = part.x();
            double y = part.y();
            double kx = + field.normal(order);
            double ky = - field.skew(order);

            while(--order > 0) {
                double kxt = x * kx - y * ky;
                double kyt = x * ky + y * kx;
                kx = kxt + field.normal(order);
                ky = kyt - field.skew(order);
            }

            part.px() -= kx * scale - angle * (1.0 + part.pt());
            part.py() += ky * scale;
            part.t()  -= angle * x;
        }
        itsBunch_m->set_part(part, i);
    }

    // Drift to out-plane.
    applyDrift(length / 2.0);
}


void ThinTracker::visitSeptum(const Septum &sept) {
    // Assume the particles go through the magnet's window.
    applyDrift(flip_s * sept.getElementLength());
}


void ThinTracker::visitSolenoid(const Solenoid &solenoid) {
    double length = flip_s * solenoid.getElementLength();

    if(length) {
        double ks = (flip_B * itsReference.getQ() * solenoid.getBz() * Physics::c) /
                    (2.0 * itsReference.getP());

        if(ks) {
            double kin = itsReference.getM() / itsReference.getP();
            double ref = kin * kin;

            for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
                OpalParticle part = itsBunch_m->get_part(i);

                double pt = part.pt() + 1.0;
                double px = part.px() + ks * part.y();
                double py = part.py() - ks * part.x();
                double pz = std::sqrt(pt * pt - px * px - py * py);

                double k = ks / pz;
                double C = std::cos(k * length);
                double S = std::sin(k * length);

                double xt  = C * part.x()  + S * part.y();
                double yt  = C * part.y()  - S * part.x();
                double pxt = C * part.px() + S * part.py();
                double pyt = C * part.py() - S * part.px();

                part.x()  = C * xt  + (S / k) * pxt;
                part.y()  = C * yt  + (S / k) * pyt;
                part.px() = C * pxt - (S * k) * xt;
                part.py() = C * pyt - (S * k) * yt;
                part.t() += length * (pt * ref - (px * px + py * py + 3.0 * pt * pt * ref) / 2.0);
                itsBunch_m->set_part(part, i);
            }
        } else {
            applyDrift(length);
        }
    }
}


void ThinTracker::applyDrift(double length) {
    double   kin  = itsReference.getM() / itsReference.getP();
    double   ref  = kin * kin;

    for(unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        OpalParticle part = itsBunch_m->get_part(i);

        double px = part.px();
        double py = part.py();
        double pt = part.pt();
        double lByPz = length / (1.0 + pt);
        part.x() += px * lByPz;
        part.y() += py * lByPz;
        part.t() += length * (pt * ref - (px * px + py * py + 3.0 * pt * pt * ref) / 2.0);
        itsBunch_m->set_part(part, i);
    }
}
