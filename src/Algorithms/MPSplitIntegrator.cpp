// ------------------------------------------------------------------------
// $RCSfile: MPSplitIntegrator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MPSplitIntegrator
//   A MPSplitIntegrator performs integration through an element using two
//   thin lenses of force 1/2, one placed at 1/6 and the other at 5/6 of
//   the length respectively.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/MPSplitIntegrator.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/Multipole.h"
#include "Algorithms/PartBunchBase.h"
#include "Algorithms/PartData.h"
#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/FVps.h"
#include "Physics/Physics.h"


// Class MPSplitIntegrator
// ------------------------------------------------------------------------

MPSplitIntegrator::MPSplitIntegrator(Multipole *mult, int slices):
    MapIntegrator(mult), itsMultipole(mult), itsSlices(slices)
{}


MPSplitIntegrator::MPSplitIntegrator(const MPSplitIntegrator &rhs):
    MapIntegrator(rhs), itsSlices(rhs.itsSlices)
{}


MPSplitIntegrator::~MPSplitIntegrator()
{}


MPSplitIntegrator *MPSplitIntegrator::clone() const {
    return new MPSplitIntegrator(*this);
}


BGeometryBase &MPSplitIntegrator::getGeometry() {
    return itsElement->getGeometry();
}


const BGeometryBase &MPSplitIntegrator::getGeometry() const {
    return itsElement->getGeometry();
}


ElementBase::ElementType MPSplitIntegrator::getType() const {
    return MPSPLITINTEGRATOR;
}


void MPSplitIntegrator::getMap(FVps<double, 6> &map, const PartData &ref,
                               bool revBeam, bool revTrack) const {
    map = FVps<double, 6>();
    trackMap(map, ref, revBeam, revTrack);
}


void MPSplitIntegrator::trackMap(FVps<double, 6> &map,
                                 const PartData &ref,
                                 bool revBeam, bool revTrack) const {
    double length = itsMultipole->getElementLength();
    if(revTrack) length = - length;
    BMultipoleField field = itsMultipole->getField();
    double scale = (ref.getQ() * Physics::c) / (ref.getP());
    if(revBeam) scale = - scale;

    if(length) {
        std::vector<double> slices;
        getSlices(slices);
        scale *= length / double(itsSlices);

        // Drift to first slice position.
        applyDrift(map, length * slices[0], ref);

        for(int i = 0; i < itsSlices; ++i) {
            // Apply first thin multipole kick.
            applyMultipole(map, field, scale);

            // Drift to next slice position or to end.
            applyDrift(map, length * (slices[i+1] - slices[i]), ref);
        }
    } else {
        // Length == 0, slicing not possible.
        applyMultipole(map, field, scale);
    }
}


void MPSplitIntegrator::trackParticle(OpalParticle &part, const PartData &ref,
                                      bool revBeam, bool revTrack) const {
    double length = itsMultipole->getElementLength();
    if(revTrack) length = - length;
    BMultipoleField field = itsMultipole->getField();
    double scale = (ref.getQ() * Physics::c) / (ref.getP());
    if(revBeam) scale = - scale;

    if(length) {
        std::vector<double> slices;
        getSlices(slices);
        scale *= length / double(itsSlices);

        // Drift to first slice position.
        applyDrift(part, length * slices[0], ref);

        for(int i = 0; i < itsSlices; ++i) {
            // Apply first thin multipole kick.
            applyMultipole(part, field, scale);

            // Drift to next slice position or to end.
            applyDrift(part, length * (slices[i+1] - slices[i]), ref);
        }
    } else {
        // Length == 0, slicing not possible.
        applyMultipole(part, field, scale);
    }
}


void MPSplitIntegrator::trackBunch(PartBunchBase<double, 3> *bunch,
                                   const PartData &ref,
                                   bool revBeam, bool revTrack) const {
    double length = itsMultipole->getElementLength();
    if(revTrack) length = - length;
    BMultipoleField field = itsMultipole->getField();
    double scale = (ref.getQ() * Physics::c) / (ref.getP());
    if(revBeam) scale = - scale;

    if(length) {
        std::vector<double> slices;
        getSlices(slices);
        scale *= length / double(itsSlices);
        for(unsigned int i = 0; i < bunch->getLocalNum(); i++) {
            OpalParticle part = bunch->get_part(i);

            // Drift to first slice position.
            applyDrift(part, length * slices[0], ref);

            for(int s = 0; s < itsSlices; ++s) {
                // Apply first thin multipole kick.
                applyMultipole(part, field, scale);

                // Drift to next slice position or to end.
                applyDrift(part, length * (slices[s+1] - slices[s]), ref);
            }
            bunch->set_part(part, i);
        }
    } else {
        // Length == 0, slicing not possible.
        for(unsigned int i = 0; i < bunch->getLocalNum(); i++) {
            OpalParticle part = bunch->get_part(i);
            applyMultipole(part, field, scale);
            bunch->set_part(part, i);
        }
    }
}


void MPSplitIntegrator::applyDrift(FVps<double, 6> &map, double length,
                                   const PartData &reference) const {
    double kin = reference.getM() / reference.getP();
    double ref  = kin * kin;
    FTps<double, 6> px = map[1];
    FTps<double, 6> py = map[3];
    FTps<double, 6> pt = map[5];
    FTps<double, 6> lByPz = length / (1.0 + pt);
    map[0] += px * lByPz;
    map[2] += py * lByPz;
    map[4] += length * (pt * ref - (px * px + py * py + (3.0 * ref) * pt * pt) / 2.0);
}


void MPSplitIntegrator::applyDrift(OpalParticle &part, double length,
                                   const PartData &reference) const {
    double kin = reference.getM() / reference.getP();
    double ref  = kin * kin;
    double px = part.px();
    double py = part.py();
    double pt = part.pt();
    double lByPz = length / (1.0 + pt);
    part.x() += px * lByPz;
    part.y() += py * lByPz;
    part.t() += length * (pt * ref - (px * px + py * py + (3.0 * ref) * pt * pt) / 2.0);
}


void MPSplitIntegrator::applyMultipole(FVps<double, 6> &map,
                                       const BMultipoleField &field,
                                       double scale) const {
    int order = field.order();

    if(order > 0) {
        FTps<double, 6> x = map[0];
        FTps<double, 6> y = map[2];
        FTps<double, 6> kx = + field.normal(order);
        FTps<double, 6> ky = - field.skew(order);

        while(--order > 0) {
            FTps<double, 6> kxt = x * kx - y * ky;
            FTps<double, 6> kyt = x * ky + y * kx;
            kx = kxt + field.normal(order);
            ky = kyt - field.skew(order);
        }

        map[1] -= kx * scale;
        map[3] += ky * scale;
    }
}


void MPSplitIntegrator::applyMultipole(OpalParticle &part,
                                       const BMultipoleField &field,
                                       double scale) const {
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

        part.px() -= kx * scale;
        part.py() += ky * scale;
    }
}


void MPSplitIntegrator::getSlices(std::vector<double> &slices) const {
    slices.clear();
    slices.reserve(itsSlices + 1);

    switch(itsSlices) {

        case 1:
            slices.push_back(1.0 / 2.0);
            break;

        case 2:
            slices.push_back(1.0 / 6.0);
            slices.push_back(5.0 / 6.0);
            break;

        case 3:
            slices.push_back(1.0 / 8.0);
            slices.push_back(4.0 / 8.0);
            slices.push_back(7.0 / 8.0);
            break;

        case 4:
            slices.push_back(3.0 / 30.0);
            slices.push_back(11.0 / 30.0);
            slices.push_back(19.0 / 30.0);
            slices.push_back(27.0 / 30.0);
            break;

        default: {
            double step = 1.0 / double(itsSlices);
            double pos = step / 2.0;
            for(int i = 1; i < itsSlices; ++i) {
                slices.push_back(pos);
                pos += step;
            }
            break;
        }
    }

    slices.push_back(1.0);
}