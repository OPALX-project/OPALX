// ------------------------------------------------------------------------
// $RCSfile: ThickTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThickTracker
//   The visitor class for building a map of given order for a beamline
//   using a finite-length lenses for all elements.
//   Multipole-like elements are done by expanding the Lie series.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include <cfloat>

#include "Algorithms/ThickTracker.h"

#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/CyclotronValley.h"

#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Beamlines/Beamline.h"

#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/FVps.h"

#include "Physics/Physics.h"

#include "Utilities/NumToStr.h"


class Beamline;
class PartData;
using Physics::c;

#define PSdim 6
typedef FVector<double, PSdim> Vector;
typedef FMatrix<double, PSdim, PSdim> Matrix;
typedef FTps<double, PSdim> Series;
typedef FVps<double, PSdim> Map, VSeries;
typedef FMatrix<FTps<double, PSdim>, PSdim, PSdim> MxSeries;

namespace {
    Vector implicitIntStep(const Vector &zin, const VSeries &f, const MxSeries gradf1, double ds,
                           int nx = 20);
    Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds,
                        int nx = 20, int cx = 4);
    Vector fixedPointInt2(const Vector &zin, const VSeries &f, double ds,
                          int nx = 50);
    Vector fixedPointInt4(const Vector &zin, const VSeries &f, double s, double ds,
                          int nx = 50, int cx = 4);
};

// Class ThickTracker
// ------------------------------------------------------------------------

ThickTracker::ThickTracker(const Beamline &beamline,
                           const PartData &reference,
                           bool revBeam, bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack)
{}


ThickTracker::ThickTracker(const Beamline &beamline,
                           const PartBunch &bunch,
                           const PartData &reference,
                           bool revBeam, bool revTrack):
    Tracker(beamline, bunch, reference, revBeam, revTrack)
{}


ThickTracker::~ThickTracker()
{}


void ThickTracker::visitBeamBeam(const BeamBeam &) {
    // *** MISSING *** Tracker for beam-beam.
}


void ThickTracker::visitCollimator(const Collimator &coll) {
    applyDrift(flip_s * coll.getElementLength());
}


void ThickTracker::visitCorrector(const Corrector &corr) {
    // This implements the ThinTracker version.
    double length = flip_s * corr.getElementLength();

    // Drift through first half of length.
    if(length != 0.0) applyDrift(length / 2.0);

    // Apply kick.
    double scale =
        (flip_s * flip_B * itsReference.getQ() * c) / itsReference.getP();
    if(length != 0.0) scale *= length; // Is this right?!?
    const BDipoleField &field = corr.getField();

    for(unsigned int i = 0; i < itsBunch.getLocalNum(); i++) {
        Particle part = itsBunch.get_part(i);
        part.px() -= field.getBy() * scale;
        part.py() += field.getBx() * scale;
        itsBunch.set_part(part, i);
    }

    // Drift through second half of length.
    if(length != 0.0) applyDrift(length / 2.0);
}


void ThickTracker::visitDiagnostic(const Diagnostic &diag) {
    // The diagnostic has no effect on particle tracking.
    applyDrift(flip_s * diag.getElementLength());
}


void ThickTracker::visitDrift(const Drift &drift) {
    applyDrift(flip_s * drift.getElementLength());
}


void ThickTracker::visitLambertson(const Lambertson &lamb) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * lamb.getElementLength());
}


void ThickTracker::visitMarker(const Marker &marker) {
    // Do nothing.
}


void ThickTracker::visitMonitor(const Monitor &corr) {
    applyDrift(flip_s * corr.getElementLength());
}


void ThickTracker::visitMultipole(const Multipole &mult) {
    //std::cerr << "==> In ThickTracker::visitMultipole(const Multipole &mult)" << std::endl;
    // Geometry and Field
    double length = flip_s * mult.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = mult.getField();
    int order = field.order();
    double beta = itsReference.getBeta();

    if(length == 0.0) {
        // Thin multipole, field coefficients are integral(B*dl).
        scale *= flip_s;
        applyThinMultipole(field, scale);
    } else {
        // Finite-length multipole, field coefficients are B.
        // Apply entrance fringe field.
        applyEntranceFringe(0.0, 0.0, field, scale);

        // Construct terms in Hamiltonian H.
        // (1) Define variables.
        static const Series px = Series::makeVariable(PX);
        static const Series py = Series::makeVariable(PY);
        static const Series pt = Series::makeVariable(PT) + 1.0; // 1 + \delta = p/p0
        static const Map ident;
        // (2) Construct kinematic terms.
        double kin = itsReference.getM() / itsReference.getP(); // 1/(beta*gamma)
        Series ptk = (pt * pt + kin * kin);
        Series pz2 = (pt * pt - px * px - py * py);
        // (3) Build vector potential in straight reference frame.
        Series As = buildMultipoleVectorPotential(field) * scale;
        As.setTruncOrder(Series::EXACT);

        // Use implicit integration to propagate particles.
        for(unsigned int i = 0; i < itsBunch.getLocalNum(); i++) {
            Particle part = itsBunch.get_part(i);
            Vector zin;
            zin[X] = part.x();
            zin[PX] = part.px();
            zin[Y] = part.y();
            zin[PY] = part.py();
            zin[T] = part.t();
            zin[PT] = part.pt();
            // Build H_trans and J.grad(H_trans).
            if(zin[X] != DBL_MAX) {
                bool ok = true;
                try {
                    Map translate = zin + ident;
                    Series H_trans = (1 / beta) * sqrt(ptk.substitute(translate), order)
                                     - sqrt(pz2.substitute(translate), order)
                                     + As.substitute(translate);
                    VSeries JgradH;
                    for(int i = 0; i < PSdim; i += 2) {
                        JgradH[ i ] =  H_trans.derivative(i + 1);
                        JgradH[i+1] = -H_trans.derivative(i);
                    }
                    static const Vector zeroV;
                    Vector zf = zin + implicitInt4(zeroV, JgradH, length, 0.5 * length);
                    part.x() = zf[X];
                    part.px() = zf[PX];
                    part.y() = zf[Y];
                    part.py() = zf[PY];
                    part.t() = zf[T];
                    part.pt() = zf[PT];
                } catch(SingularMatrixError &smx) {
                    std::cerr << " [ThickTracker::visitMultipole:] Singular matrix error caught from method " << smx.where() << "\n"
                              << "   " << smx.what() << "\n"
                              << "   Treating particle as lost." << std::endl;
                    ok = false;
                } catch(DomainError &dmx) {
                    std::cerr << " [ThickTracker::visitMultipole:] Domain error caught from method " << dmx.where() << "\n"
                              << "   Treating particle as lost." << std::endl;
                    ok = false;
                } catch(ConvergenceError &cnvx) {
                    std::cerr << " [ThickTracker::visitMultipole:] Convergence error caught from method " << cnvx.where() << "\n"
                              << "   " << cnvx.what() << "\n"
                              << "   Treating particle as lost." << std::endl;
                    ok = false;
                }
                if(!ok) {
                    part.x() = DBL_MAX;
                    part.px() = DBL_MAX;
                    part.y() = DBL_MAX;
                    part.py() = DBL_MAX;
                    part.t() = DBL_MAX;
                    part.pt() = DBL_MAX;
                }
            }
            itsBunch.set_part(part, i);
        }

        // Apply exit fringe field.
        applyExitFringe(0.0, 0.0, field, scale);
    }
    //std::cerr << "==> Leaving ThickTracker::visitMultipole(...)" << std::endl;
}

void ThickTracker::visitProbe(const Probe &prob) {
    applyDrift(flip_s * prob.getElementLength());
}


void ThickTracker::visitRBend(const RBend &bend) {
    //std::cerr << "==> In ThickTracker::visitRBend(const RBend &bend)" << std::endl;
    // Geometry and Field.
    const RBendGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = bend.getField();
    int order = field.order();
    double beta_inv = 1.0 / itsReference.getBeta();

    // Compute slices and stepsize.
    int nSlices = (int) bend.getSlices();
    double stepsize = bend.getStepsize();
    if(stepsize != 0.0) {
        int nst = (int) ceil(length / stepsize);
        nSlices = std::max(nSlices, nst);
    }
    double d_ell = length / nSlices;

    if(length == 0.0) {
        // Thin RBend.
        double half_angle = flip_s * geometry.getBendAngle() / 2.0;
        Euclid3D rotat = Euclid3D::YRotation(- half_angle);
        applyTransform(rotat, 0.0); // Transform from in-plane to mid-plane.
        applyThinMultipole(field, scale); // Apply multipole kick.
        applyTransform(rotat, 0.0); // Transform from mid-plane to out-plane.
    } else {
        // Finite-length RBend.
        if(back_path) {
            // Apply rotation global to local.
            applyTransform(Inverse(geometry.getExitPatch()));
            // Apply entrance fringe field.
            applyEntranceFringe(bend.getExitFaceRotation(), scale, field, scale);
        } else {
            // Apply rotation global to local.
            applyTransform(geometry.getEntrancePatch());
            // Apply exit fringe field.
            applyEntranceFringe(bend.getEntryFaceRotation(), scale, field, scale);
        }

        // Construct terms in Hamiltonian H.
        // (1) Define variables.
        static const Series px  = Series::makeVariable(PX);
        static const Series py  = Series::makeVariable(PY);
        static const Series de  = Series::makeVariable(PT);
        static const Series de2 = de * de;
        static const Series pxy2 = px * px + py * py;
        static const Map ident;
        // (2) Construct kinematic terms.
        double kin = itsReference.getM() / itsReference.getP(); // mc^2/pc = 1/(beta*gamma)
        double kin2 = kin * kin;
        // (3) Vector potential in straight reference.
        Series As = buildMultipoleVectorPotential(field) * scale;
        As.setTruncOrder(Series::EXACT);

        // Use implicit integration to propagate particles.

        for(unsigned int i = 0; i < itsBunch.getLocalNum(); i++) {
            Particle part = itsBunch.get_part(i);
            Vector zin, zf;
            zin[X] = part.x();
            zin[PX] = part.px();
            zin[Y] = part.y();
            zin[PY] = part.py();
            zin[T] = part.t();
            zin[PT] = part.pt();
            if(zin[X] != DBL_MAX) {
                bool ok = true;
                try {
                    for(int slice = 0; slice < nSlices; ++slice) {
                        // Build H_trans and J.grad(H_trans).
                        Map translate = zin + ident;
                        double d1 = 1.0 + zin[PT];
                        Series de12 = d1 * d1 + 2.0 * d1 * de + de2;
                        Series H_trans = beta_inv * sqrt(de12 + kin2, order)
                                         - sqrt(de12 - zin[PX] * zin[PX] - zin[PY] * zin[PY] - 2.0 * (zin[PX] * px + zin[PY] * py) - pxy2, order)
                                         + As.substitute(translate);
                        VSeries JgradH;
                        for(int i = 0; i < PSdim; i += 2) {
                            JgradH[ i ] =  H_trans.derivative(i + 1);
                            JgradH[i+1] = -H_trans.derivative(i);
                        }
                        // Propagate particle.
                        static const Vector zeroV;
                        zf = zin + implicitInt4(zeroV, JgradH, d_ell, 0.5 * d_ell);
                        //            zf = zin;
                        //            std::cerr << " [visitRBend:] zin = " << zin;
                        //            std::cerr << "                dl = " << d_ell << std::endl;
                        //            std::cerr << "                zf = " << zf     << std::endl;
                        zin = zf;
                    }
                } catch(SingularMatrixError &smx) {
                    std::cerr << " [ThickTracker::visitRBend:] Singular matrix error caught; treating particle as lost." << std::endl;
                    ok = false;
                } catch(DomainError &dmx) {
                    std::cerr << " [ThickTracker::visitRBend:] Domain error caught; treating particle as lost." << std::endl;
                    ok = false;
                } catch(ConvergenceError &cnvx) {
                    std::cerr << " [ThickTracker::visitRBend:] Convergence error caught; treating particle as lost." << std::endl;
                    ok = false;
                }
                if(ok) {
                    part.x() = zf[X];
                    part.px() = zf[PX];
                    part.y() = zf[Y];
                    part.py() = zf[PY];
                    part.t() = zf[T];
                    part.pt() = zf[PT];
                } else {
                    part.x() = DBL_MAX;
                    part.px() = DBL_MAX;
                    part.y() = DBL_MAX;
                    part.py() = DBL_MAX;
                    part.t() = DBL_MAX;
                    part.pt() = DBL_MAX;
                }
            }
            itsBunch.set_part(part, i);
        }

        if(back_path) {
            // Apply exit fringe field.
            applyExitFringe(bend.getEntryFaceRotation(), scale, field, scale);
            // Transform from local to global.
            applyTransform(Inverse(geometry.getEntrancePatch()));
        } else {
            // Apply entrance fringe field.
            applyEntranceFringe(bend.getExitFaceRotation(), scale, field, scale);
            // Transform from local to global.
            applyTransform(geometry.getExitPatch());
        }
    }
    //std::cerr << "==> Leaving ThickTracker::visitRBend(...)" << std::endl;
}


void ThickTracker::visitRFCavity(const RFCavity &as) {
    // This implements the ThinTracker version.
    double length = flip_s * as.getElementLength();

    // Drift through half length.
    if(length != 0.0) applyDrift(length / 2.0);

    // Apply accelerating voltage.
    double freq = as.getFrequency();
    double kin = itsReference.getM() / itsReference.getP();
    double peak = flip_s * as.getAmplitude() / itsReference.getP();

    for(unsigned int i = 0; i < itsBunch.getLocalNum(); i++) {
        Particle part = itsBunch.get_part(i);
        double pt    = (part.pt() + 1.0);
        double speed = (c * pt) / sqrt(pt * pt + kin * kin);
        double phase = as.getPhase() + (freq * part.t()) / speed;
        part.pt() += peak * sin(phase) / pt;
        itsBunch.set_part(part, i);
    }

    // Drift through half length.
    if(length != 0.0) applyDrift(length / 2.0);
}


void ThickTracker::visitRFQuadrupole(const RFQuadrupole &rfq) {
    // *** MISSING *** Tracking for RF quadrupole.
    applyDrift(flip_s * rfq.getElementLength());
}


void ThickTracker::visitSBend(const SBend &bend) {
    //std::cerr << "==> In ThickTracker::visitSBend(const SBend &bend)" << std::endl;
    // Geometry and Field.
    const PlanarArcGeometry &geometry = bend.getGeometry();
    double length = flip_s * geometry.getElementLength();
    double scale = (flip_B * itsReference.getQ() * c) / itsReference.getP();
    const BMultipoleField &field = bend.getField();
    int order = field.order();
    double beta_inv = 1.0 / itsReference.getBeta();

    // Compute slices and stepsize.
    int nSlices = (int) bend.getSlices();
    double stepsize = bend.getStepsize();
    if(stepsize != 0.0) {
        int nst = (int) ceil(length / stepsize);
        nSlices = std::max(nSlices, nst);
    }
    double d_ell = length / nSlices;

    if(length == 0.0) {
        // Thin SBend.
        double half_angle = geometry.getBendAngle() / 2.0;
        Euclid3D rotat = Euclid3D::YRotation(- half_angle);
        applyTransform(rotat, 0.0); // Transform from in-plane to mid-plane.
        applyThinSBend(field, scale, 0.0); // Apply multipole kick.
        applyTransform(rotat, 0.0); // Transform from mid-plane to out-plane.
    } else {
        // Finite-length SBend.
        // Apply entrance fringe field.
        applyEntranceFringe(bend.getEntryFaceRotation(),
                            bend.getEntryFaceCurvature(), field, scale);

        // Construct terms in Hamiltonian H.
        // (1) Define variables.
        static const Series x   = Series::makeVariable(X);
        static const Series px  = Series::makeVariable(PX);
        static const Series py  = Series::makeVariable(PY);
        static const Series de  = Series::makeVariable(PT);
        static const Series de2 = de * de;
        static const Series pxy2 = px * px + py * py;
        static const Map ident;
        // (2) Construct kinematic terms.
        double h = geometry.getCurvature();
        double kin = itsReference.getM() / itsReference.getP(); // mc^2/pc = 1/(beta*gamma)
        double kin2 = kin * kin;
        Series hx1 = (1.0 + h * x);
        // (3) Build vector potential, actually -(1+h*x)*As in reference frame of curvature h.
        Series As = buildSBendVectorPotential(field, h) * scale;
        As.setTruncOrder(Series::EXACT);

        // Use implicit integration to propagate particles.


        for(unsigned int i = 0; i < itsBunch.getLocalNum(); i++) {
            Particle part = itsBunch.get_part(i);

            Vector zin, zf;
            zin[X] = part.x();
            zin[PX] = part.px();
            zin[Y] = part.y();
            zin[PY] = part.py();
            zin[T] = part.t();
            zin[PT] = part.pt();
            if(zin[X] != DBL_MAX) {
                bool ok = true;
                try {
                    for(int slice = 0; slice < nSlices; ++slice) {
                        // Build H_trans and J.grad(H_trans).
                        Map translate = zin + ident;
                        double d1 = 1.0 + zin[PT];
                        Series de12 = d1 * d1 + 2.0 * d1 * de + de2;
                        Series H_trans = beta_inv * sqrt(de12 + kin2, order)
                                         - (hx1 + h * zin[X])
                                         * sqrt(de12 - zin[PX] * zin[PX] - zin[PY] * zin[PY] - 2.0 * (zin[PX] * px + zin[PY] * py) - pxy2, order)
                                         + As.substitute(translate);
                        VSeries JgradH;
                        for(int i = 0; i < PSdim; i += 2) {
                            JgradH[ i ] =  H_trans.derivative(i + 1);
                            JgradH[i+1] = -H_trans.derivative(i);
                        }
                        // Propagate particle.
                        static const Vector zeroV;
                        zf = zin + fixedPointInt4(zeroV, JgradH, d_ell, 0.5 * d_ell);
                        //zf = zin + implicitInt4(zeroV,JgradH,d_ell,0.5*d_ell);
                        //            zf = zin;
                        //            std::cerr << " [visitSBend:] zin = " << zin;
                        //            std::cerr << "                dl = " << d_ell << std::endl;
                        //            std::cerr << "                zf = " << zf     << std::endl;
                        zin = zf;
                    }
                } catch(SingularMatrixError &smx) {
                    std::cerr << " [ThickTracker::visitSBend:] Singular matrix error caught from method " << smx.where() << "\n"
                              << "   " << smx.what() << "\n"
                              << "   Treating particle as lost." << std::endl;
                    ok = false;
                } catch(DomainError &dmx) {
                    std::cerr << " [ThickTracker::visitSBend:] Domain error caught from method " << dmx.where() << "\n"
                              "   Treating particle as lost." << std::endl;
                    ok = false;
                } catch(ConvergenceError &cnvx) {
                    std::cerr << " [ThickTracker::visitSBend:] Convergence error caught from method " << cnvx.where() << "\n"
                              << "   " << cnvx.what() << "\n"
                              << "   Treating particle as lost." << std::endl;
                    ok = false;
                }
                if(ok) {
                    part.x() = zf[X];
                    part.px() = zf[PX];
                    part.y() = zf[Y];
                    part.py() = zf[PY];
                    part.t() = zf[T];
                    part.pt() = zf[PT];
                } else {
                    part.x() = DBL_MAX;
                    part.px() = DBL_MAX;
                    part.y() = DBL_MAX;
                    part.py() = DBL_MAX;
                    part.t() = DBL_MAX;
                    part.pt() = DBL_MAX;
                }
            }
            itsBunch.set_part(part, i);
        }

        // Apply exit fringe field.
        applyExitFringe(bend.getExitFaceRotation(),
                        bend.getExitFaceCurvature(), field, scale);
    }
    //std::cerr << "==> Leaving ThickTracker::visitSBend(...)" << std::endl;
}


void ThickTracker::visitSeparator(const Separator &sep) {
    // This implements the ThinTracker version.
    double length = flip_s * sep.getElementLength();

    if(length != 0.0) {
        // Drift through first half of length.
        applyDrift(length / 2.0);

        // Apply kick.
        double scale = (length * itsReference.getQ()) / itsReference.getP();
        double Ex = scale * sep.getEx();
        double Ey = scale * sep.getEy();

        for(unsigned int i = 0; i < itsBunch.getLocalNum(); i++) {
            Particle part = itsBunch.get_part(i);
            double pt = 1.0 + part.pt();
            part.px() += Ex / pt;
            part.py() += Ey / pt;
            itsBunch.set_part(part, i);
        }

        // Drift through second half of length.
        applyDrift(length / 2.0);
    }
}


void ThickTracker::visitSeptum(const Septum &sept) {
    // Assume the particle go through the magnet's window.
    applyDrift(flip_s * sept.getElementLength());
}


void ThickTracker::visitSolenoid(const Solenoid &solenoid) {
    // This needs to be checked!
    double length = flip_s * solenoid.getElementLength();

    if(length != 0.0) {
        double ks = (flip_B * itsReference.getQ() * solenoid.getBz() * c) /
                    (2.0 * itsReference.getP());

        if(ks) {
            double kin = itsReference.getM() / itsReference.getP();
            double refTime = length / itsReference.getBeta();


            for(unsigned int i = 0; i < itsBunch.getLocalNum(); i++) {
                Particle part = itsBunch.get_part(i);
                double pt = part.pt() + 1.0;
                double px = part.px() + ks * part.y();
                double py = part.py() - ks * part.x();
                double pz = sqrt(pt * pt - px * px - py * py);
                double E = sqrt(pt * pt + kin * kin);

                double k = ks / pz;
                double C = cos(k * length);
                double S = sin(k * length);

                double xt  = C * part.x()  + S * part.y();
                double yt  = C * part.y()  - S * part.x();
                double pxt = C * part.px() + S * part.py();
                double pyt = C * part.py() - S * part.px();

                part.x()  = C * xt  + (S / k) * pxt;
                part.y()  = C * yt  + (S / k) * pyt;
                part.px() = C * pxt - (S * k) * xt;
                part.py() = C * pyt - (S * k) * yt;
                part.t() += pt * (refTime / E - length / pz);
                itsBunch.set_part(part, i);
            }
        } else {
            applyDrift(length);
        }
    }
}


void ThickTracker::visitParallelPlate(const ParallelPlate &pplate) {
    //do nothing
}

void ThickTracker::visitCyclotronValley(const CyclotronValley &cv) {
    // Do nothing here.
}


void ThickTracker::applyEntranceFringe(double angle, double curve,
                                       const BMultipoleField &field, double scale) {
    // *** MISSING *** Higher order terms for entrance fringe.
    double ca = cos(angle);
    double sa = sin(angle);
    double ta = tan(angle);

    int order = field.order();
    static const Series x = Series::makeVariable(X);
    Series by = field.normal(order);
    for(int i = order; --i >= 1;) by = by * x + field.normal(i);
    by *= scale;

    //Track particles.

    for(unsigned int i = 0; i < itsBunch.getLocalNum(); i++) {
        Particle part = itsBunch.get_part(i);
        // Read particle coordinates.
        double x = part.x();
        double px = part.px();
        double y = part.y();
        double py = part.py();
        double t = part.t();
        double pt = part.pt();
        if(x != DBL_MAX) {
            double p = pt + 1.0;
            // rotate to magnet face
            double ps = sqrt(p * p - px * px - py * py);
            x = x / (ca * (1.0 - ta * px / ps));
            px = ca * px + sa * ps;
            double ellovpp = sa * x / ps;
            y += ellovpp * py;
            t -= ellovpp * p;
            // fringe
            Vector z;
            z[X] = x;
            z[PX] = px;
            z[Y] = y;
            z[PY] = py;
            z[T] = t;
            z[PT] = pt;
            double psy = sqrt(p * p - px * px);
            py -= by.evaluate(z) * y * px / psy;
            // rotate from magnet face
            ps = sqrt(p * p - px * px - py * py);
            x = x / (ca * (1.0 + ta * px / ps));
            px = ca * px - sa * ps;
            ellovpp = sa * x / ps;
            y -= ellovpp * py;
            t += ellovpp * p;
            // edge effect
            z[X] = x;
            z[PX] = px;
            z[Y] = y;
            z[PY] = py;
            z[T] = t;  //z[PT] = pt;
            px += by.evaluate(z) * ta * x;
            // Load particle coordinates.
            part.x() = x;
            part.px() = px;
            part.y() = y;
            part.py() = py;
            part.t() = t;
            part.pt() = pt;
        }
        itsBunch.set_part(part, i);
    }
}


void ThickTracker::applyExitFringe(double angle, double curve,
                                   const BMultipoleField &field, double scale) {
    // *** MISSING *** Higher order terms for exit fringe.
    double ca = cos(angle);
    double sa = sin(angle);
    double ta = tan(angle);

    int order = field.order();
    static const Series x = Series::makeVariable(X);
    Series by = field.normal(order);
    for(int i = order; --i >= 1;) by = by * x + field.normal(i);
    by *= scale;

    //Track particles.
    for(unsigned int i = 0; i < itsBunch.getLocalNum(); i++) {
        Particle part = itsBunch.get_part(i);
        // Read particle coordinates.
        double x = part.x();
        double px = part.px();
        double y = part.y();
        double py = part.py();
        double t = part.t();
        double pt = part.pt();
        if(x != DBL_MAX) {
            double p = pt + 1.0;
            // edge effect
            Vector z;
            z[X] = x;
            z[PX] = px;
            z[Y] = y;
            z[PY] = py;
            z[T] = t;
            z[PT] = pt;
            px += by.evaluate(z) * ta * x;
            // rotate to magnet face
            double ps = sqrt(p * p - px * px - py * py);
            x = x / (ca * (1.0 + ta * px / ps));
            px = ca * px - sa * ps;
            double ellovpp = sa * x / ps;
            y -= ellovpp * py;
            t += ellovpp * p;
            // fringe
            z[X] = x;
            z[PX] = px;
            z[Y] = y;  //z[PY] = py;
            z[T] = t;  //z[PT] = pt;
            double psy = sqrt(p * p - px * px);
            py += by.evaluate(z) * y * px / psy;
            // rotate from magnet face
            ps = sqrt(p * p - px * px - py * py);
            x = x / (ca * (1.0 - ta * px / ps));
            px = ca * px + sa * ps;
            ellovpp = sa * x / ps;
            y += ellovpp * py;
            t -= ellovpp * p;
            // Load particle coordinates.
            part.x() = x;
            part.px() = px;
            part.y() = y;
            part.py() = py;
            part.t() = t;
            part.pt() = pt;
        }
        itsBunch.set_part(part, i);
    }
}

namespace {
    Vector implicitInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx) {
        //std::cerr << "==> In implicitInt4(zin,f,s,ds,nx,cx) ..." << std::endl;
        // Default: nx = 20, cx = 4

        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a distance s, in steps of size ds.  It uses a "Yoshida-fied"
        // version of implicitInt2 to obtain zf accurate through fourth-order in the
        // step-size ds.  When f derives from a Hamiltonian---i.e., f = J.grad(H)---
        // then this routine performs symplectic integration.  The optional arguments
        // nx and cx have the same meaning as in implicitInt2().

        // Convergence warning flag.
        static bool cnvWarn = false;

        // The Yoshida constants: 2ya+yb=1; 2ya^3+yb^3=0.
        static const double yt = pow(2., 1 / 3.);
        static const double ya = 1 / (2. - yt);
        static const double yb = -yt * ya;

        // Build matrix grad(f).
        MxSeries gradf;
        for(int i = 0; i < PSdim; ++i)
            for(int j = 0; j < PSdim; ++j)
                gradf[i][j] = f[i].derivative(j);

        // Initialize accumulated length, current step-size, and number of cuts.
        double as = std::abs(s), st = 0., dsc = std::abs(ds);
        if(s < 0.) dsc = -dsc;
        int ci = 0;

        // Integrate each step.
        Vector zf = zin;
        while(std::abs(st) < as) {
            Vector zt;
            bool ok = true;
            try {
                if(std::abs(st + dsc) > as) dsc = s - st;
                zt = ::implicitIntStep(zf, f, gradf, ya * dsc, nx);
                zt = ::implicitIntStep(zt, f, gradf, yb * dsc, nx);
                zt = ::implicitIntStep(zt, f, gradf, ya * dsc, nx);
            } catch(ConvergenceError &cnverr) {
                if(++ci > cx) {
                    string msg = "Convergence not achieved within " + NumToStr(cx) + " cuts of step-size!";
                    throw ConvergenceError("ThickTracker::implicitInt4()", msg);
                }
                if(!cnvWarn) {
                    std::cerr << " <***WARNING***> [ThickTracker::implicitInt4()]:\n"
                              << "   Cutting step size, a probable violation of the symplectic condition."
                              << std::endl;
                    cnvWarn = true;
                }
                dsc *= 0.5;
                ok = false;
            }
            if(ok) {zf = zt; st += dsc;}
        }

        //std::cerr << "==> Leaving implicitInt4(...)" << std::endl;
        return zf;
    }

    Vector implicitIntStep(const Vector &zin, const VSeries &f, const MxSeries gradf, double ds, int nx) {
        //std::cerr << "==> In implicitIntStep(zin,f,gradf,ds,nx) ..." << std::endl;
        //std::cerr << "  ds = " << ds << std::endl;
        //std::cerr << " zin =\n" << zin << std::endl;
        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a single step of size ds, using Newton's method to solve the
        // implicit equation zf = zin + ds*f((zin+zf)/2).  For reasons of efficiency,
        // its arguments include the matrix gradf = grad(f).  The (optional) argument
        // nx limits the number of Newton iterations.  This routine returns a result
        // zf accurate through second-order in the step-size ds.  When f derives from
        // a Hamiltonian---i.e., f=J.grad(H)---then this routine performs symplectic
        // integration.

        // Set up flags, etc., for convergence (bounce) test.
        FVector<bool, PSdim> bounce(false);
        Vector dz, dz_old;
        int bcount = 0;
        static const double thresh = 1.e-8;

        // Use second-order Runge-Kutta integration to determine a good initial guess.
        double ds2 = 0.5 * ds;
        Vector z = f.constantTerm(zin);
        z = zin + ds2 * (z + f.constantTerm(zin + ds * z));

        // Newton iterations:
        //   z :-> [I-ds/2.grad(f)]^{-1}.[zin+ds.f((zin+z)/2)-ds/2.grad(f).z]
        // (A possible method for speeding up this computation would
        //  be to recompute grad(f) every n-th step, where n > 1!)
        Vector zf;
        int ni = 0;
        while(bcount < PSdim) {
            if(ni == nx) {
                string msg = "Convergence not achieved within " + NumToStr(nx) + " iterations!";
                throw ConvergenceError("ThickTracker::implicitIntStep()", msg);
            }

            // Build gf = -ds/2.grad(f)[(zin+z)/2] and idgf_inv = [I-ds/2.grad(f)]^{-1}[(zin+z)/2].
            Vector zt = 0.5 * (zin + z);
            Matrix gf, idgf, idgf_inv;
            for(int i = 0; i < PSdim; ++i)
                for(int j = 0; j < PSdim; ++j)
                    gf[i][j] = -ds2 * gradf[i][j].evaluate(zt);
            idgf = gf;
            for(int i = 0; i < PSdim; ++i) idgf[i][i] += 1.;
            FLUMatrix<double, PSdim> lu(idgf);
            idgf_inv = lu.inverse();

            // Execute Newton step.
            zf = idgf_inv * (zin + ds * f.constantTerm(zt) + gf * z);

            //std::cerr << " -(ds/2)grad(f) =\n" << gf << std::endl;
            //std::cerr << " f =\n" << f.constantTerm(zt) << std::endl;
            //std::cerr << "zk =\n" << zf << std::endl;

            // Test for convergence ("bounce" test).
            dz_old = dz;
            dz = zf - z;
            if(ni) { // (we need at least two iterations before testing makes sense)
                for(int i = 0; i < PSdim; ++i) {
                    if(!bounce[i] && (dz[i] == 0. || (std::abs(dz[i]) < thresh && std::abs(dz[i]) >= std::abs(dz_old[i]))))
                        {bounce[i] = true; ++bcount;}
                }
            }
            z = zf;
            ++ni;
        }

        //std::cerr << "  zf =\n" << zf << std::endl;
        //std::cerr << "==> Leaving implicitIntStep(zin,f,gradf,ds,nx)" << std::endl;
        return zf;
    }

    Vector fixedPointInt2(const Vector &zin, const VSeries &f, double ds, int nx) {
        //std::cerr << "==> In fixedPointInt2(zin,f,ds,nx) ..." << std::endl;
        // Default: nx = 50

        //std::cerr << "  ds = " << ds << std::endl;
        //std::cerr << " zin =\n" << zin << std::endl;
        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a single step of size ds by iterating the equation
        //         z = zin + ds * f((zin+z)/2)
        // to find a fixed-point zf for z.  It is accurate through second order in the
        // step size ds.

        // Set up flags, etc., for convergence (bounce) test.
        FVector<bool, PSdim> bounce(false);
        Vector dz, dz_old;
        int bcount = 0;
        static const double thresh = 1.e-8;

        // Iterate z :-> zin + ds * f( (zin + z)/2 ).
        Vector zf;
        Vector z = zin;
        int ni = 0;
        while(bcount < PSdim) {
            if(ni == nx) {
                string msg = "Convergence not achieved within " + NumToStr(nx) + " iterations!";
                throw ConvergenceError("ThickTracker::fixedPointInt2()", msg);
            }

            // Do iteration.
            zf = zin + ds * f.constantTerm((zin + z) / 2.0);

            // Test for convergence.
            dz_old = dz;
            dz = zf - z;
            if(ni) { // (we need at least two iterations before testing makes sense)
                for(int i = 0; i < PSdim; ++i) {
                    if(!bounce[i] && (dz[i] == 0. || (std::abs(dz[i]) < thresh && std::abs(dz[i]) >= std::abs(dz_old[i]))))
                        {bounce[i] = true; ++bcount;}
                }
            }
            z = zf;
            ++ni;
        }
        //std::cerr << "  zf =\n" << zf << std::endl;
        //std::cerr << "==> Leaving fixedPointInt2(...)" << std::endl;
        return zf;
    }

    Vector fixedPointInt4(const Vector &zin, const VSeries &f, double s, double ds, int nx, int cx) {
        //std::cerr << "==> In fixedPointInt4(zin,f,s,ds,nx,cx) ..." << std::endl;
        // Default: nx = 50, cx = 4

        // This routine integrates the N-dimensional autonomous differential equation
        // z' = f(z) for a distance s, in steps of size ds.  It uses a "Yoshida-fied"
        // version of fixedPointInt2 to obtain zf accurate through fourth-order in the
        // step-size ds.  The optional arguments nx and cx have the same meaning as in
        // implicitInt2().

        // Convergence warning flag.
        static bool cnvWarn = false;

        // The Yoshida constants: 2ya+yb=1; 2ya^3+yb^3=0.
        static const double yt = pow(2., 1 / 3.);
        static const double ya = 1 / (2. - yt);
        static const double yb = -yt * ya;

        // Initialize accumulated length, current step-size, and number of cuts.
        double as = std::abs(s), st = 0., dsc = std::abs(ds);
        if(s < 0.) dsc = -dsc;
        int ci = 0;

        // Integrate each step.
        Vector zf = zin;
        while(std::abs(st) < as) {
            Vector zt;
            bool ok = true;
            try {
                if(std::abs(st + dsc) > as) dsc = s - st;
                zt = ::fixedPointInt2(zf, f, ya * dsc, nx);
                zt = ::fixedPointInt2(zt, f, yb * dsc, nx);
                zt = ::fixedPointInt2(zt, f, ya * dsc, nx);
            } catch(ConvergenceError &cnverr) {
                if(++ci > cx) {
                    string msg = "Convergence not achieved within " + NumToStr(cx) + " cuts of step-size!";
                    throw ConvergenceError("ThickTracker::fixedPointInt4()", msg);
                }
                if(!cnvWarn) {
                    std::cerr << " <***WARNING***> [ThickTracker::fixedPointInt4()]:\n"
                              << "   Cutting step size, a probable violation of the symplectic condition."
                              << std::endl;
                    cnvWarn = true;
                }
                dsc *= 0.5;
                ok = false;
            }
            if(ok) {zf = zt; st += dsc;}
        }

        //std::cerr << "==> Leaving fixedPointInt4(...)" << std::endl;
        return zf;
    }
};
