/**
 * @class TravelingWave
 * @brief Interface for traveling wave cavities.
 *
 * Class TravelingWave defines the abstract interface for Traveling Wave RF cavities.
 *   - TW: Traveling Wave Cavity
 */
#include "AbsBeamline/TravelingWave.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Astra1DDynamic.h"
#include "Fields/Fieldmap.h"
#include "PartBunch/PartBunch.h"
#include "Physics/Units.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/GSLSpline.h"

#include <fstream>
#include <iostream>

extern Inform* gmsg;

TravelingWave::TravelingWave() : TravelingWave("") {
}

TravelingWave::TravelingWave(const TravelingWave& right)
    : RFCavity(right),
      scaleCore_m(right.scaleCore_m),
      scaleCoreError_m(right.scaleCoreError_m),
      phaseCore1_m(right.phaseCore1_m),
      phaseCore2_m(right.phaseCore2_m),
      phaseExit_m(right.phaseExit_m),
      startCoreField_m(right.startCoreField_m),
      startExitField_m(right.startExitField_m),
      mappedStartExitField_m(right.mappedStartExitField_m),
      periodLength_m(right.periodLength_m),
      numCells_m(right.numCells_m),
      cellLength_m(right.cellLength_m),
      mode_m(right.mode_m) {
}

TravelingWave::TravelingWave(const std::string& name)
    : RFCavity(name),
      scaleCore_m(1.0),
      scaleCoreError_m(0.0),
      phaseCore1_m(0.0),
      phaseCore2_m(0.0),
      phaseExit_m(0.0),
      startCoreField_m(0.0),
      startExitField_m(0.0),
      mappedStartExitField_m(0.0),
      periodLength_m(0.0),
      numCells_m(1),
      cellLength_m(0.0),
      mode_m(1) {
}

TravelingWave::~TravelingWave() {
}

void TravelingWave::accept(BeamlineVisitor& visitor) const {
    visitor.visitTravelingWave(*this);
}

/* ========================================================================== */
/* ============================== Apply Functions =========================== */
/**
 * @brief Applies the Traveling Wave RF Cavity field to all particles inside the RF cavity
 * Note: The field is applied in three regions: the core region, the exit region, and the 
 * transition region in between. The field in the transition region is calculated as the 
 * superposition of two fields from the core region, which are phase-shifted by 90 degrees 
 * and spatially shifted by one cell length. This approach allows for a smooth transition 
 * of the field from the core to the exit region, while maintaining the correct phase 
 * relationship between the fields.
 */
bool TravelingWave::apply(const std::shared_ptr<ParticleContainer_t>& pc)
{
    // RF parameters (copied to device)
    const double freq       = frequency_m;
    const double scaleEntry = scale_m + scaleError_m;
    const double scaleCore  = scaleCore_m + scaleCoreError_m;
    const double phaseError = phaseError_m;

    const double startField      = startField_m;
    const double startCoreField  = startCoreField_m;
    const double startExitField  = startExitField_m;
    const double mappedStartExit = mappedStartExitField_m;
    const double periodLength    = periodLength_m;
    const double cellLength      = cellLength_m;
    const double elementLength   = getElementLength();

    // Reference particle time -> to be consistent with OPAL
    const double t = RefPartBunch_m->getT() + 0.5 * RefPartBunch_m->getdT();
    const double omega_t = freq * t;

    const double cosEntry = std::cos(omega_t + phase_m + phaseError);
    const double sinEntry = std::sin(omega_t + phase_m + phaseError);

    const double cosCore1 = std::cos(omega_t + phaseCore1_m + phaseError);
    const double sinCore1 = std::sin(omega_t + phaseCore1_m + phaseError);

    const double cosCore2 = std::cos(omega_t + phaseCore2_m + phaseError);
    const double sinCore2 = std::sin(omega_t + phaseCore2_m + phaseError);

    const double cosExit = std::cos(omega_t + phaseExit_m + phaseError);
    const double sinExit = std::sin(omega_t + phaseExit_m + phaseError);

    auto* dynamicFieldmap = dynamic_cast<Astra1DDynamic*>(fieldmap_m);
    if (dynamicFieldmap == nullptr) {
        throw GeneralClassicException(
            "TravelingWave::apply",
            "TravelingWave particle application currently requires an Astra1DDynamic field map.");
    }

    dynamicFieldmap->applyTravelingWave(
        pc,
        scaleEntry * cosEntry,
        -scaleEntry * sinEntry,
        scaleCore * cosCore1,
        -scaleCore * sinCore1,
        scaleCore * cosCore2,
        -scaleCore * sinCore2,
        scaleEntry * cosExit,
        -scaleEntry * sinExit,
        startField,
        startCoreField,
        startExitField,
        mappedStartExit,
        periodLength,
        cellLength,
        elementLength);

    return false;
}

bool TravelingWave::apply(
    const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B) 
{
    std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    auto Rview = pc->R.getView();
    auto Pview = pc->P.getView();

    const Vector_t<double, 3> R = Rview(i);
    const Vector_t<double, 3> P = Pview(i);

    return apply(R, P, t, E, B);
}

bool TravelingWave::apply(
    const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& t,
    Vector_t<double, 3>& E, Vector_t<double, 3>& B) 
{
    const double omega_t = frequency_m * t;
    
    if (R(2) < -0.5 * periodLength_m || R(2) + 0.5 * periodLength_m >= getElementLength())
        return false;

    Vector_t<double, 3> tmpR({R(0), R(1), R(2) + 0.5 * periodLength_m});
    Vector_t<double, 3> tmpE({0.0, 0.0, 0.0}), tmpB({0.0, 0.0, 0.0});
    double tmpcos = 0.0, tmpsin = 0.0;

    if (tmpR(2) < startCoreField_m) {
        if (!fieldmap_m->isInside(tmpR))
            return getFlagDeleteOnTransverseExit();

        tmpcos = (scale_m + scaleError_m) * std::cos(omega_t + phase_m + phaseError_m);
        tmpsin = -(scale_m + scaleError_m) * std::sin(omega_t + phase_m + phaseError_m);

    } else if (tmpR(2) < startExitField_m) {
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2)        = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;

        if (!fieldmap_m->isInside(tmpR))
            return getFlagDeleteOnTransverseExit();

        tmpcos = (scaleCore_m + scaleCoreError_m)
                 * std::cos(omega_t + phaseCore1_m + phaseError_m);
        tmpsin = -(scaleCore_m + scaleCoreError_m)
                 * std::sin(omega_t + phaseCore1_m + phaseError_m);

        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        E += tmpcos * tmpE;
        B += tmpsin * tmpB;

        tmpE = 0.0;
        tmpB = 0.0;

        tmpR(2) = z + cellLength_m;
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos = (scaleCore_m + scaleCoreError_m)
                 * std::cos(omega_t + phaseCore2_m + phaseError_m);
        tmpsin = -(scaleCore_m + scaleCoreError_m)
                 * std::sin(omega_t + phaseCore2_m + phaseError_m);

    } else {
        tmpR(2) -= mappedStartExitField_m;
        if (!fieldmap_m->isInside(tmpR))
            return getFlagDeleteOnTransverseExit();

        tmpcos = (scale_m + scaleError_m) * std::cos(omega_t + phaseExit_m + phaseError_m);
        tmpsin = -(scale_m + scaleError_m) * std::sin(omega_t + phaseExit_m + phaseError_m);
    }

    fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
    E += tmpcos * tmpE;
    B += tmpsin * tmpB;

    return false;
}

bool TravelingWave::applyToReferenceParticle(
    const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& t,
    Vector_t<double, 3>& E, Vector_t<double, 3>& B) 
{
    const double omega_t = frequency_m * t;

    if (R(2) < -0.5 * periodLength_m || R(2) + 0.5 * periodLength_m >= getElementLength())
        return false;

    Vector_t<double, 3> tmpR({R(0), R(1), R(2) + 0.5 * periodLength_m});
    Vector_t<double, 3> tmpE({0.0, 0.0, 0.0}), tmpB({0.0, 0.0, 0.0});
    double tmpcos = 0.0, tmpsin = 0.0;

    if (tmpR(2) < startCoreField_m) {
        if (!fieldmap_m->isInside(tmpR))
            return true;

        tmpcos = scale_m * std::cos(omega_t + phase_m);
        tmpsin = -scale_m * std::sin(omega_t + phase_m);

    } else if (tmpR(2) < startExitField_m) {
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2)        = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;

        if (!fieldmap_m->isInside(tmpR))
            return true;

        tmpcos = scaleCore_m * std::cos(omega_t + phaseCore1_m);
        tmpsin = -scaleCore_m * std::sin(omega_t + phaseCore1_m);

        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        E += tmpcos * tmpE;
        B += tmpsin * tmpB;

        tmpE = 0.0;
        tmpB = 0.0;

        tmpR(2) = z + cellLength_m;
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos = scaleCore_m * std::cos(omega_t + phaseCore2_m);
        tmpsin = -scaleCore_m * std::sin(omega_t + phaseCore2_m);

    } else {
        tmpR(2) -= mappedStartExitField_m;
        if (!fieldmap_m->isInside(tmpR))
            return true;

        tmpcos = scale_m * std::cos(omega_t + phaseExit_m);
        tmpsin = -scale_m * std::sin(omega_t + phaseExit_m);
    }

    fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
    E += tmpcos * tmpE;
    B += tmpsin * tmpB;

    return false;
}

void TravelingWave::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    if (bunch == nullptr) {
        startField = -0.5 * periodLength_m;
        endField   = startExitField_m;
        return;
    }

    Inform msg("TravelingWave ", *gmsg);

    RefPartBunch_m = bunch;
    double zBegin = 0.0, zEnd = 0.0;
    RFCavity::initialise(bunch, zBegin, zEnd);
    if (std::abs(startField_m) > 0.0) {
        throw GeneralClassicException(
            "TravelingWave::initialise",
            "The field map of a traveling wave structure has to begin at 0.0");
    }

    periodLength_m = (zEnd - zBegin) / 2.0;
    cellLength_m   = periodLength_m * mode_m;
    startField_m   = -0.5 * periodLength_m;

    startCoreField_m       = periodLength_m / 2.0;
    startExitField_m       = startCoreField_m + (numCells_m - 1) * cellLength_m;
    mappedStartExitField_m = startExitField_m - 3.0 * periodLength_m / 2.0;

    startField = -periodLength_m / 2.0;
    endField   = startField + startExitField_m + periodLength_m / 2.0;
    setElementLength(endField - startField);

    scaleCore_m      = scale_m / std::sin(Physics::two_pi * mode_m);
    scaleCoreError_m = scaleError_m / std::sin(Physics::two_pi * mode_m);
    phaseCore1_m     = phase_m + Physics::pi * mode_m / 2.0;
    phaseCore2_m     = phase_m + Physics::pi * mode_m * 1.5;
    phaseExit_m =
        phase_m
        - Physics::two_pi * ((numCells_m - 1) * mode_m - std::floor((numCells_m - 1) * mode_m));
}

void TravelingWave::initialise(PartBunch_t* /*bunch*/, std::shared_ptr<AbstractTimeDependence> /*freq_atd*/,
                               std::shared_ptr<AbstractTimeDependence> /*ampl_atd*/,
                               std::shared_ptr<AbstractTimeDependence> /*phase_atd*/) {
  *gmsg << "not implemented:: file: " << __FILE__ << " line: " << __LINE__ << " function: " << __func__ << endl;
}


void TravelingWave::finalise() {
}

bool TravelingWave::bends() const {
    return false;
}

void TravelingWave::goOnline(const double&) {
    Fieldmap::readMap(filename_m);
    online_m = true;
}

void TravelingWave::goOffline() {
    Fieldmap::freeMap(filename_m);
}

void TravelingWave::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = -0.5 * periodLength_m;
    zEnd   = zBegin + getElementLength();
}

void TravelingWave::getElementDimensions(double& begin, double& end) const {
    begin = -0.5 * periodLength_m;
    end   = begin + getElementLength();
}

ElementType TravelingWave::getType() const {
    return ElementType::TRAVELINGWAVE;
}

double TravelingWave::getAutoPhaseEstimate(
    const double& E0, const double& t0, const double& q, const double& mass) {
    std::vector<double> t, E, t2, E2;
    std::vector<std::pair<double, double> > F;
    double phi = 0.0, tmp_phi, dphi = 0.5 * Units::deg2rad;
    double phaseC1 = phaseCore1_m - phase_m;
    double phaseC2 = phaseCore2_m - phase_m;
    double phaseE  = phaseExit_m - phase_m;

    fieldmap_m->getOnaxisEz(F);
    if (F.size() == 0)
        return 0.0;

    int N1    = static_cast<int>(std::floor(F.size() / 4.)) + 1;
    int N2    = F.size() - 2 * N1 + 1;
    int N3    = 2 * N1 + static_cast<int>(std::floor((numCells_m - 1) * N2 * mode_m)) - 1;
    int N4    = static_cast<int>(std::round(N2 * mode_m));
    double Dz = F[N1 + N2].first - F[N1].first;

    t.resize(N3, t0);
    t2.resize(N3, t0);
    E.resize(N3, E0);
    E2.resize(N3, E0);
    for (int i = 1; i < N1; ++i) {
        E[i]  = E0 + (F[i].first + F[i - 1].first) / 2. * scale_m / mass;
        E2[i] = E[i];
    }
    for (int i = N1; i < N3 - N1 + 1; ++i) {
        int I    = (i - N1) % N2 + N1;
        double Z = (F[I].first + F[I - 1].first) / 2. + std::floor((i - N1) / N2) * Dz;
        E[i]     = E0 + Z * scaleCore_m / mass;
        E2[i]    = E[i];
    }
    for (int i = N3 - N1 + 1; i < N3; ++i) {
        int I    = i - N3 - 1 + 2 * N1 + N2;
        double Z = (F[I].first + F[I - 1].first) / 2. + ((numCells_m - 1) * mode_m - 1) * Dz;
        E[i]     = E0 + Z * scale_m / mass;
        E2[i]    = E[i];
    }

    for (int iter = 0; iter < 10; ++iter) {
        double A = 0.0;
        double B = 0.0;
        for (int i = 1; i < N1; ++i) {
            t[i]  = t[i - 1] + getdT(i, i, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, i, E2, F, mass);
            A += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdA(i, i, t, 0.0, F);
            B += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdB(i, i, t, 0.0, F);
        }
        for (int i = N1; i < N3 - N1 + 1; ++i) {
            int I = (i - N1) % N2 + N1;
            int J = (i - N1 + N4) % N2 + N1;
            t[i]  = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            A += scaleCore_m * (1. + frequency_m * (t2[i] - t[i]) / dphi)
                 * (getdA(i, I, t, phaseC1, F) + getdA(i, J, t, phaseC2, F));
            B += scaleCore_m * (1. + frequency_m * (t2[i] - t[i]) / dphi)
                 * (getdB(i, I, t, phaseC1, F) + getdB(i, J, t, phaseC2, F));
        }
        for (int i = N3 - N1 + 1; i < N3; ++i) {
            int I = i - N3 - 1 + 2 * N1 + N2;
            t[i]  = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            A += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdA(i, I, t, phaseE, F);
            B += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdB(i, I, t, phaseE, F);
        }

        if (std::abs(B) > 0.0000001) {
            tmp_phi = std::atan(A / B);
        } else {
            tmp_phi = Physics::pi / 2;
        }
        if (q * (A * std::sin(tmp_phi) + B * std::cos(tmp_phi)) < 0) {
            tmp_phi += Physics::pi;
        }

        if (std::abs(phi - tmp_phi) < frequency_m * (t[N3 - 1] - t[0]) / N3) {
            for (int i = 1; i < N1; ++i) {
                E[i] = E[i - 1] + q * scale_m * getdE(i, i, t, phi, F);
            }
            for (int i = N1; i < N3 - N1 + 1; ++i) {
                int I = (i - N1) % N2 + N1;
                int J = (i - N1 + N4) % N2 + N1;
                E[i] =
                    E[i - 1]
                    + q * scaleCore_m
                          * (getdE(i, I, t, phi + phaseC1, F) + getdE(i, J, t, phi + phaseC2, F));
            }
            for (int i = N3 - N1 + 1; i < N3; ++i) {
                int I = i - N3 - 1 + 2 * N1 + N2;
                E[i]  = E[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE, F);
            }

            const int prevPrecision = ippl::Info->precision(8);
            *ippl::Info << level2 << "estimated phase= " << tmp_phi
                        << " rad = " << tmp_phi * Units::rad2deg << " deg,\n"
                        << "Ekin= " << E[N3 - 1] << " MeV" << std::setprecision(prevPrecision)
                        << "\n"
                        << endl;
            return tmp_phi;
        }
        phi = tmp_phi - std::round(tmp_phi / Physics::two_pi) * Physics::two_pi;

        for (int i = 1; i < N1; ++i) {
            E[i]  = E[i - 1] + q * scale_m * getdE(i, i, t, phi, F);
            E2[i] = E2[i - 1]
                    + q * scale_m * getdE(i, i, t, phi + dphi, F);  // should I use here t or t2?
            t[i]  = t[i - 1] + getdT(i, i, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, i, E2, F, mass);
            E[i]  = E[i - 1] + q * scale_m * getdE(i, i, t, phi, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, i, t2, phi + dphi, F);
        }
        for (int i = N1; i < N3 - N1 + 1; ++i) {
            int I = (i - N1) % N2 + N1;
            int J = (i - N1 + N4) % N2 + N1;
            E[i]  = E[i - 1]
                   + q * scaleCore_m
                         * (getdE(i, I, t, phi + phaseC1, F) + getdE(i, J, t, phi + phaseC2, F));
            E2[i] = E2[i - 1]
                    + q * scaleCore_m
                          * (getdE(i, I, t, phi + phaseC1 + dphi, F)
                             + getdE(i, J, t, phi + phaseC2 + dphi, F));  // concerning t: see above
            t[i]  = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            E[i]  = E[i - 1]
                   + q * scaleCore_m
                         * (getdE(i, I, t, phi + phaseC1, F) + getdE(i, J, t, phi + phaseC2, F));
            E2[i] = E2[i - 1]
                    + q * scaleCore_m
                          * (getdE(i, I, t2, phi + phaseC1 + dphi, F)
                             + getdE(i, J, t2, phi + phaseC2 + dphi, F));
        }
        for (int i = N3 - N1 + 1; i < N3; ++i) {
            int I = i - N3 - 1 + 2 * N1 + N2;
            E[i]  = E[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE, F);
            E2[i] =
                E2[i - 1]
                + q * scale_m * getdE(i, I, t, phi + phaseE + dphi, F);  // concerning t: see above
            t[i]  = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            E[i]  = E[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, I, t2, phi + phaseE + dphi, F);
        }
        //             msg << ", Ekin= " << E[N3-1] << " MeV" << endl;
    }

    const int prevPrecision = ippl::Info->precision(8);
    *ippl::Info << level2 << "estimated phase= " << tmp_phi << " rad = " << tmp_phi * Units::rad2deg
                << " deg,\n"
                << "Ekin= " << E[N3 - 1] << " MeV" << std::setprecision(prevPrecision) << "\n"
                << endl;

    return phi;
}

bool TravelingWave::isInside(const Vector_t<double, 3>& r) const {
    return (isInsideTransverse(r) && r(2) >= -0.5 * periodLength_m && r(2) < startExitField_m);
}
