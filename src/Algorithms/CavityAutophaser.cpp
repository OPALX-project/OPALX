#include "Algorithms/CavityAutophaser.h"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "Utilities/Options.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"
#include "AbstractObjects/OpalData.h"

extern Inform *gmsg;

CavityAutophaser::CavityAutophaser(const PartData &ref,
                                   std::shared_ptr<Component> cavity):
    itsReference_m(ref),
    itsCavity_m(cavity)
{
    double zbegin = 0.0, zend = 0.0;
    PartBunch *fakeBunch = NULL;
    cavity->initialise(fakeBunch, zbegin, zend);
    initialR_m = Vector_t(0, 0, zbegin);
}

CavityAutophaser::~CavityAutophaser() {

}

double CavityAutophaser::getPhaseAtMaxEnergy(const Vector_t &R,
                                             const Vector_t &P,
                                             double t,
                                             double dt) {
    initialP_m = Vector_t(0, 0, euclidian_norm(P));
    double tErr  = (initialR_m(2) - R(2)) * sqrt(dot(P,P) + 1.0) / (P(2) * Physics::c);

    double initialEnergy = Util::getEnergy(P, itsReference_m.getM()) * 1e-6;
    double originalPhase = 0.0, newPhase = 0.0, AstraPhase = 0.0;
    double initialPhase = guessCavityPhase(t + tErr);
    double optimizedPhase = 0.0;
    double amplitude = 0.0;
    double designEnergy = 0.0;
    double finalEnergy = 0.0;
    const double length = itsCavity_m->getElementLength();

    if(!(itsCavity_m->getType() == ElementBase::TRAVELINGWAVE ||
         itsCavity_m->getType() == ElementBase::RFCAVITY)) {
        throw OpalException("CavityAutophaser::getPhaseAtMaxEnergy()",
                            "given element is not a cavity");
    }

    RFCavity *element = static_cast<RFCavity *>(itsCavity_m.get());
    amplitude = element->getAmplitudem();
    designEnergy = element->getDesignEnergy();
    originalPhase = element->getPhasem();

    if (amplitude == 0.0) {
        if (designEnergy <= 0.0 || length <= 0.0) {
            throw OpalException("CavityAutophaser::getPhaseAtMaxEnergy()",
                                "neither amplitude or design energy given to cavity " + element->getName());
        }

        amplitude = 2 * (designEnergy - initialEnergy) / (std::abs(itsReference_m.getQ()) * length);

        element->setAmplitudem(amplitude);
    }

    if (designEnergy > 0.0) {
        int count = 0;
        while (count < 1000) {
            initialPhase = guessCavityPhase(t + tErr);
            auto status = optimizeCavityPhase(initialPhase, t + tErr, dt);

            optimizedPhase = status.first;
            finalEnergy = status.second;

            if (std::abs(designEnergy - finalEnergy) < 1e-7) break;

            amplitude *= std::abs(designEnergy / finalEnergy);
            element->setAmplitudem(amplitude);
            initialPhase = optimizedPhase;

            ++ count;
        }
    }
    auto status = optimizeCavityPhase(initialPhase, t + tErr, dt);

    optimizedPhase = status.first;
    finalEnergy = status.second;


    AstraPhase = std::fmod(optimizedPhase + Physics::pi / 2, Physics::two_pi);
    newPhase = std::fmod(originalPhase + optimizedPhase + Physics::two_pi, Physics::two_pi);
    element->setPhasem(newPhase);
    element->setAutophaseVeto();


    INFOMSG(itsCavity_m->getName() << "_phi = "  << optimizedPhase * Physics::rad2deg <<  " [deg], "
            << "corresp. in Astra = " << AstraPhase * Physics::rad2deg << " [deg],\n"
            << "E = " << finalEnergy << " [MeV], " << "phi_nom = " << originalPhase * Physics::rad2deg << " [deg]\n"
            << "Ez_0 = " << amplitude << " [MV/m]" << "\n"
            << "time = " << (t + tErr) * 1e9 << " [ns], dt = " << dt * 1e12 << " [ps]" << endl);

    OpalData::getInstance()->setMaxPhase(itsCavity_m->getName(), optimizedPhase);


    return optimizedPhase;
}

double CavityAutophaser::guessCavityPhase(double t) {
    const Vector_t &refP = initialP_m;
    double Phimax = 0.0;
    bool apVeto;
    RFCavity *element = static_cast<RFCavity *>(itsCavity_m.get());
    double orig_phi = element->getPhasem();
    apVeto = element->getAutophaseVeto();
    if (apVeto) {
        INFOMSG(level1 << " ----> APVETO -----> "
                << element->getName() << endl);
        return orig_phi;
    }

    Phimax = element->getAutoPhaseEstimate(getEnergyMeV(refP),
                                           t,
                                           itsReference_m.getQ(),
                                           itsReference_m.getM() * 1e-6);

    return std::fmod(Phimax + Physics::two_pi, Physics::two_pi);
}

std::pair<double, double> CavityAutophaser::optimizeCavityPhase(double initialPhase,
                                                                double t,
                                                                double dt) {

    double originalPhase = 0.0;
    RFCavity *element = static_cast<RFCavity *>(itsCavity_m.get());
    originalPhase = element->getPhasem();

    if (element->getAutophaseVeto()) {
        std::pair<double, double> status(originalPhase, Util::getEnergy(initialP_m, itsReference_m.getM() * 1e-6));
        return status;
    }

    double Phimax = initialPhase;
    double phi = initialPhase;
    double dphi = Physics::pi / 360.0;
    const int numRefinements = Options::autoPhase;

    int j = -1;
    double E = track(initialR_m, initialP_m, t, dt, phi);
    double Emax = E;

    do {
        j ++;
        Emax = E;
        initialPhase = phi;
        phi -= dphi;
        E = track(initialR_m, initialP_m, t, dt, phi);
    } while(E > Emax);

    if(j == 0) {
        phi = initialPhase;
        E = Emax;
        j = -1;
        do {
            j ++;
            Emax = E;
            initialPhase = phi;
            phi += dphi;
            E = track(initialR_m, initialP_m, t, dt, phi);
        } while(E > Emax);
    }

    for(int rl = 0; rl < numRefinements; ++ rl) {
        dphi /= 2.;
        phi = initialPhase - dphi;
        E = track(initialR_m, initialP_m, t, dt, phi);
        if(E > Emax) {
            initialPhase = phi;
            Emax = E;
        } else {
            phi = initialPhase + dphi;
            E = track(initialR_m, initialP_m, t, dt, phi);
            if(E > Emax) {
                initialPhase = phi;
                Emax = E;
            }
        }
    }
    Phimax = std::fmod(initialPhase + Physics::two_pi, Physics::two_pi);

    E = track(initialR_m, initialP_m, t, dt, Phimax + originalPhase);
    std::pair<double, double> status(Phimax, E);
    return status;
}

double CavityAutophaser::track(Vector_t R,
                               Vector_t P,
                               double t,
                               const double dt,
                               const double phase) const {
    const Vector_t &refP = initialP_m;
    double initialPhase = 0.0;
    double finalMomentum = 0.0;

    RFCavity *rfc = static_cast<RFCavity *>(itsCavity_m.get());
    initialPhase = rfc->getPhasem();
    rfc->setPhasem(phase);

    std::pair<double, double> pe = rfc->trackOnAxisParticle(refP(2),
                                                            t,
                                                            dt,
                                                            itsReference_m.getQ(),
                                                            itsReference_m.getM() * 1e-6);
    finalMomentum = pe.first;
    rfc->setPhasem(initialPhase);

    double finalGamma = sqrt(1.0 + finalMomentum * finalMomentum);
    double finalKineticEnergy = (finalGamma - 1.0) * itsReference_m.getM() * 1e-6;

    return finalKineticEnergy;
}