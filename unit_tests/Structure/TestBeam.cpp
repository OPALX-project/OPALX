#include <gtest/gtest.h>

#include "Attributes/Attributes.h"
#include "Physics/ParticleProperties.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <vector>

namespace {

    class BeamPhotonTest : public ::testing::Test {
    protected:
        void setParticle(Beam& beam, const std::string& name) {
            Attributes::setPredefinedString(*beam.findAttribute("PARTICLE"), name);
        }

        void setEnergy(Beam& beam, double value) {
            Attributes::setReal(*beam.findAttribute("ENERGY"), value);
        }

        void setGamma(Beam& beam, double value) {
            Attributes::setReal(*beam.findAttribute("GAMMA"), value);
        }

        void setPC(Beam& beam, double value) {
            Attributes::setReal(*beam.findAttribute("PC"), value);
        }

        void setMass(Beam& beam, double value) {
            Attributes::setReal(*beam.findAttribute("MASS"), value);
        }

        void setCharge(Beam& beam, double value) {
            Attributes::setReal(*beam.findAttribute("CHARGE"), value);
        }

        void setSources(Beam& beam, const std::string& value) {
            Attributes::setString(*beam.findAttribute("SOURCES"), value);
        }

        void setNALLOC(Beam& beam, double value) {
            Attributes::setReal(*beam.findAttribute("NALLOC"), value);
        }

        void makeValidPhotonBeam(Beam& beam) {
            setParticle(beam, "PHOTON");
            setEnergy(beam, 5.0);
            setNALLOC(beam, 1.0);
        }
    };

    TEST_F(BeamPhotonTest, PhotonParticlePropertiesAreRegistered) {
        EXPECT_EQ(ParticleProperties::getParticleType("PHOTON"), ParticleType::PHOTON);
        EXPECT_EQ(ParticleProperties::getParticleTypeString(ParticleType::PHOTON), "PHOTON");
        EXPECT_DOUBLE_EQ(ParticleProperties::getParticleMass(ParticleType::PHOTON), 0.0);
        EXPECT_DOUBLE_EQ(ParticleProperties::getParticleCharge(ParticleType::PHOTON), 0.0);
    }

    TEST_F(BeamPhotonTest, PhotonBeamAcceptsEnergyOnlyDefinition) {
        Beam beam;
        makeValidPhotonBeam(beam);

        EXPECT_NO_THROW(beam.execute());
        EXPECT_TRUE(beam.isPhoton());
        EXPECT_DOUBLE_EQ(beam.getMass(), 0.0);
        EXPECT_DOUBLE_EQ(beam.getCharge(), 0.0);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRequiresEnergy) {
        Beam beam;
        setParticle(beam, "PHOTON");
        setNALLOC(beam, 1.0);

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRejectsGamma) {
        Beam beam;
        makeValidPhotonBeam(beam);
        setGamma(beam, 10.0);

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRejectsMomentum) {
        Beam beam;
        makeValidPhotonBeam(beam);
        setPC(beam, 1.0);

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRejectsExplicitMassAndCharge) {
        Beam beam;
        makeValidPhotonBeam(beam);
        setMass(beam, 1.0);

        EXPECT_THROW(beam.execute(), OpalException);

        Beam beam2;
        makeValidPhotonBeam(beam2);
        setCharge(beam2, 1.0);

        EXPECT_THROW(beam2.execute(), OpalException);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRejectsSources) {
        Beam beam;
        makeValidPhotonBeam(beam);
        setSources(beam, "srcs");

        EXPECT_THROW(beam.execute(), OpalException);
    }

    class BeamPolarizationTest : public BeamPhotonTest {
    protected:
        void setPolarization(Beam& beam, const std::vector<double>& pol) {
            Attributes::setRealArray(*beam.findAttribute("POLARIZATION"), pol);
        }
        void makeMuonBeam(Beam& beam) {
            setParticle(beam, "MUON");
            setEnergy(beam, 0.5);  // GeV, well above muon mass
            setNALLOC(beam, 1.0);
            setSources(beam, "src");
        }
    };

    TEST_F(BeamPolarizationTest, MuonBeamDefaultPolarizationIsZero) {
        Beam beam;
        makeMuonBeam(beam);  // POLARIZATION unset: getPolarization() defaults to {0,0,0}.

        EXPECT_NO_THROW(beam.execute());
        const std::vector<double> pol = beam.getPolarization();
        ASSERT_EQ(pol.size(), 3u);
        EXPECT_DOUBLE_EQ(pol[0], 0.0);
        EXPECT_DOUBLE_EQ(pol[1], 0.0);
        EXPECT_DOUBLE_EQ(pol[2], 0.0);
    }

    TEST_F(BeamPolarizationTest, MuonNonzeroPolarizationAccepted) {
        // Setting POLARIZATION on a muon is accepted and enables spin tracking.
        Beam beam;
        makeMuonBeam(beam);
        setPolarization(beam, {0.0, 0.0, 1.0});

        EXPECT_NO_THROW(beam.execute());
        EXPECT_TRUE(beam.hasPolarization());
    }

    TEST_F(BeamPolarizationTest, MuonZeroPolarizationEnablesSpinTracking) {
        // An explicit {0,0,0} still counts as "set" and enables spin tracking.
        Beam beam;
        makeMuonBeam(beam);
        setPolarization(beam, {0.0, 0.0, 0.0});

        EXPECT_NO_THROW(beam.execute());
        EXPECT_TRUE(beam.hasPolarization());
    }

    TEST_F(BeamPolarizationTest, MuonWithoutPolarizationHasSpinOff) {
        // Leaving POLARIZATION unset is valid and disables spin tracking.
        Beam beam;
        makeMuonBeam(beam);

        EXPECT_NO_THROW(beam.execute());
        EXPECT_FALSE(beam.hasPolarization());
    }

    TEST_F(BeamPolarizationTest, PolarizationMagnitudeGreaterThanOneIsRejected) {
        Beam beam;
        makeMuonBeam(beam);
        setPolarization(beam, {0.8, 0.8, 0.8});  // |P| ~ 1.39

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPolarizationTest, PolarizationWrongLengthIsRejected) {
        Beam beam;
        makeMuonBeam(beam);
        setPolarization(beam, {0.0, 1.0});  // length 2, not 3

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPolarizationTest, NonMuonSpeciesRejectsNonzeroPolarization) {
        Beam beam;
        setParticle(beam, "PION");
        setEnergy(beam, 1.0);
        setNALLOC(beam, 1.0);
        setSources(beam, "src");
        setPolarization(beam, {0.0, 0.0, 1.0});

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPolarizationTest, NonMuonSpeciesRejectsZeroPolarization) {
        // POLARIZATION is not valid on non-muon species, even at magnitude 0.
        Beam beam;
        setParticle(beam, "PION");
        setEnergy(beam, 1.0);
        setNALLOC(beam, 1.0);
        setSources(beam, "src");
        setPolarization(beam, {0.0, 0.0, 0.0});

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPolarizationTest, NonMuonSpeciesAcceptedWhenPolarizationUnset) {
        // Leaving POLARIZATION unspecified is fine on any species.
        Beam beam;
        setParticle(beam, "PION");
        setEnergy(beam, 1.0);
        setNALLOC(beam, 1.0);
        setSources(beam, "src");

        EXPECT_NO_THROW(beam.execute());
    }

}  // namespace
