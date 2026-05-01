#include <gtest/gtest.h>

#include "Attributes/Attributes.h"
#include "Physics/ParticleProperties.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"

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

}  // namespace
