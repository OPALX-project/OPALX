#include <gtest/gtest.h>

#include "Attributes/Attributes.h"
#include "Physics/ParticleProperties.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <vector>
#include "../KeywordObjects.h"

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
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeValidPhotonBeam(beam);

        EXPECT_NO_THROW(beam.execute());
        EXPECT_TRUE(beam.isPhoton());
        EXPECT_DOUBLE_EQ(beam.getMass(), 0.0);
        EXPECT_DOUBLE_EQ(beam.getCharge(), 0.0);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRequiresEnergy) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        setParticle(beam, "PHOTON");
        setNALLOC(beam, 1.0);

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRejectsGamma) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeValidPhotonBeam(beam);
        setGamma(beam, 10.0);

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRejectsMomentum) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeValidPhotonBeam(beam);
        setPC(beam, 1.0);

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRejectsExplicitMassAndCharge) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeValidPhotonBeam(beam);
        setMass(beam, 1.0);

        EXPECT_THROW(beam.execute(), OpalException);

        auto beam2Object = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam2      = *beam2Object;
        makeValidPhotonBeam(beam2);
        setCharge(beam2, 1.0);

        EXPECT_THROW(beam2.execute(), OpalException);
    }

    TEST_F(BeamPhotonTest, PhotonBeamRejectsSources) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeValidPhotonBeam(beam);
        setSources(beam, "srcs");

        EXPECT_THROW(beam.execute(), OpalException);
    }

    // A DAUGHTERBEAM has no GAMMA/ENERGY/PC and no SOURCES. Beam::execute must accept it;
    // TrackRun requires SOURCES on every beam that is not a DAUGHTERBEAM.
    TEST_F(BeamPhotonTest, BeamWithoutEnergyOrSourcesParses) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        setParticle(beam, "ELECTRON");
        setNALLOC(beam, 1.0);

        EXPECT_NO_THROW(beam.execute());
        EXPECT_FALSE(beam.hasExplicitEnergy());
        EXPECT_FALSE(beam.hasEmissionSources());

        setSources(beam, "src");
        EXPECT_TRUE(beam.hasEmissionSources());
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
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
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
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeMuonBeam(beam);
        setPolarization(beam, {0.0, 0.0, 1.0});

        EXPECT_NO_THROW(beam.execute());
        EXPECT_TRUE(beam.hasPolarization());
    }

    TEST_F(BeamPolarizationTest, MuonZeroPolarizationEnablesSpinTracking) {
        // An explicit {0,0,0} still counts as "set" and enables spin tracking.
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeMuonBeam(beam);
        setPolarization(beam, {0.0, 0.0, 0.0});

        EXPECT_NO_THROW(beam.execute());
        EXPECT_TRUE(beam.hasPolarization());
    }

    TEST_F(BeamPolarizationTest, MuonWithoutPolarizationHasSpinOff) {
        // Leaving POLARIZATION unset is valid and disables spin tracking.
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeMuonBeam(beam);

        EXPECT_NO_THROW(beam.execute());
        EXPECT_FALSE(beam.hasPolarization());
    }

    TEST_F(BeamPolarizationTest, PolarizationMagnitudeGreaterThanOneIsRejected) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeMuonBeam(beam);
        setPolarization(beam, {0.8, 0.8, 0.8});  // |P| ~ 1.39

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPolarizationTest, PolarizationWrongLengthIsRejected) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeMuonBeam(beam);
        setPolarization(beam, {0.0, 1.0});  // length 2, not 3

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPolarizationTest, NonMuonSpeciesRejectsNonzeroPolarization) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        setParticle(beam, "PION");
        setEnergy(beam, 1.0);
        setNALLOC(beam, 1.0);
        setSources(beam, "src");
        setPolarization(beam, {0.0, 0.0, 1.0});

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPolarizationTest, NonMuonSpeciesRejectsZeroPolarization) {
        // POLARIZATION is not valid on non-muon species, even at magnitude 0.
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        setParticle(beam, "PION");
        setEnergy(beam, 1.0);
        setNALLOC(beam, 1.0);
        setSources(beam, "src");
        setPolarization(beam, {0.0, 0.0, 0.0});

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamPolarizationTest, NonMuonSpeciesAcceptedWhenPolarizationUnset) {
        // Leaving POLARIZATION unspecified is fine on any species.
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        setParticle(beam, "PION");
        setEnergy(beam, 1.0);
        setNALLOC(beam, 1.0);
        setSources(beam, "src");

        EXPECT_NO_THROW(beam.execute());
    }

    TEST_F(BeamPhotonTest, ZStopAndUntilMaxStepsAreUnsetByDefault) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        EXPECT_FALSE(beam.getZStop().has_value());
        EXPECT_FALSE(beam.getUntilMaxSteps());

        Attributes::setReal(*beam.findAttribute("ZSTOP"), 12.5);
        Attributes::setBool(*beam.findAttribute("UNTILMAXSTEPS"), true);
        ASSERT_TRUE(beam.getZStop().has_value());
        EXPECT_DOUBLE_EQ(*beam.getZStop(), 12.5);
        EXPECT_TRUE(beam.getUntilMaxSteps());
    }

    class BeamReferenceDirectionTest : public BeamPhotonTest {
    protected:
        void setReferenceDirection(Beam& beam, const std::vector<double>& dir) {
            Attributes::setRealArray(*beam.findAttribute("REFDIRECTION"), dir);
        }
        void makeMuonBeam(Beam& beam) {
            setParticle(beam, "MUON");
            setPC(beam, 0.028);
            setNALLOC(beam, 1.0);
            setSources(beam, "src");
        }
    };

    TEST_F(BeamReferenceDirectionTest, DefaultIsZ) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeMuonBeam(beam);

        EXPECT_NO_THROW(beam.execute());
        EXPECT_FALSE(beam.hasReferenceDirection());
        EXPECT_EQ(beam.getReferenceDirection(), (std::vector<double>{0.0, 0.0, 1.0}));
    }

    TEST_F(BeamReferenceDirectionTest, IsNormalised) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeMuonBeam(beam);
        setReferenceDirection(beam, {3.0, 0.0, 4.0});

        EXPECT_NO_THROW(beam.execute());
        const std::vector<double> dir = beam.getReferenceDirection();
        ASSERT_EQ(dir.size(), 3u);
        EXPECT_DOUBLE_EQ(dir[0], 0.6);
        EXPECT_DOUBLE_EQ(dir[1], 0.0);
        EXPECT_DOUBLE_EQ(dir[2], 0.8);
    }

    TEST_F(BeamReferenceDirectionTest, WrongLengthIsRejected) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeMuonBeam(beam);
        setReferenceDirection(beam, {0.0, 1.0});

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamReferenceDirectionTest, ZeroIsRejected) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        makeMuonBeam(beam);
        setReferenceDirection(beam, {0.0, 0.0, 0.0});

        EXPECT_THROW(beam.execute(), OpalException);
    }

    TEST_F(BeamReferenceDirectionTest, NeedsEnergy) {
        auto beamObject = opalx::test::makeKeywordObject<Beam>("BEAM", "BEAM");
        Beam& beam      = *beamObject;
        setParticle(beam, "MUON");
        setNALLOC(beam, 1.0);
        setReferenceDirection(beam, {0.0, 0.0, 1.0});

        EXPECT_THROW(beam.execute(), OpalException);
    }

}  // namespace
