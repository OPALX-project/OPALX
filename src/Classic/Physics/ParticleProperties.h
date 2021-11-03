//
// Class ParticleProperties
//   Base class for representing particle properties
//
// Copyright (c) 2021, Pedro Calvo, CIEMAT, Spain
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
#ifndef PARTICLEPROPERTIES_H
#define PARTICLEPROPERTIES_H

#include <boost/bimap.hpp>

#include <map>
#include <string>

enum class ParticleType:short {
    UNNAMED = -1,
    ELECTRON,
    POSITRON,
    MUON,
    PROTON,
    ANTIPROTON,
    DEUTERON,
    HMINUS,
    HYDROGEN,
    H2P,
    H3P,
    ALPHA,
    CARBON,
    XENON,
    URANIUM
};

enum class ParticleOrigin:short {
    REGULAR,
    SECONDARY,
    STRIPPED
};


class ParticleProperties {
public:
    static ParticleType getParticleType (const std::string& str);

    static std::string getParticleTypeString(ParticleType type);

    static double getParticleMass(ParticleType type);

    static double getParticleCharge(ParticleType type);
    static double getParticleChargeInCoulomb(ParticleType type);

private:
    static const boost::bimap<ParticleType, std::string> bmParticleType_m;

    static const std::map<ParticleType, double> particleMass_m;
    static const std::map<ParticleType, double> particleCharge_m;
};

#endif /* PARTICLEPROPERTIES_H */
