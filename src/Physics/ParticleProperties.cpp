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
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"

#include "Utilities/BiMap.h"

ParticleType ParticleProperties::getParticleType(const std::string& str) {
    auto it = bmParticleType_s.right.find(str);
    if (it != bmParticleType_s.right.end()) {
        return it->second;
    } else {
        return ParticleType::UNNAMED;
    }
}

std::string ParticleProperties::getParticleTypeString(const ParticleType& type) {
    return bmParticleType_s.left.at(type);
}

double ParticleProperties::getParticleMass(const ParticleType& type) {
    return particleMass_m.at(type);
}

double ParticleProperties::getParticleCharge(const ParticleType& type) {
    return particleCharge_m.at(type);
}

double ParticleProperties::getParticleChargeInCoulomb(const ParticleType& type) {
    return getParticleCharge(type) * Physics::q_e;
}

double ParticleProperties::getParticleLifetime(const ParticleType& type) {
    auto it = particleLifetime_m.find(type);
    if (it != particleLifetime_m.end()) {
        return it->second;
    }
    throw std::runtime_error(
            "ParticleProperties::getParticleLifetime: no lifetime for particle \""
            + getParticleTypeString(type) + "\" (stable or unsupported).");
}

const BiMap<ParticleType, std::string> ParticleProperties::bmParticleType_s = []() {
    BiMap<ParticleType, std::string> bimap;
    bimap.insert(ParticleType::UNNAMED, "UNNAMED");
    bimap.insert(ParticleType::PHOTON, "PHOTON");
    bimap.insert(ParticleType::ELECTRON, "ELECTRON");
    bimap.insert(ParticleType::POSITRON, "POSITRON");
    bimap.insert(ParticleType::MUON, "MUON");
    bimap.insert(ParticleType::PION, "PION");
    bimap.insert(ParticleType::PROTON, "PROTON");
    bimap.insert(ParticleType::ANTIPROTON, "ANTIPROTON");
    bimap.insert(ParticleType::DEUTERON, "DEUTERON");
    bimap.insert(ParticleType::HMINUS, "HMINUS");
    bimap.insert(ParticleType::HYDROGEN, "HYDROGEN");
    bimap.insert(ParticleType::H2P, "H2P");
    bimap.insert(ParticleType::H3P, "H3P");
    bimap.insert(ParticleType::ALPHA, "ALPHA");
    bimap.insert(ParticleType::CARBON, "CARBON");
    bimap.insert(ParticleType::XENON, "XENON");
    bimap.insert(ParticleType::URANIUM, "URANIUM");
    return bimap;
}();

const std::map<ParticleType, double> ParticleProperties::particleMass_m = {
        {ParticleType::PHOTON, 0.0},
        {ParticleType::ELECTRON, Physics::m_e},
        {ParticleType::POSITRON, Physics::m_e},
        {ParticleType::MUON, Physics::m_mu},
        {ParticleType::PION, Physics::m_pi},
        {ParticleType::PROTON, Physics::m_p},
        {ParticleType::ANTIPROTON, Physics::m_p},
        {ParticleType::DEUTERON, Physics::m_d},
        {ParticleType::HMINUS, Physics::m_hm},
        {ParticleType::HYDROGEN, Physics::m_h},
        {ParticleType::H2P, Physics::m_h2p},
        {ParticleType::H3P, Physics::m_h3p},
        {ParticleType::ALPHA, Physics::m_alpha},
        {ParticleType::CARBON, Physics::m_c},
        {ParticleType::XENON, Physics::m_xe},
        {ParticleType::URANIUM, Physics::m_u}};

const std::map<ParticleType, double> ParticleProperties::particleCharge_m = {
        {ParticleType::PHOTON, 0.0},      {ParticleType::ELECTRON, -1.0},
        {ParticleType::POSITRON, 1.0},    {ParticleType::MUON, -1.0},
        {ParticleType::PION, 1.0},        {ParticleType::PROTON, 1.0},
        {ParticleType::ANTIPROTON, -1.0}, {ParticleType::DEUTERON, 1.0},
        {ParticleType::HMINUS, -1.0},     {ParticleType::HYDROGEN, 0.0},
        {ParticleType::H2P, 1.0},         {ParticleType::H3P, 1.0},
        {ParticleType::ALPHA, 2.0},       {ParticleType::CARBON, 6.0},
        {ParticleType::XENON, 54.0},      {ParticleType::URANIUM, 92.0}};

const std::map<ParticleType, double> ParticleProperties::particleLifetime_m = {
        {ParticleType::MUON, Physics::tau_mu}, {ParticleType::PION, Physics::tau_pi}};
