//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPAL_LOCALPROCESS_H
#define OPAL_LOCALPROCESS_H

#include "Algorithms/CoordinateSystemTrafo.h"

#include <cstddef>
#include <string>

template <typename T, unsigned Dim>
class ParticleContainer;
class Component;

/**
 * @brief Information likely required by local processes 
 */
struct LocalProcessContext {
    double dt;                         
    long long globalTrackStep;         
    std::size_t containerIdx;          
    CoordinateSystemTrafo refToLocal;  
};

/**
 * @brief Base class for physics processes attached to individual beamline elements.
 *
 * Local processes complement global processes (see GlobalProcess): instead of acting on
 * all particles everywhere, they are attached to specific beamline elements in the input
 * file (e.g. synchrotron radiation on a bend, particle-matter interaction on a collimator)
 * and act only while the bunch overlaps that element.
 * 
 * @note Structure: 
 *   - apply() is called once per iteration for every process on every element and on every
 *   particle container. 
 *   - Selecting the particles that are actually inside the
 *   element is the process's responsibility.
 *   - The return value is the number of particles marked invalid. 
 *   
 */
class LocalProcess {
public:
    explicit LocalProcess(const std::string& name);
    virtual ~LocalProcess();

    /// Definition name from the input file (for diagnostics and output).
    const std::string& getName() const;

    /**
     * @brief Apply the process to the particles of one container for one step.
     *
     * @param element The beamline element this process is attached to (element-local frame).
     * @param pc      Particle container, transformed into the element-local frame.
     * @param context Step context (dt, global step, container index, frame transform).
     * @return Number of particles marked invalid via the container's InvalidMask.
     */
    virtual std::size_t apply(
            Component& element, ParticleContainer<double, 3>& pc,
            const LocalProcessContext& context) = 0;

private:
    /// Name of the definition that created this process.
    std::string name_m;
};

#endif
