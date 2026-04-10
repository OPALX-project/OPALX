//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef ALGORITHMS_ABSTRACTTIMEDEPENDENCE_H_
#define ALGORITHMS_ABSTRACTTIMEDEPENDENCE_H_

#include <map>
#include <memory>
#include <string>

/** @class AbstractTimeDependence
 *
 *  Time dependence abstraction for field parameters that vary slowly with time;
 *  for example, in RF cavities for synchrotrons the RF frequency varies to
 *  match the time-of-flight of the particles
 *
 *  This base class stores a mapping of string to time dependence. At Visit
 *  time, we do a map look-up to assign the appropriate TimeDependence to the
 *  relevant field Elements.
 */
class AbstractTimeDependence {
    // Unit tests for this class are in with PolynomialTimeDependence
public:
    /** Destructor does nothing */
    virtual ~AbstractTimeDependence() = default;

    /** Inheritable copy constructor
     *
     *  @returns new AbstractTimeDependence that is a copy of this. User owns
     *  returned memory.
     */
    virtual AbstractTimeDependence* clone() = 0;

    /** getValue(time) returns the value as a function of time.
     *
     *  This could represent RF voltage or frequency, magnetic field strength
     *  etc
     */
    virtual double getValue(double time) = 0;

    /** getIntegral(time) returns the integral from 0 to time of the function.
     *
     */
    virtual double getIntegral(double time) = 0;

    /** Look up the time dependence that has a given name
     *
     *  @param name name of the time dependence
     *
     *  @returns shared_ptr to the appropriate time dependence.
     *  @throws GeneralClassicException if name is not recognised
     */
    static std::shared_ptr<AbstractTimeDependence> getTimeDependence(const std::string& name);

    /** Add a value to the lookup table
     *
     *  @param name name of the time dependence. If name already exists in the
     *  map, it is overwritten with the new value.
     *  @param time_dep shared_ptr to the time dependence.
     */
    static void setTimeDependence(
            const std::string& name, std::shared_ptr<AbstractTimeDependence> time_dep);

    /** Get the name corresponding to a given time_dep
     *
     *  @param time_dep time dependence to lookup
     *
     *  @returns name corresponding to the time dependence. Note that this
     *  just does a dumb loop over the stored map values; so O(N).
     *  @throws GeneralClassicException if time_dep is not recognised
     */
    static std::string getName(const std::shared_ptr<AbstractTimeDependence>& time_dep);

private:
    static std::map<std::string, std::shared_ptr<AbstractTimeDependence> > td_map;
};

#endif
