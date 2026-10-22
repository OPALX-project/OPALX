//
// Class IndexMap
//
// This class stores and prints the sequence of elements that the referenc particle passes.
// Each time the reference particle enters or leaves an element an entry is added to the map.
// With help of this map one can determine which element can be found at a given position.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020 Christof Metzger-Kraus
//
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
#ifndef OPAL_INDEXMAP_H
#define OPAL_INDEXMAP_H

#include <map>
#include <ostream>

#include "AbsBeamline/ElementBase.h"
#include "Utilities/OpalException.h"

#include <set>
#include <utility>

/**
 * @brief Path-coordinate index of support-selected runtime element sets.
 *
 * OrbitThreader records which element supports were encountered along the reference,
 * so tracking can obtain candidate elements for a bunch-length interval. This is not
 * a map of nominal body owners, logical BeamlineMembership, or 6x6 transfer matrices.
 * A selected set may include field-free drifts and overlapping fringe supports.
 *
 * A periodic index stores one measured reference-return interval and wraps later
 * queries onto it. Its period is a path length [m], not the sum of support extents
 * or necessarily Ring::getLength(). Reusing the index does not establish orbit closure.
 */
class IndexMap {
public:
    struct Range {
        typedef double first_type;
        typedef double second_type;
        first_type begin;
        second_type end;
    };
    typedef Range key_t;
    typedef std::set<std::shared_ptr<ElementBase>> value_t;

    IndexMap();

    /// Configure query wrapping onto [origin, origin + period), both lengths in metres.
    /// Does not change stored ranges; the caller must supply matching one-period data.
    /// @throws OpalException If origin/period is nonfinite or period is not positive.
    void setPeriod(double origin, double period);

    /// Record support candidates on a path interval [m]; reversed endpoints are swapped
    /// and zero-length intervals ignored. Consecutive ranges with identical sets are
    /// coalesced; noncontiguous re-entries remain separate. Reverse lookup keeps the first visit.
    void add(key_t::first_type initialStep, key_t::second_type finalStep, const value_t& val);

    /// Union of candidates intersecting [s-ds, s+ds], where ds is a nonnegative half-width [m].
    /// Periodic queries crossing the seam join tail/head sets; widths spanning a full
    /// period return all candidates. This selects supports, not nonzero field values.
    /// @throws OutOfBounds For negative ds or a nonperiodic query beyond the stored range.
    value_t query(key_t::first_type s, key_t::second_type ds);

    /// floor((s-origin)/period), with turn zero starting at origin and negative turns allowed.
    /// @throws OpalException If no period was configured.
    long long getTurnNumber(double s) const;

    /// Path fraction times 2*pi in [0, 2*pi), not a geometrical azimuth or betatron phase.
    /// @param s Absolute reference-path coordinate [m].
    /// @throws OpalException If no period was configured.
    double getPhase(double s) const;

    void tidyUp(double sStop);

    void print(std::ostream&) const;
    void saveSDDS(double startS) const;
    size_t size() const;

    /// Reverse lookup: the path-length range of the FIRST crossing of @p element. Later
    /// crossings (ring) are recorded only in the range-to-elements map used by query().
    key_t getRange(const IndexMap::value_t::value_type& element) const;

    class OutOfBounds : public OpalException {
    public:
        OutOfBounds(const std::string& meth, const std::string& msg) : OpalException(meth, msg) {}

        OutOfBounds(const OutOfBounds& rhs) : OpalException(rhs) {}

        virtual ~OutOfBounds() {}

    private:
        OutOfBounds();
    };

private:
    class myCompare {
    public:
        bool operator()(const key_t& x, const key_t& y) const {
            if (x.begin < y.begin) return true;

            if (x.begin == y.begin) {
                if (x.end < y.end) return true;
            }

            return false;
        }
    };

    typedef std::map<key_t, value_t, myCompare> map_t;
    // Reverse map: one range per element — the first crossing (see getRange).
    typedef std::map<value_t::value_type, key_t, std::owner_less<value_t::value_type>>
            invertedMap_t;
    map_t mapRange2Element_m;
    invertedMap_t mapElement2Range_m;

    double totalPathLength_m;
    double periodOrigin_m;
    double period_m;

    value_t queryRange(double lowerLimit, double upperLimit);
    double wrap(double s) const;

    static bool almostEqual(double, double);
    static const double oneMinusEpsilon_m;
};

inline size_t IndexMap::size() const { return mapRange2Element_m.size(); }

inline std::ostream& operator<<(std::ostream& out, const IndexMap& im) {
    im.print(out);
    return out;
}

inline Inform& operator<<(Inform& out, const IndexMap& im) {
    im.print(out.getStream());
    return out;
}

#endif
