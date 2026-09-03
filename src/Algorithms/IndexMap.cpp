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
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <tuple>

#include "AbsBeamline/Multipole.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/IndexMap.h"
#include "Physics/Physics.h"
#include "Structure/ElementPositionWriter.h"
#include "Utilities/Util.h"

extern Inform* gmsg;

const double IndexMap::oneMinusEpsilon_m = 1.0 - std::numeric_limits<double>::epsilon();
namespace {
    void insertFlags(std::vector<double>& flags, std::shared_ptr<ElementBase> element);
}

IndexMap::IndexMap()
    : mapRange2Element_m(),
      mapElement2Range_m(),
      totalPathLength_m(0.0),
      periodOrigin_m(0.0),
      period_m(0.0) {}

void IndexMap::setPeriod(double origin, double period) {
    if (!(period > 0.0) || !std::isfinite(origin) || !std::isfinite(period)) {
        throw OpalException("IndexMap::setPeriod", "period must be finite and positive");
    }
    periodOrigin_m = origin;
    period_m       = period;
}

void IndexMap::print(std::ostream& out) const {
    if (mapRange2Element_m.empty()) return;

    out << "* Size of map " << mapRange2Element_m.size() << " sections " << std::endl;
    out << std::fixed << std::setprecision(6);
    auto mapIti = mapRange2Element_m.begin();
    auto mapItf = mapRange2Element_m.end();

    double totalLength     = (*mapRange2Element_m.rbegin()).first.end;
    unsigned int numDigits = std::floor(std::max(0.0, log(totalLength) / log(10.0))) + 1;

    for (; mapIti != mapItf; mapIti++) {
        const key_t key   = (*mapIti).first;
        const value_t val = (*mapIti).second;
        out << "* Key: (" << std::setw(numDigits + 7) << std::right << key.begin << " - "
            << std::setw(numDigits + 7) << std::right << key.end
            << ") number of overlapping elements " << val.size() << "\n";

        for (auto element : val) {
            out << "* " << std::setw(25 + 2 * numDigits) << " " << element->getName() << "\n";
        }
    }
}

IndexMap::value_t IndexMap::query(key_t::first_type s, key_t::second_type ds) {
    if (ds < 0.0) {
        throw OutOfBounds("IndexMap::query", "query half-width must not be negative");
    }
    if (period_m == 0.0) {
        return queryRange(s - ds, s + ds);
    }

    if (2.0 * ds >= period_m) {
        value_t elements;
        for (const auto& entry : mapRange2Element_m) {
            const auto& rangeElements = entry.second;
            elements.insert(rangeElements.begin(), rangeElements.end());
        }
        return elements;
    }

    const double lower        = s - ds;
    const double upper        = s + ds;
    const double wrappedLower = wrap(lower);
    const double wrappedUpper = wrap(upper);
    const auto lowerTurn = static_cast<long long>(std::floor((lower - periodOrigin_m) / period_m));
    const auto upperTurn = static_cast<long long>(std::floor((upper - periodOrigin_m) / period_m));

    if (lowerTurn == upperTurn && wrappedLower <= wrappedUpper) {
        return queryRange(wrappedLower, wrappedUpper);
    }

    value_t elements = queryRange(wrappedLower, periodOrigin_m + period_m);
    value_t head     = queryRange(periodOrigin_m, wrappedUpper);
    elements.insert(head.begin(), head.end());
    return elements;
}

IndexMap::value_t IndexMap::queryRange(double lowerLimit, double upperLimit) {
    upperLimit = std::min(totalPathLength_m, upperLimit);
    value_t elementSet;

    map_t::reverse_iterator rit = mapRange2Element_m.rbegin();
    if (rit != mapRange2Element_m.rend() && lowerLimit > (*rit).first.end) {
        throw OutOfBounds("IndexMap::query", "out of bounds");
    }

    // Start the scan at the first entry that can overlap [lowerLimit, upperLimit]. Entries before
    // lower_bound have begin < lowerLimit; step back over those still reaching into the interval.
    map_t::iterator it        = mapRange2Element_m.lower_bound(key_t{lowerLimit, lowerLimit});
    const map_t::iterator end = mapRange2Element_m.end();
    while (it != mapRange2Element_m.begin() && std::prev(it)->first.end > lowerLimit) {
        --it;
    }

    for (; it != end; ++it) {
        const double low  = (*it).first.begin;
        const double high = (*it).first.end;

        if (upperLimit < low) return elementSet;  // all following entries begin even later
        if (lowerLimit < high) break;
    }

    if (it == end) return elementSet;

    map_t::iterator last = std::next(it);
    for (; last != end; ++last) {
        const double low = (*last).first.begin;

        if (upperLimit < low) break;
    }

    for (; it != last; ++it) {
        const value_t& a = (*it).second;
        elementSet.insert(a.cbegin(), a.cend());
    }

    return elementSet;
}

double IndexMap::wrap(double s) const {
    double offset = std::fmod(s - periodOrigin_m, period_m);
    if (offset < 0.0) {
        offset += period_m;
    }
    return periodOrigin_m + offset;
}

long long IndexMap::getTurnNumber(double s) const {
    if (period_m == 0.0) {
        throw OpalException("IndexMap::getTurnNumber", "map is not periodic");
    }
    return static_cast<long long>(std::floor((s - periodOrigin_m) / period_m));
}

double IndexMap::getPhase(double s) const {
    if (period_m == 0.0) {
        throw OpalException("IndexMap::getPhase", "map is not periodic");
    }
    const double twoPi = 2.0 * std::acos(-1.0);
    return twoPi * (wrap(s) - periodOrigin_m) / period_m;
}

void IndexMap::add(key_t::first_type initialS, key_t::second_type finalS, const value_t& val) {
    if (initialS > finalS) {
        std::swap(initialS, finalS);
    }
    key_t key{initialS, finalS * oneMinusEpsilon_m};

    mapRange2Element_m.insert(std::pair<key_t, value_t>(key, val));
    totalPathLength_m = (*mapRange2Element_m.rbegin()).first.end;

    // Build each element's first-crossing range in the reverse map. The reference orbit threads
    // through an element over several consecutive steps, each a separate add(); merge them into
    // one contiguous interval. A non-contiguous add for a known element means the orbit
    // re-entered it (ring); the first crossing is kept and the new range is only recorded in
    // mapRange2Element_m, which query() reads.
    for (value_t::iterator setIt = val.begin(); setIt != val.end(); ++setIt) {
        auto res = mapElement2Range_m.insert(std::make_pair(*setIt, key));
        if (!res.second) {
            key_t& currentRange = res.first->second;
            if (almostEqual(key.begin, currentRange.end / oneMinusEpsilon_m)) {
                currentRange.end = key.end;  // extend the element's contiguous range
            }
        }
    }
}

void IndexMap::tidyUp(double sStop) {
    map_t::reverse_iterator rit = mapRange2Element_m.rbegin();

    if (rit != mapRange2Element_m.rend() && (*rit).second.empty() && sStop > (*rit).first.begin) {
        key_t key{(*rit).first.begin, sStop};
        value_t val;

        mapRange2Element_m.erase(std::next(rit).base());
        mapRange2Element_m.insert(std::pair<key_t, value_t>(key, val));
    }
}

enum elements {
    DIPOLE = 0,
    QUADRUPOLE,
    SEXTUPOLE,
    OCTUPOLE,
    DECAPOLE,
    MULTIPOLE,
    SOLENOID,
    RFCAVITY,
    MONITOR,
    OTHER,
    SIZE
};

void IndexMap::saveSDDS(double initialPathLength) const {
    if (mapRange2Element_m.empty()) return;

    std::vector<std::tuple<double, std::vector<double>, std::string> > sectors;

    // add for each sector four rows:
    // (s_i, 0)
    // (s_i, 1)
    // (s_f, 1)
    // (s_f, 0)
    // to the file, where
    // s_i is the start of the range and
    // s_f is the end of the range.
    // Each range-map entry {range, elements} is one sector of the reference path. On a ring
    // only the first turn is written: the reverse map holds first-crossing ranges only, so a
    // sector containing an element whose stored range ends before the sector begins is a
    // re-entry (second turn) — stop there.
    for (const auto& [range, sectorElements] : mapRange2Element_m) {
        if (sectorElements.empty()) continue;

        const bool reentered =
                std::any_of(sectorElements.begin(), sectorElements.end(), [&](const auto& element) {
                    return mapElement2Range_m.at(element).end < range.begin;
                });
        if (reentered) {
            *gmsg << level2
                  << "* IndexMap: ring detected (element re-entered at s = " << range.begin
                  << " m); the element position file contains the first pass only" << endl;
            break;
        }

        double sectorBegin = range.begin;
        double sectorEnd   = range.end;

        std::vector<std::tuple<double, std::vector<double>, std::string> > currentSector(4);
        std::get<0>(currentSector[0]) = sectorBegin;
        std::get<0>(currentSector[1]) = sectorBegin;
        std::get<0>(currentSector[2]) = sectorEnd;
        std::get<0>(currentSector[3]) = sectorEnd;

        for (unsigned short i = 0; i < 4; ++i) {
            auto& flags = std::get<1>(currentSector[i]);
            flags.resize(SIZE, 0);
        }

        for (auto element : sectorElements) {
            const auto& elementRange = mapElement2Range_m.at(element);
            if (elementRange.begin < sectorBegin) {
                ::insertFlags(std::get<1>(currentSector[0]), element);
                std::get<2>(currentSector[0]) += element->getName() + ", ";
            }

            ::insertFlags(std::get<1>(currentSector[1]), element);
            std::get<2>(currentSector[1]) += element->getName() + ", ";

            ::insertFlags(std::get<1>(currentSector[2]), element);
            std::get<2>(currentSector[2]) += element->getName() + ", ";

            if (elementRange.end > sectorEnd) {
                ::insertFlags(std::get<1>(currentSector[3]), element);
                std::get<2>(currentSector[3]) += element->getName() + ", ";
            }
        }

        for (unsigned short i = 0; i < 4; ++i) {
            sectors.push_back(currentSector[i]);
        }
    }

    // make the entries of the rf cavities a zigzag line
    const unsigned int numEntries = sectors.size();
    auto it                       = mapElement2Range_m.begin();
    auto end                      = mapElement2Range_m.end();
    for (; it != end; ++it) {
        auto element = (*it).first;
        auto name    = element->getName();
        auto type    = element->getType();
        if (type != ElementType::RFCAVITY && type != ElementType::TRAVELINGWAVE) {
            continue;
        }

        auto range = (*it).second;

        unsigned int i = 0;
        for (; i < numEntries; ++i) {
            if (std::get<0>(sectors[i]) >= range.begin) {
                break;
            }
        }

        if (i == numEntries) continue;

        unsigned int j = ++i;
        while (j < numEntries && std::get<0>(sectors[j]) < range.end) {
            ++j;
        }
        // Bound fix since we can go beyond the range when threading through
        // elements multiple times
        if (j == numEntries) {
            --j;
        }

        double length = range.end - range.begin;
        for (; i <= j; ++i) {
            double pos  = std::get<0>(sectors[i]);
            auto& items = std::get<1>(sectors[i]);

            items[RFCAVITY] = 1.0 - 2 * (pos - range.begin) / length;
        }
    }

    // add row if range of first sector starts after initialPathLength
    if (!sectors.empty() && std::get<0>(sectors[0]) > initialPathLength) {
        auto tmp = sectors;
        sectors  = std::vector<std::tuple<double, std::vector<double>, std::string> >(1);
        std::get<0>(sectors[0]) = initialPathLength;
        std::get<1>(sectors[0]).resize(SIZE, 0.0);

        sectors.insert(sectors.end(), tmp.begin(), tmp.end());
    }

    std::string fileName = Util::combineFilePath(
            {OpalData::getInstance()->getAuxiliaryOutputDirectory(),
             OpalData::getInstance()->getInputBasename() + "_ElementPositions.sdds"});
    ElementPositionWriter writer(fileName);

    for (auto sector : sectors) {
        std::string names = std::get<2>(sector);
        if (!names.empty()) {
            names = names.substr(0, names.length() - 2);
        }
        names = "\"" + names + "\"";
        writer.addRow(std::get<0>(sector), std::get<1>(sector), names);
    }
}

namespace {
    void insertFlags(std::vector<double>& flags, std::shared_ptr<ElementBase> element) {
        switch (element->getType()) {
            case ElementType::MULTIPOLE: {
                const Multipole* mult = static_cast<const Multipole*>(element.get());
                switch (mult->getMaxNormalComponentIndex()) {
                    case 1:
                        flags[DIPOLE] = (mult->isFocusing(0) ? 1 : -1);
                        break;
                    case 2:
                        flags[QUADRUPOLE] = (mult->isFocusing(1) ? 1 : -1);
                        break;
                    case 3:
                        flags[SEXTUPOLE] = (mult->isFocusing(2) ? 1 : -1);
                        break;
                    case 4:
                        flags[OCTUPOLE] = (mult->isFocusing(3) ? 1 : -1);
                        break;
                    case 5:
                        flags[DECAPOLE] = (mult->isFocusing(4) ? 1 : -1);
                        break;
                    default:
                        flags[MULTIPOLE] = 1;
                }
            } break;
            case ElementType::RFCAVITY:
            case ElementType::TRAVELINGWAVE:
                flags[RFCAVITY] = 1;
                break;
            case ElementType::SOLENOID:
                flags[SOLENOID] = 1;
                break;
            case ElementType::CONSTANTEFIELDCAVITY:
            case ElementType::CONSTANTFOCUSING:
                flags[OTHER] = 1;
                break;
            default:
                flags[OTHER] = 1;
                break;
        }
    }
}  // namespace

IndexMap::key_t IndexMap::getRange(const IndexMap::value_t::value_type& element) const {
    invertedMap_t::const_iterator it = mapElement2Range_m.find(element);
    if (it == mapElement2Range_m.end()) {
        throw OpalException(
                "IndexMap::getRange()", "Element \"" + element->getName() + "\" not registered");
    }
    return it->second;
}

bool IndexMap::almostEqual(double x, double y) {
    return (std::abs(x - y) < std::numeric_limits<double>::epsilon() * std::abs(x + y) * 2
            || std::abs(x - y) < std::numeric_limits<double>::min());
}
