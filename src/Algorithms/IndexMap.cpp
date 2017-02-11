#include "Algorithms/IndexMap.h"
#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Bend.h"
#include "Utilities/Util.h"
#include "Physics/Physics.h"
#include "OPALconfig.h"

#include <boost/filesystem.hpp>

#include <limits>
#include <iostream>
#include <fstream>

extern Inform *gmsg;

IndexMap::IndexMap():
    map_m(),
    totalPathLength_m(0.0),
    oneMinusEpsilon_m(1.0 - std::numeric_limits<double>::epsilon())
{ }

void IndexMap::print(std::ostream &out) const {
    out << "Size of map " << map_m.size() << " sections " << std::endl;
    out << std::fixed << std::setprecision(6);
    auto mapIti = map_m.begin();
    auto mapItf = map_m.end();
    for (; mapIti != mapItf; mapIti++) {
        const auto key = (*mapIti).first;
        const auto val = (*mapIti).second;
        out << "Key: (" << std::setw(14) << key.first << "-" << std::setw(14) << key.second << ") number of overlapping elements " << val.size() << "\n";
        for (auto element: val) {
            out << std::setw(37) << " " << element->getName() << "\n";
        }
    }
}

IndexMap::value_t IndexMap::query(key_t::first_type s, key_t::second_type ds) {
    const double lowerLimit = (ds < s? s - ds: 0);
    const double upperLimit = std::min(totalPathLength_m, s + ds);
    value_t elementSet;

    map_t::reverse_iterator rit = map_m.rbegin();
    if (rit != map_m.rend() && lowerLimit > (*rit).first.second) {
        throw OutOfBounds("IndexMap::query", "out of bounds");
    }

    map_t::iterator it = map_m.begin();
    const map_t::iterator end = map_m.end();

    for (; it != end; ++ it) {
        const double low = (*it).first.first;
        const double high = (*it).first.second;

        if (lowerLimit < high && upperLimit >= low) break;
    }

    if (it == end) return elementSet;

    map_t::iterator last = std::next(it);
    for (; last != end; ++ last) {
        const double low = (*last).first.first;

        if (upperLimit < low) break;
    }

    for (; it != last; ++ it) {
        const value_t &a = (*it).second;
        elementSet.insert(a.cbegin(), a.cend());
    }

    return elementSet;
}

void IndexMap::tidyUp() {
    map_t::reverse_iterator rit = map_m.rbegin();

    if (rit != map_m.rend() && (*rit).second.size() == 0) {
        map_m.erase(std::next(rit).base());
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
    SIZE
};

void IndexMap::saveSDDS(double startS) const {
    std::string fileName(OpalData::getInstance()->getInputBasename() + "_ElementPositions.sdds");
    std::ofstream sdds;
    if (OpalData::getInstance()->hasPriorTrack() && boost::filesystem::exists(fileName)) {
        Util::rewindLinesSDDS(fileName, startS, false);
        sdds.open(fileName, std::ios::app);
    } else {
        sdds.open(fileName);

        std::string indent("        ");

        sdds << "SDDS1\n"
             << "&description \n"
             << indent << "text=\"Element positions '" << OpalData::getInstance()->getInputFn() << "'\", \n"
             << indent << "contents=\"element positions\" \n"
             << "&end\n";
        sdds << "&parameter\n"
             << indent << "name=revision, \n"
             << indent << "type=string,\n"
             << indent << "description=\"git revision of opal\"\n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=s, \n"
             << indent << "type=double, \n"
             << indent << "units=m, \n"
             << indent << "description=\"1 longitudinal position\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=dipole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"2 dipole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=quadrupole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"3 quadrupole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=sextupole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"4 sextupole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=octupole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"5 octupole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=decapole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"6 decapole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=multipole, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"7 higher multipole field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=solenoid, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"8 solenoid field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=rfcavity, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"9 RF field present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=monitor, \n"
             << indent << "type=float, \n"
             << indent << "units=1, \n"
             << indent << "description=\"10 monitor present\" \n"
             << "&end\n";
        sdds << "&column \n"
             << indent << "name=element_names, \n"
             << indent << "type=string, \n"
             << indent << "description=\"names of elements\" \n"
             << "&end\n";
        sdds << "&data \n"
             << indent << "mode=ascii, \n"
             << indent << "no_row_counts=1 \n"
             << "&end\n";

        sdds << PACKAGE_NAME << " " << PACKAGE_VERSION_STR << " git rev. " << Util::getGitRevision() << std::endl;
    }

    std::vector<std::vector<int> > allItems(SIZE);
    std::vector<double> allPositions;
    std::vector<std::string> allNames;
    std::vector<double> typeMultipliers = {3.3333e-1, 1.0, 0.5, 0.25, 1.0, 1.0, 1.0, 1.0, 1.0};
    unsigned int counter = 0;

    auto mapIti = map_m.begin();
    auto mapItf = map_m.end();
    for (; mapIti != mapItf; mapIti++) {
        const auto key = (*mapIti).first;
        const auto val = (*mapIti).second;

        std::vector<int> items(SIZE, 0);
        std::string names("");
        for (auto element: val) {
            switch (element->getType()) {
            case ElementBase::RBEND:
            case ElementBase::SBEND:
                {
                    const Bend* bend = static_cast<const Bend*>(element.get());
                    if (bend->getRotationAboutZ() > 0.5 * Physics::pi &&
                        bend->getRotationAboutZ() < 1.5 * Physics::pi) {
                        items[DIPOLE] = -1;
                    } else {
                        items[DIPOLE] = 1;
                    }
                }
                break;
            case ElementBase::MULTIPOLE:
                {
                    const Multipole* mult = static_cast<const Multipole*>(element.get());
                    switch(mult->getMaxNormalComponentIndex()) {
                    case 1:
                        items[DIPOLE] = (mult->isFocusing(0)? 1: -1);
                        break;
                    case 2:
                        items[QUADRUPOLE] = (mult->isFocusing(1)? 1: -1);
                        break;
                    case 3:
                        items[SEXTUPOLE] = (mult->isFocusing(2)? 1: -1);
                        break;
                    case 4:
                        items[OCTUPOLE] = (mult->isFocusing(3)? 1: -1);
                        break;
                    case 5:
                        items[DECAPOLE] = (mult->isFocusing(4)? 1: -1);
                        break;
                    default:
                        items[MULTIPOLE] = 1;
                    }
                }
                break;
            case ElementBase::SOLENOID:
                items[SOLENOID] = 1;
                break;
            case ElementBase::RFCAVITY:
            case ElementBase::TRAVELINGWAVE:
                items[RFCAVITY] = 1;
                break;
            case ElementBase::MONITOR:
                items[MONITOR] = 1;
                break;
            default:
                continue;
            }

            names += element->getName() + ", ";
        }

        if (names.length() == 0) continue;

        allPositions.push_back(key.first);
        allPositions.push_back(key.first);
        allPositions.push_back(key.second);
        allPositions.push_back(key.second);

        allNames.push_back("");
        allNames.push_back(names.substr(0, names.length() - 2));
        allNames.push_back(allNames.back());
        allNames.push_back("");

        for (unsigned int i = 0; i < SIZE; ++ i) {
            allItems[i].push_back(0);
            allItems[i].push_back(items[i]);
            allItems[i].push_back(items[i]);
            allItems[i].push_back(0);
        }

        ++ counter;
    }

    if (counter == 0) return;

    const unsigned int totalSize = counter;
    for (unsigned int i = 0; i < SIZE; ++ i) {
        for (unsigned int j = 0; j < totalSize - 1; ++ j) {
            if (allItems[i][4 * j + 1] == 0) continue;
            if (allItems[i][4 * (j + 1) + 1] == 0) continue;
            if (allItems[i][4 * j + 1] * allItems[i][4 * (j + 1) + 1] < 0) continue;
            if (std::abs(allPositions[4 * j + 3] - allPositions[4 * (j + 1)]) > 1e-6) continue;

            allItems[i][4 * j + 3] = allItems[i][4 * (j + 1)] = allItems[i][4 * j + 1];
        }
    }

    for (unsigned int j = 0; j < totalSize - 1; ++ j) {
        allItems[RFCAVITY][4 * j + 2] = -allItems[RFCAVITY][4 * j + 1];
    }

    for (unsigned int i = 0; i < 4 * totalSize; ++ i) {
        sdds << std::setw(16) << std::setprecision(8) << allPositions[i];

        sdds << std::setw(10) << std::setprecision(5) << allItems[0][i] * typeMultipliers[0];
        for (unsigned int j = 1; j < SIZE; ++ j)
            sdds << std::setw(6) << std::setprecision(3) << allItems[j][i] * typeMultipliers[j];
        sdds << " \"" << allNames[i] << "\"\n";
    }

    sdds.close();
}