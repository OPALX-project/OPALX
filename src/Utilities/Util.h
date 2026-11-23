//
// Namespace Util
//   This namespace contains useful global methods.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef UTIL
#define UTIL

#include "Algorithms/Quaternion.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <string>
#include "OPALTypes.h"

// ------- DON'T DELETE: start --------
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define __DBGMSG__ __FILENAME__ << ": " << __LINE__ << "\t"
// ------- DON'T DELETE: end   --------

namespace Util {
    /// Directory, relative to the working directory, for the auxiliary output files of a run.
    inline constexpr const char* auxiliaryOutputDirectory = "data";

    std::string getGitRevision();

    inline double getGamma(ippl::Vector<double, 3> p) {
        double dotP = 0.0;  // \todo dot(p, p);  // .apply();
        for (unsigned i = 0; i < 3; i++)
            dotP += p(i) * p(i);
        return std::sqrt(dotP + 1.0);
    }

    inline ippl::Vector<double, 3> getBeta(ippl::Vector<double, 3> p) { return p / getGamma(p); }

    inline double getKineticEnergy(ippl::Vector<double, 3> p, double mass) {
        return (getGamma(p) - 1.0) * mass;
    }

    inline double getBetaGamma(double Ekin, double mass) {
        double value = std::sqrt(std::pow(Ekin / mass + 1.0, 2) - 1.0);
        if (value < std::numeric_limits<double>::epsilon()) value = std::sqrt(2 * Ekin / mass);
        return value;
    }

    inline std::string getTimeString(double time, unsigned int precision = 3) {
        std::string timeUnit(" [ps]");

        time *= 1e12;
        if (std::abs(time) > 1000) {
            time /= 1000;
            timeUnit = std::string(" [ns]");

            if (std::abs(time) > 1000) {
                time /= 1000;
                timeUnit = std::string(" [ms]");
            }
        } else if (std::abs(time) < 1.0) {
            time *= 1000;
            timeUnit = std::string(" [fs]");
        }

        std::stringstream timeOutput;
        timeOutput << std::fixed << std::setw(precision + 2) << std::setprecision(precision) << time
                   << timeUnit;
        return timeOutput.str();
    }

    inline std::string getLengthString(double spos, unsigned int precision = 3) {
        std::string sposUnit(" [m]");

        if (std::abs(spos) < 1.0) {
            spos *= 1000.0;
            sposUnit = std::string(" [mm]");
        }

        if (std::abs(spos) < 1.0) {
            spos *= 1000.0;
            sposUnit = std::string(" [um]");
        }

        std::stringstream positionOutput;
        positionOutput << std::fixed << std::setw(precision + 2) << std::setprecision(precision)
                       << spos << sposUnit;
        return positionOutput.str();
    }

    inline std::string getLengthString(ippl::Vector<double, 3> spos, unsigned int precision = 3) {
        std::string sposUnit(" [m]");
        double maxPos = std::abs(spos(0));
        for (unsigned int i = 1; i < 3u; ++i) {
            maxPos = std::max(maxPos, std::abs(spos(i)));
        }

        std::stringstream positionOutput;

        if (maxPos < 1.0) {
            maxPos *= 1000.0;
            spos     = spos * 1000.0;
            sposUnit = std::string(" [mm]");
        }

        if (maxPos < 1.0) {
            maxPos *= 1000.0;
            spos     = spos * 1000.0;
            sposUnit = std::string(" [um]");
        }

        positionOutput << std::fixed << std::setprecision(precision) << "( "
                       << std::setw(precision + 7) << spos(0) << " , " << std::setw(precision + 7)
                       << spos(1) << " , " << std::setw(precision + 7) << spos(2) << " )"
                       << sposUnit;
        return positionOutput.str();
    }

    inline std::string getEnergyString(double energyInMeV, unsigned int precision = 3) {
        std::string energyUnit(" [MeV]");
        double energy = energyInMeV;

        if (energy > 1000.0) {
            energy /= 1000.0;
            energyUnit = std::string(" [GeV]");
        } else if (energy < 1.0) {
            energy *= 1000.0;
            energyUnit = std::string(" [keV]");
            if (energy < 1.0) {
                energy *= 1000.0;
                energyUnit = std::string(" [eV]");
            }
        }

        std::stringstream energyOutput;
        energyOutput << std::fixed << std::setw(precision + 2) << std::setprecision(precision)
                     << energy << energyUnit;

        return energyOutput.str();
    }

    ippl::Vector<double, 3> getTaitBryantAngles(
            Quaternion rotation, const std::string& elementName = "");

    std::string toUpper(const std::string& str);

    std::string combineFilePath(std::initializer_list<std::string>);

    struct KahanAccumulation {
        long double sum;
        long double correction;
        KahanAccumulation();

        KahanAccumulation& operator+=(double value);
    };

    unsigned int rewindLinesSDDS(
            const std::string& fileName, double maxSPos, bool checkForTime = true);

    std::string base64_encode(
            const std::string& string_to_encode);  // unsigned char const* , unsigned int len);

    template <typename T, typename A>
    T* c_data(std::vector<T, A>& v) {
        return v.empty() ? static_cast<T*>(0) : &(v[0]);
    }

    template <typename T, typename A>
    T const* c_data(std::vector<T, A> const& v) {
        return v.empty() ? static_cast<T const*>(0) : &(v[0]);
    }
}  // namespace Util

#endif
