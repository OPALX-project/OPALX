#include "Fields/G4BL2DMagnetoStatic.h"
#include "Fields/FM2DMagnetoStatic.h"  // reuse of computeField(), see the class doc
#include "PartBunch/PartBunch.h"
#include "Physics/Units.h"
#include "Utilities/GeneralOpalException.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

namespace {

    /**
     * @brief Read the next line that carries content.
     *
     * Strips '#' comments and surrounding whitespace and skips lines that are
     * left empty.
     *
     * @note This deliberately does not use Fieldmap::getLine(). That reads
     * through a 256 character buffer (READ_BUFFER_LENGTH), while a data row
     * here is nR values wide -- 51 values, about 660 characters, for the muE4
     * WSX map -- and would be silently truncated.
     *
     * @return false at end of file
     */
    bool nextLine(std::istream& in, std::string& line) {
        while (std::getline(in, line)) {
            const size_t comment = line.find('#');
            if (comment != std::string::npos) {
                line.erase(comment);
            }
            const size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) {
                continue;
            }
            const size_t last = line.find_last_not_of(" \t\r\n");
            line              = line.substr(first, last - first + 1);
            return true;
        }
        return false;
    }

    /// @brief First whitespace separated token of a line, the section keyword
    std::string firstToken(const std::string& line) {
        std::istringstream is(line);
        std::string token;
        is >> token;
        return token;
    }

    /// @brief The `key=value` pairs of a line, skipping the leading keyword
    std::map<std::string, std::string> keyValues(const std::string& line) {
        std::istringstream is(line);
        std::string token;
        is >> token;  // discard the section keyword

        std::map<std::string, std::string> pairs;
        while (is >> token) {
            const size_t equals = token.find('=');
            if (equals != std::string::npos) {
                pairs[token.substr(0, equals)] = token.substr(equals + 1);
            }
        }
        return pairs;
    }

}  // namespace

/**
 * @brief Constructor. Parses and validates the whole file without storing it.
 *
 * The grid bounds have to be known this early because Solenoid::initialise()
 * calls getFieldDimensions() straight after the factory, before goOnline()
 * triggers readMap().
 */
G4BL2DMagnetoStatic::G4BL2DMagnetoStatic(std::string aFilename, bool zReverse)
    : Fieldmap(aFilename),
      rbegin_m(0.0),
      rend_m(0.0),
      zbegin_m(0.0),
      zend_m(-1e-3),
      hz_m(0.0),
      hr_m(0.0),
      num_gridpr_m(0),
      num_gridpz_m(0),
      fieldScale_m(1.0) {
    Type        = TG4BL2DMagnetoStatic;
    zReverse_m  = zReverse;
    normalize_m = false;  // absolute Tesla, so KS is a plain multiplier

    parseFile(false);
}

G4BL2DMagnetoStatic::~G4BL2DMagnetoStatic() { freeMap(); }

void G4BL2DMagnetoStatic::parseFile(bool loadData) {
    std::ifstream in(Filename_m.c_str());
    if (!in.good()) {
        noFieldmapWarning();
        throw GeneralOpalException(
                "G4BL2DMagnetoStatic::parseFile", "Could not open fieldmap '" + Filename_m + "'");
    }

    std::string line;
    if (!nextLine(in, line)) {
        throw GeneralOpalException(
                "G4BL2DMagnetoStatic::parseFile", "Fieldmap '" + Filename_m + "' is empty");
    }

    // ---- optional `param` line ------------------------------------------
    fieldScale_m        = 1.0;
    std::string keyword = firstToken(line);
    if (keyword == "param") {
        const std::map<std::string, std::string> pairs = keyValues(line);
        auto optional = [&](const std::string& key, double fallback) -> double {
            const auto entry = pairs.find(key);
            if (entry == pairs.end()) {
                return fallback;
            }
            try {
                return std::stod(entry->second);
            } catch (const std::exception&) {
                throw GeneralOpalException(
                        "G4BL2DMagnetoStatic::parseFile", "Could not read " + key + "='"
                                                                  + entry->second + "' in '"
                                                                  + Filename_m + "'");
            }
        };

        // G4beamline evaluates B * normB * current_deck / current_file, where
        // current_file is the coil current the map was computed at. Fold the
        // file's half of that in here, so KS on the element is exactly
        // G4beamline's deck-side `current=` and KS = 1 matches a bare
        // `fieldmap ... ` placement. Both default to 1 in G4beamline.
        const double normB   = optional("normB", 1.0);
        const double current = optional("current", 1.0);
        if (current == 0.0) {
            throw GeneralOpalException(
                    "G4BL2DMagnetoStatic::parseFile",
                    "The 'param' line of '" + Filename_m
                            + "' sets current=0, which G4beamline would divide by");
        }
        fieldScale_m = normB / current;

        if (!nextLine(in, line)) {
            throw GeneralOpalException(
                    "G4BL2DMagnetoStatic::parseFile",
                    "Fieldmap '" + Filename_m + "' has no section after its 'param' line");
        }
        keyword = firstToken(line);
    }

    // ---- `cylinder` line -------------------------------------------------
    if (keyword != "cylinder") {
        throw GeneralOpalException(
                "G4BL2DMagnetoStatic::parseFile",
                "Expected a 'cylinder' section in '" + Filename_m + "', found '" + keyword + "'");
    }

    const std::map<std::string, std::string> grid = keyValues(line);
    auto require                                  = [&](const std::string& key) -> double {
        const auto entry = grid.find(key);
        if (entry == grid.end()) {
            throw GeneralOpalException(
                    "G4BL2DMagnetoStatic::parseFile",
                    "The 'cylinder' line of '" + Filename_m + "' has no " + key);
        }
        try {
            return std::stod(entry->second);
        } catch (const std::exception&) {
            throw GeneralOpalException(
                    "G4BL2DMagnetoStatic::parseFile",
                    "Could not read " + key + "='" + entry->second + "' in '" + Filename_m + "'");
        }
    };

    // There is no R0 in this format, r always starts at 0. X0, Y0 and
    // tolerance are accepted and ignored.
    const double Z0 = require("Z0");
    const double dZ = require("dZ");
    const double dR = require("dR");
    num_gridpz_m    = static_cast<int>(require("nZ"));
    num_gridpr_m    = static_cast<int>(require("nR"));

    if (num_gridpz_m < 2 || num_gridpr_m < 2) {
        throw GeneralOpalException(
                "G4BL2DMagnetoStatic::parseFile",
                "'" + Filename_m + "' needs at least 2 grid points in z and in r");
    }
    if (dZ <= 0.0 || dR <= 0.0) {
        throw GeneralOpalException(
                "G4BL2DMagnetoStatic::parseFile",
                "'" + Filename_m + "' has a non-positive dZ or dR");
    }

    hz_m     = dZ * Units::mm2m;
    hr_m     = dR * Units::mm2m;
    rbegin_m = 0.0;
    rend_m   = (num_gridpr_m - 1) * hr_m;
    zbegin_m = Z0 * Units::mm2m;
    zend_m   = zbegin_m + (num_gridpz_m - 1) * hz_m;

    if (zReverse_m) {
        const double zbegin = zbegin_m;
        zbegin_m            = -zend_m;
        zend_m              = -zbegin;
    }

    // ---- field blocks ----------------------------------------------------
    auto Bz = FieldstrengthBz_m.view_host();
    auto Br = FieldstrengthBr_m.view_host();

    bool haveBz = false;
    while (nextLine(in, line)) {
        const std::string label = firstToken(line);

        if (label != "Bz" && label != "Br") {
            std::string reason;
            if (label == "data") {
                reason = "the one-point-per-row 'data' block is not supported";
            } else if (label.compare(0, 6, "extend") == 0) {
                reason = "symmetry declarations ('" + label + "') are not supported";
            } else if (label == "Er" || label == "Ez" || label == "Bphi") {
                reason = "'" + label + "' blocks are not supported, only Bz and Br";
            } else {
                reason = "unexpected line '" + line + "'";
            }
            throw GeneralOpalException(
                    "G4BL2DMagnetoStatic::parseFile",
                    "In fieldmap '" + Filename_m + "': " + reason);
        }

        const bool isBz = (label == "Bz");
        haveBz          = haveBz || isBz;
        // rotation=Y180 mirrors the map in z and negates Bz; Br is unchanged.
        const double sign = (isBz && zReverse_m) ? -1.0 : 1.0;

        for (int j = 0; j < num_gridpz_m; ++j) {
            if (!nextLine(in, line)) {
                throw GeneralOpalException(
                        "G4BL2DMagnetoStatic::parseFile",
                        "Fieldmap '" + Filename_m + "' ends inside its " + label + " block, after "
                                + std::to_string(j) + " of " + std::to_string(num_gridpz_m)
                                + " rows");
            }

            const int jj = zReverse_m ? (num_gridpz_m - 1 - j) : j;

            std::istringstream row(line);
            for (int i = 0; i < num_gridpr_m; ++i) {
                double value = 0.0;
                if (!(row >> value)) {
                    throw GeneralOpalException(
                            "G4BL2DMagnetoStatic::parseFile",
                            "Row " + std::to_string(j) + " of the " + label + " block of '"
                                    + Filename_m + "' has fewer than "
                                    + std::to_string(num_gridpr_m) + " values");
                }
                if (loadData) {
                    const size_t index = static_cast<size_t>(i) * num_gridpz_m + jj;
                    if (isBz) {
                        Bz(index) = sign * fieldScale_m * value;
                    } else {
                        Br(index) = sign * fieldScale_m * value;
                    }
                }
            }

            std::string extra;
            if (row >> extra) {
                throw GeneralOpalException(
                        "G4BL2DMagnetoStatic::parseFile",
                        "Row " + std::to_string(j) + " of the " + label + " block of '" + Filename_m
                                + "' has more than " + std::to_string(num_gridpr_m) + " values");
            }
        }
    }

    if (!haveBz) {
        throw GeneralOpalException(
                "G4BL2DMagnetoStatic::parseFile", "Fieldmap '" + Filename_m + "' has no Bz block");
    }
    // A missing Br block leaves Br at zero, which is what G4beamline does too.
}

void G4BL2DMagnetoStatic::readMap() {
    if (FieldstrengthBz_m.extent(0) == 0) {
        const size_t size = static_cast<size_t>(num_gridpz_m) * num_gridpr_m;
        FieldstrengthBz_m = Kokkos::DualView<double*>("FieldstrengthBz", size);
        FieldstrengthBr_m = Kokkos::DualView<double*>("FieldstrengthBr", size);

        parseFile(true);

        FieldstrengthBz_m.modify<typename decltype(FieldstrengthBz_m)::host_mirror_space>();
        FieldstrengthBz_m.sync<typename decltype(FieldstrengthBz_m)::t_dev::device_type>();

        FieldstrengthBr_m.modify<typename decltype(FieldstrengthBr_m)::host_mirror_space>();
        FieldstrengthBr_m.sync<typename decltype(FieldstrengthBr_m)::t_dev::device_type>();

        *ippl::Info << level3 << typeset_msg("read in fieldmap '" + Filename_m + "'", "info")
                    << endl;
    }
}

void G4BL2DMagnetoStatic::freeMap() {
    if (FieldstrengthBz_m.extent(0) != 0) {
        FieldstrengthBz_m = Kokkos::DualView<double*>();
        FieldstrengthBr_m = Kokkos::DualView<double*>();

        *ippl::Info << level3 << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << endl;
    }
}

/**
 * @brief Apply the FM to all the particles.
 *
 * @param pc Particle container
 * @param scale Scaling factor applied to the field
 */
void G4BL2DMagnetoStatic::applyField(std::shared_ptr<ParticleContainer_t> pc, double scale) {
    // Local copies of member variables for use in the lambda function
    double zbegin  = zbegin_m;
    double zend    = zend_m;
    double rend    = rend_m;
    double hr      = hr_m;
    double hz      = hz_m;
    int num_gridpr = num_gridpr_m;
    int num_gridpz = num_gridpz_m;

    // Device accessible views
    auto Bz_device      = FieldstrengthBz_m.view_device();
    auto Br_device      = FieldstrengthBr_m.view_device();
    auto Rview          = pc->R.getView();
    auto Bview          = pc->B.getView();
    const size_t nLocal = pc->getLocalNum();

    Kokkos::parallel_for(
            "G4BL2DMagnetoStatic::applyField", nLocal, KOKKOS_LAMBDA(const size_t i) {
                // Check bounds
                if (Rview(i)(2) >= zbegin && Rview(i)(2) < zend
                    && sqrt(Rview(i)(0) * Rview(i)(0) + Rview(i)(1) * Rview(i)(1)) < rend) {
                    Vector_t<double, 3> tmpB = 0.0;
                    FM2DMagnetoStatic::computeField(
                            Rview(i), tmpB, Bz_device, Br_device, hr, hz, zbegin, num_gridpr,
                            num_gridpz);
                    Bview(i) += scale * tmpB;
                }
            });

    return;
}

/**
 * @brief Get the fieldstrength at position R
 *
 * @param R Position
 * @param E Electric Field (unused)
 * @param B Magnetic Field
 *
 * @return true if R is outside of the field map, false otherwise.
 */
bool G4BL2DMagnetoStatic::getFieldstrength(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B) const {
    if (isInside(R)) {
        FM2DMagnetoStatic::computeField(
                R, B, FieldstrengthBz_m.view_host(), FieldstrengthBr_m.view_host(), hr_m, hz_m,
                zbegin_m, num_gridpr_m, num_gridpz_m);
        return false;
    } else {
        return true;
    }
}

bool G4BL2DMagnetoStatic::getFieldDerivative(
        const Vector_t<double, 3>& /*R*/, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/,
        const DiffDirection& /*dir*/) const {
    throw GeneralOpalException("G4BL2DMagnetoStatic::getFieldDerivative", "not implemented");
}

void G4BL2DMagnetoStatic::getFieldDimensions(double& zBegin, double& zEnd) const {
    zBegin = zbegin_m;
    zEnd   = zend_m;
}

void G4BL2DMagnetoStatic::getFieldDimensions(
        double& xIni, double& xFinal, double& yIni, double& yFinal, double& zIni,
        double& zFinal) const {
    const double radius = std::max(std::abs(rbegin_m), std::abs(rend_m));
    xIni                = -radius;
    xFinal              = radius;
    yIni                = -radius;
    yFinal              = radius;
    zIni                = zbegin_m;
    zFinal              = zend_m;
}

void G4BL2DMagnetoStatic::swap() {
    // The block layout is fixed, so there is nothing to swap.
}

void G4BL2DMagnetoStatic::getInfo(Inform* msg) {
    (*msg) << Filename_m << " (G4beamline cylinder, 2D magnetostatic); zini= " << zbegin_m
           << " m; zfinal= " << zend_m << " m; rmax= " << rend_m << " m;"
           << (zReverse_m ? " reversed in z;" : "") << endl;
}

double G4BL2DMagnetoStatic::getFrequency() const {
    throw GeneralOpalException("G4BL2DMagnetoStatic::getFrequency", "not implemented");
    return 0.0;
}

void G4BL2DMagnetoStatic::setFrequency(double /*freq*/) {
    throw GeneralOpalException("G4BL2DMagnetoStatic::setFrequency", "not implemented");
    return;
}
