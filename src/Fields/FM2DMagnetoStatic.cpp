#include "Fields/FM2DMagnetoStatic.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Units.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <cmath>
#include <fstream>
#include <ios>

/**
 * @brief Constructor for 2D magnetostatic field map.
 * Parses the file header to read grid parameters:
 * - Orientation (XZ or ZX)
 * - Grid boundaries and number of points (zbegin_m, zend_m, rbegin_m, rend_m)
 * - Unit conversion from cm to m
 * - Grid spacing (hz_m, hr_m)
 * - Number of gridpoints (num_gridpz_m, num_gridpr_m)
 */
FM2DMagnetoStatic::FM2DMagnetoStatic(std::string aFilename) : Fieldmap(aFilename) {
    std::ifstream file;
    std::string tmpString;
    double tmpDouble;

    Type = T2DMagnetoStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if (file.good()) {
        bool parsing_passed = true;
        try {
            parsing_passed = interpretLine<std::string, std::string>(file, tmpString, tmpString);
        } catch (GeneralClassicException& e) {
            parsing_passed = interpretLine<std::string, std::string, std::string>(
                file, tmpString, tmpString, tmpString
            );

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" && tmpString != "FALSE")
                throw GeneralClassicException(
                    "FM2DMagnetoStatic::FM2DMagnetoStatic",
                    "The third string on the first line of 2D field "
                    "maps has to be either TRUE or FALSE"
                );

            normalize_m = (tmpString == "TRUE");
        }

        if (tmpString == "ZX") {
            swap_m = true;
            /// Parse rbegin_m, rend_m and num_gridpr_m(-1)
            parsing_passed =
                parsing_passed
                && interpretLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
            /// Parse rbegin_m, zend_m and num_gridpz_m(-1)
            parsing_passed =
                parsing_passed
                && interpretLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        } else if (tmpString == "XZ") {
            swap_m = false;
            /// Parse rbegin_m, zend_m and num_gridpz_m(-1)
            parsing_passed =
                parsing_passed
                && interpretLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
            /// Parse rbegin_m, rend_m and num_gridpr_m(-1)
            parsing_passed =
                parsing_passed
                && interpretLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        } else {
            std::cerr << "unknown orientation of 2D magnetostatic fieldmap" << std::endl;
            parsing_passed = false;
        }

        for (long i = 0; (i < (num_gridpz_m + 1) * (num_gridpr_m + 1)) && parsing_passed; ++i) {
            parsing_passed =
                parsing_passed && interpretLine<double, double>(file, tmpDouble, tmpDouble);
        }

        parsing_passed = parsing_passed && interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if (!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralClassicException(
                "FM2DMagnetoStatic::FM2DMagnetoStatic",
                "An error occured when reading the fieldmap '" + Filename_m + "'"
            );
        } else {
            // conversion from cm to m
            rbegin_m *= Units::cm2m;
            rend_m *= Units::cm2m;
            zbegin_m *= Units::cm2m;
            zend_m *= Units::cm2m;

            hr_m = (rend_m - rbegin_m) / num_gridpr_m;
            hz_m = (zend_m - zbegin_m) / num_gridpz_m;

            // num spacings -> num grid points
            num_gridpr_m++;
            num_gridpz_m++;
        }
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m   = -1e-3;
    }
}

FM2DMagnetoStatic::~FM2DMagnetoStatic() { freeMap(); }

void FM2DMagnetoStatic::readMap() {
    if (FieldstrengthBz_m.h_view.data() == nullptr) {
        // declare variables and allocate memory
        std::ifstream in;
        std::string tmpString;

        const size_t size = num_gridpz_m * num_gridpr_m;
        FieldstrengthBz_m = Kokkos::DualView<double*>("FieldstrengthBz", size);
        FieldstrengthBr_m = Kokkos::DualView<double*>("FieldstrengthBr", size);

        auto Bz = FieldstrengthBz_m.view_host();
        auto Br = FieldstrengthBr_m.view_host();

        // read in and parse field map
        in.open(Filename_m.c_str());
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);

        if (swap_m) {
            for (int i = 0; i < num_gridpz_m; i++) {
                for (int j = 0; j < num_gridpr_m; j++) {
                    interpretLine<double, double>(
                        in,                        // input stream
                        Br(j * num_gridpz_m + i),  // radial component
                        Bz(j * num_gridpz_m + i)   // longitudinal component
                    );
                }
            }
        } else {
            for (int j = 0; j < num_gridpr_m; j++) {
                for (int i = 0; i < num_gridpz_m; i++) {
                    interpretLine<double, double>(
                        in,                        // input stream
                        Bz(j * num_gridpz_m + i),  // longitudinal component
                        Br(j * num_gridpz_m + i)   // radial component
                    );
                }
            }
        }
        in.close();

        if (normalize_m) {
            double Bzmax = 0.0;
            // find maximum field
            for (int i = 0; i < num_gridpz_m; ++i) {
                if (std::abs(Bz(i)) > Bzmax) {
                    Bzmax = std::abs(Bz(i));
                }
            }

            // normalize field
            for (size_t i = 0; i < size; ++i) {
                Bz(i) /= Bzmax;
                Br(i) /= Bzmax;
            }
        }

        FieldstrengthBz_m.modify<Kokkos::HostSpace>();
        FieldstrengthBz_m.sync<Kokkos::DefaultExecutionSpace>();
        FieldstrengthBr_m.modify<Kokkos::HostSpace>();
        FieldstrengthBr_m.sync<Kokkos::DefaultExecutionSpace>();

        *ippl::Info << level3 << typeset_msg("read in fieldmap '" + Filename_m + "'", "info")
                    << endl;
    }
}

void FM2DMagnetoStatic::freeMap() {
    if (FieldstrengthBz_m.h_view.data() != nullptr) {
        FieldstrengthBz_m = Kokkos::DualView<double*>();
        FieldstrengthBr_m = Kokkos::DualView<double*>();

        *ippl::Info << level3 << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << endl;
    }
}

/**
 * @brief Apply the FM to all the particles.
 *
 * @param pc Particle container
 */
void FM2DMagnetoStatic::applyField(std::shared_ptr<ParticleContainer_t> pc) {
    // Local copies of member variables for use in the lambda function
    double zbegin  = zbegin_m;
    double zend    = zend_m;
    double rend    = rend_m;
    double hr      = hr_m;
    double hz      = hz_m;
    int num_gridpr = num_gridpr_m;
    int num_gridpz = num_gridpz_m;

    // Device accessible views
    auto Bz_device = FieldstrengthBz_m.view_device();
    auto Br_device = FieldstrengthBr_m.view_device();
    auto Rview     = pc->R.getView();
    auto Bview     = pc->B.getView();

    Kokkos::parallel_for(
        "FM2DMagnetoStatic::applyField", ippl::getRangePolicy(Rview), KOKKOS_LAMBDA(const int i) {
            // Check bounds
            if (Rview(i)(2) >= zbegin && Rview(i)(2) < zend
                && sqrt(Rview(i)(0) * Rview(i)(0) + Rview(i)(1) * Rview(i)(1)) < rend) {
                computeField(
                    Rview(i), Bview(i), Bz_device, Br_device, hr, hz, zbegin, num_gridpr, num_gridpz
                );
            }
        }
    );

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
bool FM2DMagnetoStatic::getFieldstrength(
    const Vector_t<double, 3>& R, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B
) const {
    if (isInside(R)) {
        computeField(
            R, B, FieldstrengthBz_m.h_view, FieldstrengthBr_m.h_view, hr_m, hz_m, zbegin_m,
            num_gridpr_m, num_gridpz_m
        );
        return false;
    } else {
        return true;
    }
}

/**
 * @brief Get the derivative of the field at position R
 *
 * @param R Position
 * @param E Electric Field (unused)
 * @param B Derivate of the magnetic field
 *
 * @note Not implemented yet
 *
 */
bool FM2DMagnetoStatic::getFieldDerivative(
    const Vector_t<double, 3>& /*R*/, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/,
    const DiffDirection& /*dir*/
) const {
    throw GeneralClassicException("FM2DMagnetoStatic::getFieldDerivative", "not implemented");
}

void FM2DMagnetoStatic::getFieldDimensions(double& zBegin, double& zEnd) const {
    zBegin = zbegin_m;
    zEnd   = zend_m;
}

void FM2DMagnetoStatic::getFieldDimensions(
    double& /*xIni*/, double& /*xFinal*/, double& /*yIni*/, double& /*yFinal*/, double& /*zIni*/,
    double& /*zFinal*/
) const {
    throw GeneralClassicException("FM2DMagnetoStatic::getFieldDimensions", "not implemented");
}

void FM2DMagnetoStatic::swap() {
    if (swap_m)
        swap_m = false;
    else
        swap_m = true;
}

void FM2DMagnetoStatic::getInfo(Inform* msg) {
    (*msg) << Filename_m << " (2D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m
           << " m;" << endl;
}

double FM2DMagnetoStatic::getFrequency() const {
    throw GeneralClassicException("FM2DMagnetoStatic::getFrequency", "not implemented");
    return 0.0;
}

void FM2DMagnetoStatic::setFrequency(double /*freq*/) {
    throw GeneralClassicException("FM2DMagnetoStatic::setFrequency", "not implemented");
    return;
}
