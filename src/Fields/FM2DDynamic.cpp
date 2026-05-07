#include "Fields/FM2DDynamic.h"
#include "Fields/Fieldmap.hpp"
#include "PartBunch/PartBunch.h"
#include "Physics/Units.h"
#include "Utilities/GeneralOpalException.h"
#include "Utilities/Util.h"

#include <cmath>
#include <fstream>
#include <ios>

FM2DDynamic::FM2DDynamic(const std::string& filename) : Fieldmap(filename) {
    std::ifstream file;
    std::string tmpString;
    double tmpDouble;

    Type = T2DDynamic;  // I can put ut also as Type(T2DDynamic) in the header file.

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if (file.good()) {
        bool parsing_passed = true;
        try {
            parsing_passed = interpretLine<std::string, std::string>(file, tmpString, tmpString);
        } catch (GeneralOpalException& e) {
            parsing_passed = interpretLine<std::string, std::string, std::string>(
                    file, tmpString, tmpString, tmpString);

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" && tmpString != "FALSE")
                throw GeneralOpalException(
                        "FM2DDynamic::FM2DDynamic",
                        "The third string on the first line of 2D field "
                        "maps has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }

        if (tmpString == "ZX") {
            swap_m = true;
            /// Parse rbegin_m, rend_m and num_gridpr_m(-1)
            parsing_passed =
                    parsing_passed
                    && interpretLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
            /// Parse frequency_m
            parsing_passed = parsing_passed && interpretLine<double>(file, frequency_m);
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
            /// Parse frequency_m
            parsing_passed = parsing_passed && interpretLine<double>(file, frequency_m);
            /// Parse rbegin_m, rend_m and num_gridpr_m(-1)
            parsing_passed =
                    parsing_passed
                    && interpretLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        } else {
            std::cerr << "unknown orientation of 2D dynamic fieldmap" << std::endl;
            parsing_passed = false;
            zbegin_m       = 0.0;
            zend_m         = -1e-3;
        }

        for (long i = 0; (i < (num_gridpz_m + 1) * (num_gridpr_m + 1)) && parsing_passed; ++i) {
            parsing_passed = parsing_passed
                             && interpretLine<double, double, double, double>(
                                     file, tmpDouble, tmpDouble, tmpDouble, tmpDouble);
        }

        parsing_passed = parsing_passed && interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if (!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralOpalException(
                    "FM2DDynamic::FM2DDynamic",
                    "An error occured when reading the fieldmap '" + Filename_m + "'");
        } else {
            // convert MHz to Hz and frequency to angular frequency
            frequency_m *= Physics::two_pi * Units::MHz2Hz;

            // convert cm to m
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

FM2DDynamic::~FM2DDynamic() { freeMap(); }

void FM2DDynamic::readMap() {
    if (FieldstrengthEz_m.extent(0) == 0) {
        // declare variables and allocate memory
        std::ifstream in;
        std::string tmpString;
        double tmpDouble, Ezmax = 0.0;

        const size_t size = num_gridpz_m * num_gridpr_m;
        FieldstrengthEz_m = Kokkos::DualView<double*>("FieldstrengthEz", size);
        FieldstrengthEr_m = Kokkos::DualView<double*>("FieldstrengthEr", size);
        FieldstrengthBt_m = Kokkos::DualView<double*>("FieldstrengthBt", size);

        auto Ez = FieldstrengthEz_m.view_host();
        auto Er = FieldstrengthEr_m.view_host();
        auto Bt = FieldstrengthBt_m.view_host();

        // read in field map and parse it
        in.open(Filename_m.c_str());
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);

        if (swap_m) {
            for (int i = 0; i < num_gridpz_m; i++) {
                for (int j = 0; j < num_gridpr_m; j++) {
                    interpretLine<double, double, double, double>(
                            in, Er(j * num_gridpz_m + i), Ez(j * num_gridpz_m + i),
                            Bt(j * num_gridpz_m + i), tmpDouble);
                }
            }
        } else {
            for (int j = 0; j < num_gridpr_m; j++) {
                for (int i = 0; i < num_gridpz_m; i++) {
                    interpretLine<double, double, double, double>(
                            in, Ez(j * num_gridpz_m + i), Er(j * num_gridpz_m + i), tmpDouble,
                            Bt(j * num_gridpz_m + i));
                }
            }
        }
        in.close();

        if (normalize_m) {
            // find maximum field
            // Is (i < num_gridpz_m) enough or should it be (i < size)?
            for (size_t i = 0; i < size; ++i) {
                if (std::abs(Ez(i)) > Ezmax) {
                    Ezmax = std::abs(Ez(i));
                }
            }
        } else {
            Ezmax = 1.0;
        }

        // Precompute scaling factor once
        double const scaleE = 1.e6 / Ezmax;           // MV/m -> V/m conversion
        double const scaleB = Physics::mu_0 / Ezmax;  // H -> B conversion

        for (size_t i = 0; i < size; i++) {
            Ez(i) *= scaleE;
            Er(i) *= scaleE;
            Bt(i) *= scaleB;
        }

        FieldstrengthEz_m.modify<Kokkos::HostSpace>();
        FieldstrengthEz_m.sync<Kokkos::DefaultExecutionSpace>();

        FieldstrengthEr_m.modify<Kokkos::HostSpace>();
        FieldstrengthEr_m.sync<Kokkos::DefaultExecutionSpace>();

        FieldstrengthBt_m.modify<Kokkos::HostSpace>();
        FieldstrengthBt_m.sync<Kokkos::DefaultExecutionSpace>();

        Inform m("FM2DDynamic::readMap");
        m << level3 << "Read in fieldmap '" << Filename_m << "'" << endl;
    }
}

void FM2DDynamic::freeMap() {
    if (FieldstrengthEz_m.extent(0) != 0) {
        FieldstrengthEz_m = Kokkos::DualView<double*>();
        FieldstrengthEr_m = Kokkos::DualView<double*>();
        FieldstrengthBt_m = Kokkos::DualView<double*>();

        Inform m("FM2DDynamic::freeMap");
        m << level3 << "Freed fieldmap '" << Filename_m << "'" << endl;
    }
}

/**
 * @brief Apply the FM to all the particles.
 *
 * @param pc Particle container
 * param scale Scaling factor for the field (currently not used)
 */
void FM2DDynamic::applyField(std::shared_ptr<ParticleContainer_t> pc, double) {
    // Local copies of member variables for use in the lambda function
    double zbegin  = zbegin_m;
    double zend    = zend_m;
    double rend    = rend_m;
    double hr      = hr_m;
    double hz      = hz_m;
    int num_gridpr = num_gridpr_m;
    int num_gridpz = num_gridpz_m;

    // Device accessible views
    auto Ez_device = FieldstrengthEz_m.view_device();
    auto Er_device = FieldstrengthEr_m.view_device();
    auto Bt_device = FieldstrengthBt_m.view_device();

    auto Rview = pc->R.getView();
    auto Eview = pc->E.getView();
    auto Bview = pc->B.getView();

    const size_t nLocal = pc->getLocalNum();

    Kokkos::parallel_for(
            "FM2DDynamic::applyField", nLocal, KOKKOS_LAMBDA(const size_t i) {
                if (Rview(i)(2) >= zbegin && Rview(i)(2) < zend
                    && sqrt(Rview(i)(0) * Rview(i)(0) + Rview(i)(1) * Rview(i)(1)) < rend) {
                    computeField(
                            Rview(i), Eview(i), Bview(i), Ez_device, Er_device, Bt_device, hr, hz,
                            zbegin, num_gridpr, num_gridpz);
                }
            });
}

void FM2DDynamic::applyRFField(
        std::shared_ptr<ParticleContainer_t> pc, double electricScale, double magneticScale,
        double startField, double endField) {
    double zbegin  = zbegin_m;
    double zend    = zend_m;
    double rend    = rend_m;
    double hr      = hr_m;
    double hz      = hz_m;
    int num_gridpr = num_gridpr_m;
    int num_gridpz = num_gridpz_m;

    auto Ez_device = FieldstrengthEz_m.view_device();
    auto Er_device = FieldstrengthEr_m.view_device();
    auto Bt_device = FieldstrengthBt_m.view_device();

    auto Rview = pc->R.getView();
    auto Eview = pc->E.getView();
    auto Bview = pc->B.getView();

    const size_t nLocal = pc->getLocalNum();

    Kokkos::parallel_for(
            "FM2DDynamic::applyRFField", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const auto& R = Rview(i);

                if (R(2) >= startField && R(2) < endField && R(2) >= zbegin && R(2) < zend
                    && sqrt(R(0) * R(0) + R(1) * R(1)) < rend) {
                    Vector_t<double, 3> tmpE(0.0), tmpB(0.0);

                    computeField(
                            R, tmpE, tmpB, Ez_device, Er_device, Bt_device, hr, hz, zbegin,
                            num_gridpr, num_gridpz);

                    Eview(i) += electricScale * tmpE;
                    Bview(i) += magneticScale * tmpB;
                }
            });
}

/**
 * @brief Get the fieldstrength at position R
 *
 * @param R Position
 * @param E Electric Field
 * @param B Magnetic Field
 *
 * @return true if R is outside of the field map, false otherwise.
 */
bool FM2DDynamic::getFieldstrength(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& B) const {
    if (isInside(R)) {
        computeField(
                R, E, B, FieldstrengthEz_m.h_view, FieldstrengthEr_m.h_view,
                FieldstrengthBt_m.h_view, hr_m, hz_m, zbegin_m, num_gridpr_m, num_gridpz_m);

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
bool FM2DDynamic::getFieldDerivative(
        const Vector_t<double, 3>& /*R*/, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/,
        const DiffDirection& /*dir*/) const {
    throw GeneralOpalException("FM2DDynamic::getFieldDerivative", "not implemented");
}

void FM2DDynamic::getFieldDimensions(double& zBegin, double& zEnd) const {
    zBegin = zbegin_m;
    zEnd   = zend_m;
}

void FM2DDynamic::getFieldDimensions(
        double& /*xIni*/, double& /*xFinal*/, double& /*yIni*/, double& /*yFinal*/,
        double& /*zIni*/, double& /*zFinal*/) const {
    throw GeneralOpalException("FM2DDynamic::getFieldDimensions", "not implemented");
}

void FM2DDynamic::swap() { swap_m = !swap_m; }

void FM2DDynamic::getInfo(Inform* msg) {
    (*msg) << Filename_m << " (2D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;"
           << endl;
}

double FM2DDynamic::getFrequency() const { return frequency_m; }

void FM2DDynamic::setFrequency(double freq) { frequency_m = freq; }

void FM2DDynamic::getOnaxisEz(std::vector<std::pair<double, double>>& F) {
    double dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    F.resize(num_gridpz_m);

    auto Ez = FieldstrengthEz_m.view_host();  // Access host view

    for (int i = 0; i < num_gridpz_m; ++i) {
        F[i].first  = dz * i;
        F[i].second = Ez(i) / 1e6;  // Convert V/m -> MV/m
    }
}
