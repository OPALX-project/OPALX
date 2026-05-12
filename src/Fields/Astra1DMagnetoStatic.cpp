#include "Fields/Astra1DMagnetoStatic.h"
#include "Fields/Fieldmap.hpp"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Utilities/GSLSpline.h"
#include "Utilities/GeneralOpalException.h"
#include "Utilities/Util.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <ios>
#include <sstream>

/**
 * @brief Constructor for 1D magnetostatic field map.
 *
 * Parses the fieldmap header and scans the longitudinal on-axis Bz data.
 * The constructor determines:
 * - fieldmap type (`TAstraMagnetoStatic`)
 * - Fourier interpolation accuracy (`accuracy_m`)
 * - optional normalization flag (`normalize_m`)
 * - longitudinal field boundaries (`zbegin_m`, `zend_m`)
 * - number of valid strictly increasing z samples (`num_gridpz_m`)
 * - effective mirrored-periodic interpolation length (`length_m`)
 *
 * The actual Fourier coefficients are generated later in `readMap()`.
 */
Astra1DMagnetoStatic::Astra1DMagnetoStatic(const std::string& filename) : Fieldmap(filename) {
    std::ifstream file;
    std::string tmpString;
    double tmpDouble;
    // int skippedValues = 0;

    Type = TAstraMagnetoStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if (file.good()) {
        bool parsing_passed = true;
        try {
            parsing_passed = interpretLine<std::string, int>(file, tmpString, accuracy_m);
        } catch (GeneralOpalException& e) {
            parsing_passed = interpretLine<std::string, int, std::string>(
                    file, tmpString, accuracy_m, tmpString);

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" && tmpString != "FALSE")
                throw GeneralOpalException(
                        "Astra1DMagnetoStatic::Astra1DMagnetoStatic",
                        "The third string on the first line of 1D field maps "
                        "has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }

        int acceptedValues = 0;

        parsing_passed = parsing_passed && interpretLine<double, double>(file, zbegin_m, tmpDouble);

        double lastAcceptedZ = zbegin_m;

        if (parsing_passed) {
            ++acceptedValues;
            zend_m = zbegin_m;
        }

        while (true) {
            const bool line_ok = interpretLine<double, double>(file, zend_m, tmpDouble, false);
            if (!line_ok) {
                break;
            }

            if (zend_m - lastAcceptedZ > 1e-10) {
                lastAcceptedZ = zend_m;
                ++acceptedValues;
            }
        }

        num_gridpz_m = acceptedValues;
        lines_read_m = 0;

        if (!parsing_passed && !file.eof()) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralOpalException(
                    "Astra1DMagnetoStatic::Astra1DMagnetoStatic",
                    "Error reading fieldmap '" + Filename_m + "'");
        }

        if (num_gridpz_m < 2) {
            throw GeneralOpalException(
                    "Astra1DMagnetoStatic::Astra1DMagnetoStatic",
                    "Fieldmap must contain at least two valid sampling points");
        }
        length_m = 2.0 * num_gridpz_m * (zend_m - zbegin_m) / (num_gridpz_m - 1);

        file.close();
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m   = -1e-3;
    }
}

Astra1DMagnetoStatic::~Astra1DMagnetoStatic() { freeMap(); }

void Astra1DMagnetoStatic::readMap() {
    if (FourCoefs_m.extent(0) != 0) {
        return;
    }

    if (num_gridpz_m < 2) {
        throw GeneralOpalException(
                "Astra1DMagnetoStatic::readMap",
                "Fieldmap must contain at least two valid sampling points");
    }

    std::ifstream in(Filename_m.c_str());
    if (!in.good()) {
        throw GeneralOpalException(
                "Astra1DMagnetoStatic::readMap", "Cannot open fieldmap '" + Filename_m + "'");
    }

    const double dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);

    // CPU arrays - GSL (FFT + spline) is CPU-only
    double* RealValues = new double[2 * num_gridpz_m];
    double* zvals      = new double[num_gridpz_m];

    for (int i = 0; i < 2 * num_gridpz_m; ++i) {
        RealValues[i] = 0.0;
    }

    for (int i = 0; i < num_gridpz_m; ++i) {
        zvals[i] = 0.0;
    }

    gsl_spline* Bz_interpolant = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
    gsl_interp_accel* Bz_accel = gsl_interp_accel_alloc();

    const int size = 2 * accuracy_m - 1;
    FourCoefs_m    = Kokkos::DualView<double*>("FourCoefs", size);
    auto coefs     = FourCoefs_m.view_host();

    std::string tmpString;
    getLine(in, tmpString);

    double Bz_max        = 0.0;
    double lastAcceptedZ = zbegin_m - dz;
    int accepted         = 0;

    while (accepted < num_gridpz_m) {
        double ztmp  = 0.0;
        double bztmp = 0.0;

        const bool ok = interpretLine<double, double>(in, ztmp, bztmp, false);
        if (!ok) {
            break;
        }

        if (ztmp - lastAcceptedZ > 1e-10) {
            zvals[accepted]      = ztmp;
            RealValues[accepted] = bztmp;
            Bz_max               = std::max(Bz_max, std::abs(bztmp));
            lastAcceptedZ        = ztmp;
            ++accepted;
        }
    }
    in.close();

    if (accepted != num_gridpz_m) {
        gsl_spline_free(Bz_interpolant);
        gsl_interp_accel_free(Bz_accel);
        delete[] zvals;
        delete[] RealValues;

        std::ostringstream os;
        os << "Mismatch between counted and parsed fieldmap points in '" << Filename_m
           << "': expected " << num_gridpz_m << ", got " << accepted;
        throw GeneralOpalException("Astra1DMagnetoStatic::readMap", os.str());
    }

    if (Bz_max == 0.0) {
        gsl_spline_free(Bz_interpolant);
        gsl_interp_accel_free(Bz_accel);
        delete[] zvals;
        delete[] RealValues;

        throw GeneralOpalException(
                "Astra1DMagnetoStatic::readMap",
                "Maximum on-axis magnetic field is zero in fieldmap '" + Filename_m + "'");
    }

    gsl_spline_init(Bz_interpolant, zvals, RealValues, num_gridpz_m);

    int ii = num_gridpz_m;
    for (int i = 0; i < num_gridpz_m - 1; ++i, ++ii) {
        const double z = zbegin_m + dz * i;
        RealValues[ii] = gsl_spline_eval(Bz_interpolant, z, Bz_accel);
    }

    // Ensure periodicity by mirroring
    RealValues[ii++] = RealValues[num_gridpz_m - 1];
    --ii;
    for (int i = 0; i < num_gridpz_m; ++i, --ii) {
        RealValues[i] = RealValues[ii];
    }

    // Disable normalization if requested
    const double norm = normalize_m ? Bz_max : 1.0;

    // Compute Fourier coefficients explicitly from the mirrored periodic samples.
    const int M = 2 * num_gridpz_m;

    double a0 = 0.0;
    for (int j = 0; j < M; ++j) {
        a0 += RealValues[j];
    }
    coefs(0) = a0 / (norm * M);

    for (int l = 1; l < accuracy_m; ++l) {
        double a_l = 0.0;
        double b_l = 0.0;

        for (int j = 0; j < M; ++j) {
            const double theta = Physics::two_pi * double(j) / double(M);
            a_l += RealValues[j] * std::cos(l * theta);
            b_l += RealValues[j] * std::sin(l * theta);
        }

        a_l *= 2.0 / double(M);
        b_l *= 2.0 / double(M);

        const int n = 2 * l - 1;

        coefs(n)     = a_l / norm;
        coefs(n + 1) = -b_l / norm;
    }

    gsl_spline_free(Bz_interpolant);
    gsl_interp_accel_free(Bz_accel);

    delete[] zvals;
    delete[] RealValues;

    FourCoefs_m.modify<Kokkos::HostSpace>();
    FourCoefs_m.sync<Kokkos::DefaultExecutionSpace>();

    Inform m("Astra1DMagnetoStatic::readMap");
    m << level3 << "Read in fieldmap '" << Filename_m << "'" << endl;
}

void Astra1DMagnetoStatic::freeMap() {
    if (FourCoefs_m.extent(0) != 0) {
        FourCoefs_m = Kokkos::DualView<double*>();

        Inform m("Astra1DMagnetoStatic::freeMap");
        m << level3 << "Freed fieldmap '" << Filename_m << "'" << endl;
    }
}

/**
 * @brief Apply the FM to all the particles.
 *
 * @param pc Particle container
 */
void Astra1DMagnetoStatic::applyField(std::shared_ptr<ParticleContainer_t> pc, double scale) {
    const double zbegin = zbegin_m;
    const double zend   = zend_m;
    const double length = length_m;
    const int accuracy  = accuracy_m;
    const size_t nLocal = pc->getLocalNum();

    auto FourCoefs_device = FourCoefs_m.view_device();

    auto Rview = pc->R.getView();
    auto Eview = pc->E.getView();
    auto Bview = pc->B.getView();

    Kokkos::parallel_for(
            "Astra1DMagnetoStatic::applyField", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const auto& R = Rview(i);

                if (R(2) >= zbegin && R(2) < zend) {
                    Vector_t<double, 3> localE(0.0);
                    Vector_t<double, 3> localB(0.0);

                    computeField(R, localE, localB, FourCoefs_device, zbegin, length, accuracy);

                    Eview(i) += scale * localE;
                    Bview(i) += scale * localB;
                }
            });
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
bool Astra1DMagnetoStatic::getFieldstrength(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& B) const {
    if (isInside(R)) {
        computeField(R, E, B, FourCoefs_m.view_host(), zbegin_m, length_m, accuracy_m);

        return false;
    }

    return true;
}

/**
 * @brief Get the derivative of the field at position R
 *
 * @param R Position
 * @param E Electric Field
 * @param B Magnetic field
 *
 * @note Not implemented yet
 *
 */
bool Astra1DMagnetoStatic::getFieldDerivative(
        const Vector_t<double, 3>& /*R*/, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/,
        const DiffDirection& /*dir*/) const {
    return false;
}

void Astra1DMagnetoStatic::getFieldDimensions(double& zBegin, double& zEnd) const {
    zBegin = zbegin_m;
    zEnd   = zend_m;
}

void Astra1DMagnetoStatic::getFieldDimensions(
        double& /*xIni*/, double& /*xFinal*/, double& /*yIni*/, double& /*yFinal*/,
        double& /*zIni*/, double& /*zFinal*/) const {
    throw GeneralOpalException("Astra1DMagnetoStatic::getFieldDimensions", "not implemented");
}

void Astra1DMagnetoStatic::swap() {}

void Astra1DMagnetoStatic::getInfo(Inform* msg) {
    (*msg) << Filename_m << " (1D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m
           << " m;" << endl;
}

double Astra1DMagnetoStatic::getFrequency() const { return 0.0; }

void Astra1DMagnetoStatic::setFrequency(double /*freq*/) {}
