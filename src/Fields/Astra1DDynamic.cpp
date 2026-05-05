#include "Fields/Astra1DDynamic.h"
#include "Fields/Fieldmap.hpp"
#include "PartBunch/PartBunch.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/GSLFFT.h"
#include "Utilities/GSLSpline.h"
#include "Utilities/GeneralOpalException.h"
#include "Utilities/Util.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <ios>

Astra1DDynamic::Astra1DDynamic(const std::string& filename) : Fieldmap(filename) {
    std::ifstream file;
    std::string tmpString;
    // int skippedValues = 0;
    double tmpDouble;

    Type = TAstraDynamic;

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
                        "Astra1DDynamic::Astra1DDynamic",
                        "The third string on the first line of 1D field "
                        "maps has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }

        parsing_passed = parsing_passed && interpretLine<double>(file, frequency_m);

        int acceptedValues = 0;

        // Read first data point
        parsing_passed = parsing_passed && interpretLine<double, double>(file, zbegin_m, tmpDouble);

        double tmpDouble2 = zbegin_m;

        if (parsing_passed) {
            ++acceptedValues;
            zend_m = zbegin_m;
        }

        // Read remaining data points
        while (true) {
            bool line_ok = interpretLine<double, double>(file, zend_m, tmpDouble, false);
            if (!line_ok) {
                break;
            }

            if (zend_m - tmpDouble2 > 1e-10) {
                tmpDouble2 = zend_m;
                ++acceptedValues;
            }
        }

        num_gridpz_m = acceptedValues;
        lines_read_m = 0;

        if (!parsing_passed && !file.eof()) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralOpalException(
                    "Astra1DDynamic::Astra1DDynamic",
                    "Error reading fieldmap '" + Filename_m + "'");
        } else {
            // conversion from MHz to Hz and from frequency to angular frequency
            frequency_m *= Physics::two_pi * Units::MHz2Hz;
            xlrep_m = frequency_m / Physics::c;
        }
        length_m = 2.0 * num_gridpz_m * (zend_m - zbegin_m) / (num_gridpz_m - 1);

        file.close();
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m   = -1e-3;
    }
}

Astra1DDynamic::~Astra1DDynamic() { freeMap(); }

void Astra1DDynamic::readMap() {
    if (FourCoefs_m.extent(0) != 0) {
        return;
    }

    // Need at least two valid points for dz, spline, and FFT setup
    if (num_gridpz_m < 2) {
        throw GeneralOpalException(
                "Astra1DDynamic::readMap",
                "Fieldmap must contain at least two valid sampling points");
    }

    std::ifstream in(Filename_m.c_str());
    if (!in.good()) {
        throw GeneralOpalException(
                "Astra1DDynamic::readMap", "Cannot open fieldmap '" + Filename_m + "'");
    }

    const double dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);

    // CPU arrays - GSL (FFT + spline) is CPU-only
    double* RealValues = new double[2 * num_gridpz_m];
    double* zvals      = new double[num_gridpz_m];

    // Initialize to something known for easier debugging
    for (int i = 0; i < 2 * num_gridpz_m; ++i) {
        RealValues[i] = 0.0;
    }
    for (int i = 0; i < num_gridpz_m; ++i) {
        zvals[i] = 0.0;
    }

    // GSL objects (CPU only)
    gsl_spline* Ez_interpolant = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
    gsl_interp_accel* Ez_accel = gsl_interp_accel_alloc();

    // Future optimization candidate:
    // gsl_fft_real_wavetable* real =
    //     gsl_fft_real_wavetable_alloc(2 * num_gridpz_m);
    // gsl_fft_real_workspace* work =
    //     gsl_fft_real_workspace_alloc(2 * num_gridpz_m);

    // Allocate Kokkos View (final data used on GPU)
    const int size = 2 * accuracy_m - 1;
    FourCoefs_m    = Kokkos::DualView<double*>("FourCoefs", size);
    auto coefs     = FourCoefs_m.view_host();

    // ---- read file / parse fieldmap ----
    // Skip the first two header lines exactly like the original code
    std::string tmpString;
    getLine(in, tmpString);
    getLine(in, tmpString);

    double Ez_max        = 0.0;
    double lastAcceptedZ = zbegin_m - dz;
    int accepted         = 0;

    while (accepted < num_gridpz_m) {
        double ztmp  = 0.0;
        double eztmp = 0.0;

        const bool ok = interpretLine<double, double>(in, ztmp, eztmp, false);
        if (!ok) {
            break;
        }

        if (ztmp - lastAcceptedZ > 1e-10) {
            zvals[accepted]      = ztmp;
            RealValues[accepted] = eztmp;
            Ez_max               = std::max(Ez_max, std::abs(eztmp));
            lastAcceptedZ        = ztmp;
            ++accepted;
        }
    }
    in.close();

    if (accepted != num_gridpz_m) {
        gsl_spline_free(Ez_interpolant);
        gsl_interp_accel_free(Ez_accel);
        delete[] zvals;
        delete[] RealValues;

        std::ostringstream os;
        os << "Mismatch between counted and parsed fieldmap points in '" << Filename_m
           << "': expected " << num_gridpz_m << ", got " << accepted;
        throw GeneralOpalException("Astra1DDynamic::readMap", os.str());
    }

    if (Ez_max == 0.0) {
        gsl_spline_free(Ez_interpolant);
        gsl_interp_accel_free(Ez_accel);
        delete[] zvals;
        delete[] RealValues;

        throw GeneralOpalException(
                "Astra1DDynamic::readMap",
                "Maximum on-axis field is zero in fieldmap '" + Filename_m + "'");
    }

    zRaw_m.resize(num_gridpz_m);
    ezRaw_m.resize(num_gridpz_m);
    for (int i = 0; i < num_gridpz_m; ++i) {
        zRaw_m[i]  = zvals[i];
        ezRaw_m[i] = RealValues[i];
    }

    // Interpolate onto an equidistant grid and build a mirrored periodic array
    gsl_spline_init(Ez_interpolant, zvals, RealValues, num_gridpz_m);

    int ii = num_gridpz_m;
    for (int i = 0; i < num_gridpz_m - 1; ++i, ++ii) {
        const double z = zbegin_m + dz * i;
        RealValues[ii] = gsl_spline_eval(Ez_interpolant, z, Ez_accel);
    }

    // Ensure periodicity by mirroring
    RealValues[ii++] = RealValues[num_gridpz_m - 1];
    --ii;
    for (int i = 0; i < num_gridpz_m; ++i, --ii) {
        RealValues[i] = RealValues[ii];
    }

    // Disable normalization if requested
    const double norm = normalize_m ? Ez_max : 1.0;

    // Compute Fourier coefficients explicitly from the mirrored periodic samples.
    const int M = 2 * num_gridpz_m;

    double a0 = 0.0;
    for (int j = 0; j < M; ++j) {
        a0 += RealValues[j];
    }
    coefs(0) = a0 / (norm * Units::Vpm2MVpm * M);

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

        const int n  = 2 * l - 1;
        coefs(n)     = a_l / (norm * Units::Vpm2MVpm);
        coefs(n + 1) = -b_l / (norm * Units::Vpm2MVpm);
    }

    // Future optimization candidate:
    // gsl_fft_real_transform(RealValues, 1, 2 * num_gridpz_m, real, work);
    // coefs(0) = RealValues[0] / (norm * Units::Vpm2MVpm * 2.0 * num_gridpz_m);
    // for (int i = 1; i < size; ++i) {
    //     coefs(i) = RealValues[i] / (norm * Units::Vpm2MVpm * num_gridpz_m);
    // }

    gsl_spline_free(Ez_interpolant);
    gsl_interp_accel_free(Ez_accel);
    // gsl_fft_real_workspace_free(work);
    // gsl_fft_real_wavetable_free(real);

    delete[] zvals;
    delete[] RealValues;

    FourCoefs_m.modify<Kokkos::HostSpace>();
    FourCoefs_m.sync<Kokkos::DefaultExecutionSpace>();

    Inform m("Astra1DDynamic::readMap");
    m << level3 << "Read in fieldmap '" << Filename_m << "'" << endl;
}

void Astra1DDynamic::freeMap() {
    if (FourCoefs_m.extent(0) != 0) {
        // Reset View (releases memory)
        FourCoefs_m = Kokkos::DualView<double*>();

        Inform m("Astra1DDynamic::freeMap");
        m << level3 << "Freed fieldmap '" << Filename_m << "'" << endl;
    }
}

/**
 * @brief Apply the FM to all the particles.
 *
 * @param pc Particle container
 * @param scale Scaling factor for the field (currently not used)
 */
void Astra1DDynamic::applyField(std::shared_ptr<ParticleContainer_t> pc, double) {
    // Local copies of member variables for use in the lambda function
    const double zbegin = zbegin_m;
    const double zend   = zend_m;
    const double length = length_m;
    const double xlrep  = xlrep_m;
    const int accuracy  = accuracy_m;
    const size_t nLocal = pc->getLocalNum();

    // Device-accessible Fourier coefficients
    auto FourCoefs_device = FourCoefs_m.view_device();

    auto Rview = pc->R.getView();
    auto Eview = pc->E.getView();
    auto Bview = pc->B.getView();

    Kokkos::parallel_for(
            "Astra1DDynamic::applyField", nLocal, KOKKOS_LAMBDA(const size_t i) {
                const auto& R = Rview(i);

                // 1D Astra field map is only bounded in z
                if (R(2) >= zbegin && R(2) < zend) {
                    computeField(
                            R, Eview(i), Bview(i), FourCoefs_device, zbegin, length, xlrep,
                            accuracy);
                }
            });
}

void Astra1DDynamic::applyTravelingWave(
        std::shared_ptr<ParticleContainer_t> pc, double entryElectricScale,
        double entryMagneticScale, double core1ElectricScale, double core1MagneticScale,
        double core2ElectricScale, double core2MagneticScale, double exitElectricScale,
        double exitMagneticScale, double /*startField*/, double startCoreField,
        double startExitField, double mappedStartExitField, double periodLength, double cellLength,
        double elementLength) {
    const double zbegin = zbegin_m;
    const double zend   = zend_m;
    const double length = length_m;
    const double xlrep  = xlrep_m;
    const int accuracy  = accuracy_m;
    const size_t nLocal = pc->getLocalNum();

    auto FourCoefs_device = FourCoefs_m.view_device();

    auto Rview = pc->R.getView();
    auto Eview = pc->E.getView();
    auto Bview = pc->B.getView();

    Kokkos::parallel_for(
            "Astra1DDynamic::applyTravelingWave", nLocal, KOKKOS_LAMBDA(const size_t i) {
                computeTravelingWaveField(
                        Rview(i), Eview(i), Bview(i), FourCoefs_device, zbegin, zend, length, xlrep,
                        accuracy, entryElectricScale, entryMagneticScale, core1ElectricScale,
                        core1MagneticScale, core2ElectricScale, core2MagneticScale,
                        exitElectricScale, exitMagneticScale, startCoreField, startExitField,
                        mappedStartExitField, periodLength, cellLength, elementLength);
            });
}

bool Astra1DDynamic::getFieldstrength(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& B) const {
    if (isInside(R)) {
        computeField(R, E, B, FourCoefs_m.h_view, zbegin_m, length_m, xlrep_m, accuracy_m);

        return false;
    } else {
        return true;
    }
}

/**
 * @brief Get the derivative of the field at position R
 *
 * @param R Position
 * @param E Electric Field
 * @param B Derivate of the magnetic field (unused)
 *
 */
bool Astra1DDynamic::getFieldDerivative(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& /*B*/,
        const DiffDirection& /*dir*/) const {
    const double kz = Physics::two_pi * (R(2) - zbegin_m) / length_m + Physics::pi;
    double ezp      = 0.0;

    auto FourCoefs = FourCoefs_m.h_view;

    int n = 1;
    for (int l = 1; l < accuracy_m; ++l, n += 2) {
        ezp += Physics::two_pi / length_m * l
               * (-FourCoefs(n) * std::sin(kz * l) - FourCoefs(n + 1) * std::cos(kz * l));
    }

    E(2) += ezp;

    return false;
}

void Astra1DDynamic::getFieldDimensions(double& zBegin, double& zEnd) const {
    zBegin = zbegin_m;
    zEnd   = zend_m;
}

void Astra1DDynamic::getFieldDimensions(
        double& /*xIni*/, double& /*xFinal*/, double& /*yIni*/, double& /*yFinal*/,
        double& /*zIni*/, double& /*zFinal*/) const {
    throw GeneralOpalException("Astra1DDynamic::getFieldDimensions", "not implemented");
}

void Astra1DDynamic::swap() {}

void Astra1DDynamic::getInfo(Inform* msg) {
    (*msg) << Filename_m << " (1D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;"
           << endl;
}

double Astra1DDynamic::getFrequency() const { return frequency_m; }

void Astra1DDynamic::setFrequency(double freq) { frequency_m = freq; }

// This function re-reads the fieldmap file to extract the on-axis Ez values,
// which are stored in F as pairs of (z, Ez).
// Not computationally optimal, but it was like that in the original
// OPAL implementation and the on-axis field values are needed for the
// `TravelingWave` element, which uses the `Astra1DDynamic` fieldmap.
void Astra1DDynamic::getOnaxisEz(std::vector<std::pair<double, double>>& F) {
    double Ez_max = 0.0;
    double tmpDouble;
    int tmpInt;
    std::string tmpString;

    F.resize(num_gridpz_m);

    std::ifstream in(Filename_m.c_str());
    interpretLine<std::string, int>(in, tmpString, tmpInt);
    interpretLine<double>(in, tmpDouble);

    for (int i = 0; i < num_gridpz_m; ++i) {
        interpretLine<double, double>(in, F[i].first, F[i].second);
        if (std::abs(F[i].second) > Ez_max) {
            Ez_max = std::abs(F[i].second);
        }
    }
    in.close();

    for (int i = 0; i < num_gridpz_m; ++i) {
        F[i].second /= Ez_max;
    }
}
