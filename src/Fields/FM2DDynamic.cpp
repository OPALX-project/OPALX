#include "Fields/FM2DDynamic.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <fstream>
#include <ios>
#include <cmath>


FM2DDynamic::FM2DDynamic(const std::string& filename)
    : Fieldmap(filename),
      FieldstrengthEz_m(nullptr),
      FieldstrengthEr_m(nullptr),
      FieldstrengthBt_m(nullptr) {
    std::ifstream file;
    std::string tmpString;
    double tmpDouble;

    Type = T2DDynamic; // I can put ut also as Type(T2DDynamic) in the header file.

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = true;
        try {
            parsing_passed = 
            interpretLine<std::string, std::string>(file, tmpString, tmpString);
        } catch (GeneralClassicException &e) {
            parsing_passed = 
            interpretLine<std::string, std::string, std::string>(
                file, tmpString, tmpString, tmpString);

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" && tmpString != "FALSE")
                throw GeneralClassicException(
                    "FM2DDynamic::FM2DDynamic", 
                    "The third string on the first line of 2D field "
                    "maps has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }

        if(tmpString == "ZX") {
            swap_m = true;
            /// Parse rbegin_m, rend_m and num_gridpr_m(-1)
            parsing_passed = 
                parsing_passed
                && interpretLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
            /// Parse frequency_m
            parsing_passed = 
                parsing_passed
                && interpretLine<double>(file, frequency_m);
            /// Parse rbegin_m, zend_m and num_gridpz_m(-1)
            parsing_passed = 
                parsing_passed
                && interpretLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        } else if(tmpString == "XZ") {
            swap_m = false;
            /// Parse rbegin_m, zend_m and num_gridpz_m(-1)
            parsing_passed = 
                parsing_passed
                && interpretLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
            /// Parse frequency_m
            parsing_passed = 
                parsing_passed
                && interpretLine<double>(file, frequency_m);
            /// Parse rbegin_m, rend_m and num_gridpr_m(-1)
            parsing_passed = 
                parsing_passed
                && interpretLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        } else {
            std::cerr << "unknown orientation of 2D dynamic fieldmap" << std::endl;
            parsing_passed = false;
            zbegin_m = 0.0;
            zend_m = -1e-3;
        }

        for(long i = 0; (i < (num_gridpz_m + 1) * (num_gridpr_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = 
                parsing_passed 
                && interpretLine<double, double, double, double>(
                    file, tmpDouble, tmpDouble, tmpDouble, tmpDouble);
        }

        parsing_passed = parsing_passed && interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralClassicException(
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
        zend_m = -1e-3;
    }
}


FM2DDynamic::~FM2DDynamic() {
    freeMap();
}

// FM2DDynamic _FM2DDynamic::create(const std::string& filename)
// {
//     return FM2DDynamic(new _FM2DDynamic(filename));
// }

void FM2DDynamic::readMap() {
    if(FieldstrengthEz_m.h_view.data() == nullptr) {
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

        // FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
        // FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];
        // FieldstrengthBt_m = new double[num_gridpz_m * num_gridpr_m];

        // read in field map and parse it
        in.open(Filename_m.c_str());
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);

        if(swap_m) {
            for(int i = 0; i < num_gridpz_m; i++) {
                for(int j = 0; j < num_gridpr_m; j++) {
                    interpretLine<double, double, double, double>(
                        in,
                        Er(i + j * num_gridpz_m),
                        Ez(i + j * num_gridpz_m),
                        Bt(i + j * num_gridpz_m),
                        tmpDouble);
                }
            }
        } else {
            for(int j = 0; j < num_gridpr_m; j++) {
                for(int i = 0; i < num_gridpz_m; i++) {
                    interpretLine<double, double, double, double>(
                        in,
                        Ez(i + j * num_gridpz_m),
                        Er(i + j * num_gridpz_m),
                        tmpDouble,
                        Bt(i + j * num_gridpz_m)
                    );

                }
            }
        }
        in.close();

        if (normalize_m) {
            // find maximum field
            // To make it more general should replace num_gridpz_m with size?
            for(int i = 0; i < num_gridpz_m; ++ i) {
                if(std::abs(Ez(i)) > Ezmax) {
                    Ezmax = std::abs(Ez(i));
                }
            }
        } else {
            Ezmax = 1.0;
        }

        // Precompute scaling factor once
        double const scaleE = 1.e6 / Ezmax;          // MV/m -> V/m conversion
        double const scaleB = Physics::mu_0 / Ezmax; // H -> B conversion

        for(int i = 0; i < size; i++) {
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

        *ippl::Info << level3 
            << typeset_msg("read in fieldmap '" + Filename_m + "'", "info")
            << endl;
        // INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
        //         << endl);
    }
}

void FM2DDynamic::freeMap() {
    // if(FieldstrengthEz_m != nullptr) {
    //     delete[] FieldstrengthEz_m;
    //     FieldstrengthEz_m = nullptr;
    //     delete[] FieldstrengthEr_m;
    //     FieldstrengthEr_m = nullptr;
    //     delete[] FieldstrengthBt_m;
    //     FieldstrengthBt_m = nullptr;
    if(FieldstrengthEz_m.h_view.data() != nullptr) {

        FieldstrengthEz_m = Kokkos::DualView<double*>();
        FieldstrengthEr_m = Kokkos::DualView<double*>();
        FieldstrengthBt_m = Kokkos::DualView<double*>();

        *ippl::Info << level3 
            << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") 
            << endl;
    }
}

/**
 * @brief Apply the FM to all the particles.
 * 
 * @param pc Particle container
 */
void FM2DDynamic::applyField(std::shared_ptr<ParticleContainer_t> pc)
{
    // Local copies of member variables for use in the lambda function
    double zbegin = zbegin_m;
    double zend = zend_m;
    double rend = rend_m;
    double hr = hr_m;
    double hz = hz_m;
    int num_gridpr = num_gridpr_m;
    int num_gridpz = num_gridpz_m;

    // Device accessible views
    auto Ez_device = FieldstrengthEz_m.view_device();
    auto Er_device = FieldstrengthEr_m.view_device();
    auto Bt_device = FieldstrengthBt_m.view_device();

    auto Rview = pc->R.getView();
    auto Eview = pc->E.getView();
    auto Bview = pc->B.getView();

    Kokkos::parallel_for("FM2DDynamic::applyField",
    ippl::getRangePolicy(Rview),
    KOKKOS_LAMBDA(const int i)
    {
        if(Rview(i)(2) >= zbegin && Rview(i)(2) < zend &&
            sqrt(Rview(i)(0)*Rview(i)(0) + Rview(i)(1)*Rview(i)(1)) < rend) 
        {
            computeField(Rview(i),
                            Eview(i),
                            Bview(i),
                            Ez_device,
                            Er_device,
                            Bt_device,
                            hr, hz, zbegin,
                            num_gridpr, num_gridpz);
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
    const Vector_t<double,3>& R,
    Vector_t<double,3>& E,
    Vector_t<double,3>& B) const
{
    if (isInside(R)) {

        computeField(
            R,
            E,
            B,
            FieldstrengthEz_m.h_view,
            FieldstrengthEr_m.h_view,
            FieldstrengthBt_m.h_view,
            hr_m,
            hz_m,
            zbegin_m,
            num_gridpr_m,
            num_gridpz_m
        );

        return false;
    }
    else{
        return true;
    }
}

// bool _FM2DDynamic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
//     // do bi-linear interpolation
//     const double RR = std::sqrt(R(0) * R(0) + R(1) * R(1));

//     const int indexr = (int)std::floor(RR / hr_m);
//     const double leverr = RR / hr_m - indexr;

//     const int indexz = (int)std::floor((R(2) - zbegin_m) / hz_m);
//     const double leverz = (R(2) - zbegin_m) / hz_m - indexz;

//     if((indexz < 0) || (indexz + 2 > num_gridpz_m))
//         return false;
//     if(indexr + 2 > num_gridpr_m)
//         return true;

//     const int index1 = indexz + indexr * num_gridpz_m;
//     const int index2 = index1 + num_gridpz_m;

//     double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1]
//                      + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index1 + 1]
//                      + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index2]
//                      + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];

//     double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1]
//                      + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index1 + 1]
//                      + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index2]
//                      + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

//     double BfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBt_m[index1]
//                      + leverz         * (1.0 - leverr) * FieldstrengthBt_m[index1 + 1]
//                      + (1.0 - leverz) * leverr         * FieldstrengthBt_m[index2]
//                      + leverz         * leverr         * FieldstrengthBt_m[index2 + 1];

//     if(RR > 1e-10) {
//         E(0) += EfieldR * R(0) / RR;
//         E(1) += EfieldR * R(1) / RR;
//         B(0) -= BfieldT * R(1) / RR;
//         B(1) += BfieldT * R(0) / RR;
//     }
//     E(2) += EfieldZ;

//     return false;
// }


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
    const Vector_t<double, 3>& /*R*/, 
    Vector_t<double, 3>& /*E*/, 
    Vector_t<double, 3>& /*B*/,
    const DiffDirection &/*dir*/) const {
       throw GeneralClassicException(
            "FM2DDynamic::getFieldDerivative","not implemented");
}

void FM2DDynamic::getFieldDimensions(double& zBegin, double& zEnd) const {
    zBegin = zbegin_m;
    zEnd   = zend_m;
}

void FM2DDynamic::getFieldDimensions(
    double& /*xIni*/, double& /*xFinal*/, double& /*yIni*/, double& /*yFinal*/, 
    double& /*zIni*/, double& /*zFinal*/) const {
        throw GeneralClassicException(
            "FM2DDynamic::getFieldDimensions","not implemented");
    }

void FM2DDynamic::swap() {
    swap_m = !swap_m;
    // if(swap_m) swap_m = false;
    // else swap_m = true;
}

void FM2DDynamic::getInfo(Inform* msg) {
    (*msg) << Filename_m << " (2D dynamic); zini= " 
           << zbegin_m << " m; zfinal= " << zend_m 
           << " m;" << endl;
}

double FM2DDynamic::getFrequency() const {
    return frequency_m;
}

void FM2DDynamic::setFrequency(double freq) {
    frequency_m = freq;
}

void FM2DDynamic::getOnaxisEz(std::vector<std::pair<double, double>>& F) {
    double dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    F.resize(num_gridpz_m);

    auto Ez = FieldstrengthEz_m.view_host();  // Access host view

    for (int i = 0; i < num_gridpz_m; ++i) {
        F[i].first  = dz * i;
        F[i].second = Ez(i) / 1e6;  // Convert V/m -> MV/m
    }
}

// void FM2DDynamic::getOnaxisEz(std::vector<std::pair<double, double> > & F) {
//     double dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
//     F.resize(num_gridpz_m);

//     for(int i = 0; i < num_gridpz_m; ++ i) {
//         F[i].first = dz * i;
//         F[i].second = FieldstrengthEz_m[i] / 1e6;

//     }
// }
