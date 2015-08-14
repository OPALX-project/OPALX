#include "Fields/FM3DDynamic.h"
#include "Fields/Fieldmap.hpp"

#include "Physics/Physics.h"

#include <fstream>
#include <ios>

extern Inform *gmsg;

using namespace std;
using Physics::mu_0;

FM3DDynamic::FM3DDynamic(std::string aFilename):
    Fieldmap(aFilename),
    FieldstrengthEz_m(NULL),
    FieldstrengthEx_m(NULL),
    FieldstrengthEy_m(NULL),
    FieldstrengthHz_m(NULL),
    FieldstrengthHx_m(NULL),
    FieldstrengthHy_m(NULL) {

    std::string tmpString;
    double tmpDouble;

    Type = T3DDynamic;
    ifstream file(Filename_m.c_str());

    if(file.good()) {
        bool parsing_passed = interpreteLine<std::string, std::string>(file, tmpString, tmpString);
        if(tmpString == "XYZ") {
            swap_m = false;
            parsing_passed = parsing_passed &&
                             interpreteLine<double>(file, frequency_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, unsigned int>(file, xbegin_m, xend_m, num_gridpx_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, unsigned int>(file, ybegin_m, yend_m, num_gridpy_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, unsigned int>(file, zbegin_m, zend_m, num_gridpz_m);
        } else {
            cerr << "at the moment only the XYZ orientation is supported!\n"
                 << "unknown orientation of 3D dynamic fieldmap" << endl;
            parsing_passed = false;
        }
        for(unsigned long i = 0; (i < (num_gridpz_m + 1) * (num_gridpx_m + 1) * (num_gridpy_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed &&
                             interpreteLine<double>(file,
                                                    tmpDouble,
                                                    tmpDouble,
                                                    tmpDouble,
                                                    tmpDouble,
                                                    tmpDouble,
                                                    tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
        } else {
            frequency_m *= Physics::two_pi * 1e6;

            xbegin_m /= 100.0;
            xend_m /= 100.0;
            ybegin_m /= 100.0;
            yend_m /= 100.0;
            zbegin_m /= 100.0;
            zend_m /= 100.0;

            hx_m = (xend_m - xbegin_m) / num_gridpx_m;
            hy_m = (yend_m - ybegin_m) / num_gridpy_m;
            hz_m = (zend_m - zbegin_m) / num_gridpz_m;

            num_gridpx_m++;
            num_gridpy_m++;
            num_gridpz_m++;

        }
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}


FM3DDynamic::~FM3DDynamic() {
    freeMap();
}

void FM3DDynamic::readMap() {
    if(FieldstrengthEz_m == NULL) {

        ifstream in(Filename_m.c_str());
        int tmpInt;
        std::string tmpString;
        double tmpDouble, Ezmax = 0.0;

        interpreteLine<std::string, std::string>(in, tmpString, tmpString);
        interpreteLine<double>(in, tmpDouble);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        FieldstrengthEz_m = new double[num_gridpz_m * num_gridpx_m * num_gridpy_m];
        FieldstrengthEx_m = new double[num_gridpz_m * num_gridpx_m * num_gridpy_m];
        FieldstrengthEy_m = new double[num_gridpz_m * num_gridpx_m * num_gridpy_m];
        FieldstrengthHz_m = new double[num_gridpz_m * num_gridpx_m * num_gridpy_m];
        FieldstrengthHx_m = new double[num_gridpz_m * num_gridpx_m * num_gridpy_m];
        FieldstrengthHy_m = new double[num_gridpz_m * num_gridpx_m * num_gridpy_m];

        //       if (swap_m)
        //         for (int i = 0; i < num_gridpz_m; i++)
        //           {
        //             for (int j = 0; j < num_gridpr_m; j++)
        //               in >> FieldstrengthEr_m[i + j * num_gridpz_m] >> FieldstrengthEz_m[i + j * num_gridpz_m] >> FieldstrengthHt_m[i + j * num_gridpz_m] >> tmpDouble;
        //             if (fabs(FieldstrengthEz_m[i]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[i]);
        //           }
        //       else
        //         {
        long ii = 0;
        for(unsigned int i = 0; i < num_gridpx_m; ++ i) {
            for(unsigned int j = 0; j < num_gridpy_m; ++ j) {
                for(unsigned int k = 0; k < num_gridpz_m; ++ k) {
                    interpreteLine<double>(in,
                                           FieldstrengthEz_m[ii],
                                           FieldstrengthEy_m[ii],
                                           FieldstrengthEx_m[ii],
                                           FieldstrengthHz_m[ii],
                                           FieldstrengthHy_m[ii],
                                           FieldstrengthHx_m[ii]);
                    ++ ii;
                }
            }
        }
        in.close();

        int index_x = static_cast<int>(ceil(-xbegin_m / hx_m));
        double lever_x = index_x * hx_m + xbegin_m;
        if(lever_x > 0.5) {
            -- index_x;
        }

        int index_y = static_cast<int>(ceil(-ybegin_m / hy_m));
        double lever_y = index_y * hy_m + ybegin_m;
        if(lever_y > 0.5) {
            -- index_y;
        }

        ii = (index_y + index_x * num_gridpy_m) * num_gridpz_m;
        for(unsigned int i = 0; i < num_gridpz_m; i++) {
            if(std::abs(FieldstrengthEz_m[ii]) > Ezmax) {
                Ezmax = std::abs(FieldstrengthEz_m[ii]);
            }
            ++ ii;
        }

        for(unsigned long i = 0; i < num_gridpx_m * num_gridpy_m * num_gridpz_m; i++) {

            FieldstrengthEz_m[i] *= 1e6 / Ezmax;
            FieldstrengthEx_m[i] *= 1e6 / Ezmax;
            FieldstrengthEy_m[i] *= 1e6 / Ezmax;
            FieldstrengthHz_m[i] *= Physics::mu_0 / Ezmax;
            FieldstrengthHx_m[i] *= Physics::mu_0 / Ezmax;
            FieldstrengthHy_m[i] *= Physics::mu_0 / Ezmax;
        }

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
                << endl);
    }
}

void FM3DDynamic::freeMap() {
    if(FieldstrengthEz_m != NULL) {
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEx_m;
        delete[] FieldstrengthEy_m;
        delete[] FieldstrengthHz_m;
        delete[] FieldstrengthHx_m;
        delete[] FieldstrengthHy_m;

        FieldstrengthEz_m = NULL;
        FieldstrengthEx_m = NULL;
        FieldstrengthEy_m = NULL;
        FieldstrengthHz_m = NULL;
        FieldstrengthHx_m = NULL;
        FieldstrengthHy_m = NULL;

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << "\n"
                << endl);
    }
}

bool FM3DDynamic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    const unsigned int index_x = static_cast<int>(floor((R(0) - xbegin_m) / hx_m));
    const double lever_x = (R(0) - xbegin_m) / hx_m - index_x;

    const unsigned int index_y = static_cast<int>(floor((R(1) - ybegin_m) / hy_m));
    const double lever_y = (R(1) - ybegin_m) / hy_m - index_y;

    const unsigned int index_z = (int)floor((R(2) - zbegin_m) / hz_m);
    const double lever_z = (R(2) - zbegin_m) / hz_m - index_z;

    if(index_z >= num_gridpz_m - 2) {
        return false;
    }

    if(index_x >= num_gridpx_m - 2|| index_y >= num_gridpy_m - 2) {
        return true;
    }

    static unsigned int deltaX = num_gridpy_m * num_gridpz_m;
    static unsigned int deltaY = num_gridpz_m;
    static unsigned int deltaZ = 1;

    const unsigned long index1 = index_x * deltaX + index_y * deltaY + index_z * deltaZ;

    E(0) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEx_m[index1                           ]
        + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEx_m[index1 + deltaX                  ]
        + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthEx_m[index1 +          deltaY         ]
        + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthEx_m[index1 + deltaX + deltaY         ]
        + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthEx_m[index1 +                   deltaZ]
        + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthEx_m[index1 + deltaX +          deltaZ]
        + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthEx_m[index1 +          deltaY + deltaZ]
        + lever_x           * lever_y         * lever_z         * FieldstrengthEx_m[index1 + deltaX + deltaY + deltaZ];

    E(1) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEy_m[index1                           ]
        + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEy_m[index1 + deltaX                  ]
        + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthEy_m[index1 +          deltaY         ]
        + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthEy_m[index1 + deltaX + deltaY         ]
        + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthEy_m[index1 +                   deltaZ]
        + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthEy_m[index1 + deltaX +          deltaZ]
        + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthEy_m[index1 +          deltaY + deltaZ]
        + lever_x           * lever_y         * lever_z         * FieldstrengthEy_m[index1 + deltaX + deltaY + deltaZ];

    E(2) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEz_m[index1                           ]
        + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEz_m[index1 + deltaX                  ]
        + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthEz_m[index1 +          deltaY         ]
        + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthEz_m[index1 + deltaX + deltaY         ]
        + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthEz_m[index1 +                   deltaZ]
        + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthEz_m[index1 + deltaX +          deltaZ]
        + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthEz_m[index1 +          deltaY + deltaZ]
        + lever_x           * lever_y         * lever_z         * FieldstrengthEz_m[index1 + deltaX + deltaY + deltaZ];

    B(0) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHx_m[index1                           ]
        + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHx_m[index1 + deltaX                  ]
        + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthHx_m[index1 +          deltaY         ]
        + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthHx_m[index1 + deltaX + deltaY         ]
        + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthHx_m[index1 +                   deltaZ]
        + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthHx_m[index1 + deltaX +          deltaZ]
        + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthHx_m[index1 +          deltaY + deltaZ]
        + lever_x           * lever_y         * lever_z         * FieldstrengthHx_m[index1 + deltaX + deltaY + deltaZ];

    B(1) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHy_m[index1                           ]
        + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHy_m[index1 + deltaX                  ]
        + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthHy_m[index1 +          deltaY         ]
        + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthHy_m[index1 + deltaX + deltaY         ]
        + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthHy_m[index1 +                   deltaZ]
        + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthHy_m[index1 + deltaX +          deltaZ]
        + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthHy_m[index1 +          deltaY + deltaZ]
        + lever_x           * lever_y         * lever_z         * FieldstrengthHy_m[index1 + deltaX + deltaY + deltaZ];

    B(2) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHz_m[index1                           ]
        + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHz_m[index1 + deltaX                  ]
        + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthHz_m[index1 +          deltaY         ]
        + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthHz_m[index1 + deltaX + deltaY         ]
        + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthHz_m[index1 +                   deltaZ]
        + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthHz_m[index1 + deltaX +          deltaZ]
        + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthHz_m[index1 +          deltaY + deltaZ]
        + lever_x           * lever_y         * lever_z         * FieldstrengthHz_m[index1 + deltaX + deltaY + deltaZ];

    return false;
}

bool FM3DDynamic::getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const {
    return false;
}

void FM3DDynamic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = xbegin_m;
    rEnd = xend_m;
}
void FM3DDynamic::getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const {}

void FM3DDynamic::swap() {}

void FM3DDynamic::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (3D dynamic) "
           << " xini= " << xbegin_m << " xfinal= " << xend_m
           << " yini= " << ybegin_m << " yfinal= " << yend_m
           << " zini= " << zbegin_m << " zfinal= " << zend_m << " (m) " << endl;
}

double FM3DDynamic::getFrequency() const {
    return frequency_m;
}

void FM3DDynamic::setFrequency(double freq) {
    frequency_m = freq;
}

void FM3DDynamic::getOnaxisEz(std::vector<std::pair<double, double> > & F) {
    F.resize(num_gridpz_m);

    int index_x = static_cast<int>(ceil(-xbegin_m / hx_m));
    double lever_x = index_x * hx_m + xbegin_m;
    if(lever_x > 0.5) {
        -- index_x;
    }

    int index_y = static_cast<int>(ceil(-ybegin_m / hy_m));
    double lever_y = index_y * hy_m + ybegin_m;
    if(lever_y > 0.5) {
        -- index_y;
    }

    std::ofstream out("data/" + Filename_m);
    unsigned int ii = (index_y + index_x * num_gridpy_m) * num_gridpz_m;
    for(unsigned int i = 0; i < num_gridpz_m; ++ i) {
        F[i].first = hz_m * i;
        F[i].second = FieldstrengthEz_m[ii ++] / 1e6;


        Vector_t R(0,0,zbegin_m + F[i].first), B(0.0), E(0.0);
        getFieldstrength(R, E, B);
        out << std::setw(16) << std::setprecision(8) << F[i].first
            << std::setw(16) << std::setprecision(8) << E(0)
            << std::setw(16) << std::setprecision(8) << E(1)
            << std::setw(16) << std::setprecision(8) << E(2)
            << std::setw(16) << std::setprecision(8) << B(0)
            << std::setw(16) << std::setprecision(8) << B(1)
            << std::setw(16) << std::setprecision(8) << B(2) << "\n";
    }
    out.close();
}