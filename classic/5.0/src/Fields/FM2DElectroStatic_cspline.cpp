#include <fstream>
#include <ios>

#include "Fields/FM2DElectroStatic_cspline.hh"
#include "Fields/Fieldmap.icc"

using namespace std;

FM2DElectroStatic_cspline::FM2DElectroStatic_cspline(string aFilename):
    Fieldmap(aFilename),
    Ez_interpolants_m(NULL),
    Er_interpolants_m(NULL),
    Ez_accel_m(NULL),
    Er_accel_m(NULL),
    Ez_values_m(NULL),
    Er_values_m(NULL),
    zvals_m(NULL),
    rvals_m(NULL) {
    Inform msg("FM2DES ");
    ifstream file;
    string tmpString;
    double tmpDouble;

    Type =  T2DElectroStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = interpreteLine<string, string>(file, tmpString, tmpString);
        if(tmpString == "ZX") {
            swap_m = true;
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        } else if(tmpString == "XZ") {
            swap_m = false;
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        } else {
            cerr << "unknown orientation of 2D electrostatic fieldmap" << endl;
            parsing_passed = false;
        }

        for(long i = 0; (i < (num_gridpz_m + 1) * (num_gridpr_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed && interpreteLine<double, double>(file, tmpDouble, tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
        } else {
            // conversion from cm to m
            rbegin_m /= 100.;
            rend_m /= 100.;
            zbegin_m /= 100.;
            zend_m /= 100.;

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

FM2DElectroStatic_cspline::~FM2DElectroStatic_cspline() {
    if(Ez_interpolants_m != NULL) {
        for(int j = 0; j < num_gridpr_m; ++ j) {
            gsl_spline_free(Ez_interpolants_m[j]);
            gsl_spline_free(Er_interpolants_m[j]);

            gsl_interp_accel_free(Ez_accel_m[j]);
            gsl_interp_accel_free(Er_accel_m[j]);
        }
        delete[] Ez_interpolants_m;
        delete[] Er_interpolants_m;

        delete[] Ez_accel_m;
        delete[] Er_accel_m;

        delete[] Ez_values_m;
        delete[] Er_values_m;

        delete[] zvals_m;
        delete[] rvals_m;
    }
}

void FM2DElectroStatic_cspline::readMap() {
    if(Ez_interpolants_m == NULL) {
        // declare variables and allocate memory
        Inform msg("FM2DES ");
        ifstream in;
        int tmpInt;
        string tmpString;
        double tmpDouble, Ezmax = 0.0;

        double *FieldstrengthEz = new double[num_gridpz_m * num_gridpr_m];
        double *FieldstrengthEr = new double[num_gridpz_m * num_gridpr_m];

        // read in and parse field map
        in.open(Filename_m.c_str());
        interpreteLine<string, string>(in, tmpString, tmpString);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        if(swap_m) {
            for(int i = 0; i < num_gridpz_m; i++) {
                for(int j = 0; j < num_gridpr_m; j++) {
                    interpreteLine<double, double>(in,
                                                   FieldstrengthEr[i + j * num_gridpz_m],
                                                   FieldstrengthEz[i + j * num_gridpz_m]);
                }
                if(fabs(FieldstrengthEz[i]) > Ezmax) Ezmax = fabs(FieldstrengthEz[i]);
            }
        } else {
            for(int j = 0; j < num_gridpr_m; j++) {
                for(int i = 0; i < num_gridpz_m; i++) {
                    interpreteLine<double, double>(in,
                                                   FieldstrengthEz[i + j * num_gridpz_m],
                                                   FieldstrengthEr[i + j * num_gridpz_m]);
                }
            }

            for(int i = 0; i < num_gridpz_m; i++) {
                if(fabs(FieldstrengthEz[i]) > Ezmax) {
                    Ezmax = fabs(FieldstrengthEz[i]);
                }
            }
        }
        in.close();

        // normalize the field such  that Ez_max = 1 MV/m
        for(int i = 0; i < num_gridpr_m * num_gridpz_m; ++ i) {
            FieldstrengthEz[i] *= 1e6 / Ezmax;  // factor 1e6 due to conversion from MV/m to V/m
            FieldstrengthEr[i] *= 1e6 / Ezmax;
        }
        // calculate the sampling positions
        zvals_m = new double[num_gridpz_m];
        for(int i = 0; i < num_gridpz_m; ++ i) {
            zvals_m[i] = (zend_m - zbegin_m) * i / (num_gridpz_m - 1);
        }
        rvals_m = new double[num_gridpr_m + 3];
        for(int j = -3; j < num_gridpr_m; ++ j) {
            rvals_m[j + 3] = (rend_m - rbegin_m) * j / (num_gridpr_m - 1);
        }

        // allocate array of splines and accelerators
        Ez_interpolants_m = new gsl_spline*[num_gridpr_m];
        Er_interpolants_m = new gsl_spline*[num_gridpr_m];

        Ez_accel_m = new gsl_interp_accel*[num_gridpr_m];
        Er_accel_m = new gsl_interp_accel*[num_gridpr_m];

        // allocate memory to store the results of the 1D spline interpolations
        Ez_values_m = new double[num_gridpr_m + 3];
        Er_values_m = new double[num_gridpr_m + 3];

        // allocate the memory for the splines and accelerators and initialize them
        for(int j = 0; j < num_gridpr_m; ++ j) {
            Ez_interpolants_m[j] = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
            Er_interpolants_m[j] = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);

            gsl_spline_init(Ez_interpolants_m[j], zvals_m, &FieldstrengthEz[j * num_gridpz_m], num_gridpz_m);
            gsl_spline_init(Er_interpolants_m[j], zvals_m, &FieldstrengthEr[j * num_gridpz_m], num_gridpz_m);

            Ez_accel_m[j] = gsl_interp_accel_alloc();
            Er_accel_m[j] = gsl_interp_accel_alloc();
        }

        delete[] FieldstrengthEz;
        delete[] FieldstrengthEr;

        msg << typeset_msg("read in field map '" + Filename_m + "'", "info") << "\n"
            << endl;
    }
}

void FM2DElectroStatic_cspline::freeMap() {
    if(Ez_interpolants_m != NULL) {
        Inform msg("FM2DES ");

        for(int j = 0; j < num_gridpr_m; ++ j) {
            gsl_spline_free(Ez_interpolants_m[j]);
            gsl_spline_free(Er_interpolants_m[j]);

            gsl_interp_accel_free(Ez_accel_m[j]);
            gsl_interp_accel_free(Er_accel_m[j]);
        }
        delete[] Ez_interpolants_m;
        delete[] Er_interpolants_m;

        delete[] Ez_accel_m;
        delete[] Er_accel_m;

        delete[] Ez_values_m;
        delete[] Er_values_m;

        delete[] zvals_m;
        delete[] rvals_m;


        msg << typeset_msg("freed field map '" + Filename_m + "'", "info") << "\n"
            << endl;
    }
}

bool FM2DElectroStatic_cspline::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    const double RR = sqrt(R(0) * R(0) + R(1) * R(1));

    // do 1D cubic spline interpolation along z direction for R = j * dr
    for(int j = 0; j < num_gridpr_m; ++ j) {
        Ez_values_m[j + 3] = gsl_spline_eval(Ez_interpolants_m[j], R(2), Ez_accel_m[j]);
        Er_values_m[j + 3] = gsl_spline_eval(Er_interpolants_m[j], R(2), Er_accel_m[j]);
    }
    // field should be rotational symmetric therefore extend the field accordingly
    // is 4 enough? don't know...
    for(int j = 1; j < 4; ++ j) {
        Ez_values_m[3 - j] = Ez_values_m[3 + j];
        Er_values_m[3 - j] = -Er_values_m[3 + j];
    }

    // allocate spline and accelerators for the spline in the 2nd dimension and initialize them
    gsl_interp_accel *Ez_acc = gsl_interp_accel_alloc();
    gsl_interp_accel *Er_acc = gsl_interp_accel_alloc();

    gsl_spline *Ez_spline = gsl_spline_alloc(gsl_interp_cspline, num_gridpr_m + 3);
    gsl_spline *Er_spline = gsl_spline_alloc(gsl_interp_cspline, num_gridpr_m + 3);

    gsl_spline_init(Ez_spline, rvals_m, Ez_values_m, num_gridpr_m + 3);
    gsl_spline_init(Er_spline, rvals_m, Er_values_m, num_gridpr_m + 3);

    // do interpolation in 2nd dimension
    double EfieldZ = gsl_spline_eval(Ez_spline, RR, Ez_acc);
    double EfieldR = gsl_spline_eval(Er_spline, RR, Er_acc);

    gsl_spline_free(Ez_spline);
    gsl_spline_free(Er_spline);

    gsl_interp_accel_free(Ez_acc);
    gsl_interp_accel_free(Er_acc);

    if(RR > 1e-10) {
        E(0) += EfieldR * R(0) / RR;
        E(1) += EfieldR * R(1) / RR;
    }
    E(2) += EfieldZ;
    return false;
}

bool FM2DElectroStatic_cspline::getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const {
    return false;
}

void FM2DElectroStatic_cspline::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM2DElectroStatic_cspline::getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const {}

void FM2DElectroStatic_cspline::swap() {
    if(swap_m) swap_m = false;
    else swap_m = true;
}

void FM2DElectroStatic_cspline::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (2D electrostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM2DElectroStatic_cspline::getFrequency() const {
    return 0.0;
}

void FM2DElectroStatic_cspline::setFrequency(double freq)
{ ;}
