#include <fstream>
#include <ios>

#include "Fields/FM2DMagnetoStatic_cspline.hh"
#include "Fields/Fieldmap.icc"

using namespace std;

FM2DMagnetoStatic_cspline::FM2DMagnetoStatic_cspline(string aFilename):
    Fieldmap(aFilename),
    Bz_interpolants_m(NULL),
    Br_interpolants_m(NULL),
    Bz_accel_m(NULL),
    Br_accel_m(NULL),
    Bz_values_m(NULL),
    Br_values_m(NULL),
    zvals_m(NULL),
    rvals_m(NULL)
{
    Inform msg("FM2DMS ");
    int tmpInt;
    string tmpString;
    double tmpDouble;

    Type = T2DMagnetoStatic_cspline;
    ifstream file(Filename_m.c_str());

    if (file.good()) {
        bool parsing_passed = interpreteLine<string, string>(file, tmpString, tmpString);
        if (tmpString == "ZX") {
            swap_m = true;
            parsing_passed = parsing_passed && 
                interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
            parsing_passed = parsing_passed && 
                interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        } else if (tmpString == "XZ") {
            swap_m = false;
            parsing_passed = parsing_passed && 
                interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
            parsing_passed = parsing_passed && 
                interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        } else {
            cerr << "unknown orientation of 2D magnetostatic fieldmap" << endl;
            parsing_passed = false;
        }

        for (long i = 0; (i < (num_gridpz_m + 1) * (num_gridpr_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed && interpreteLine<double, double>(file, tmpDouble, tmpDouble);
        }
   
        parsing_passed = parsing_passed && 
            interpreteEOF(file);

        file.close();

        if (!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
        } else {            
            rbegin_m /= 100.;
            rend_m /= 100.;
            zbegin_m /= 100.;
            zend_m /= 100.;
            
            hr_m = (rend_m - rbegin_m)/num_gridpr_m;
            hz_m = (zend_m - zbegin_m)/num_gridpz_m;
            
            num_gridpr_m++;
            num_gridpz_m++;
        }
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

FM2DMagnetoStatic_cspline::~FM2DMagnetoStatic_cspline()
{
    if (Bz_interpolants_m != NULL) {
        for (int j = 0; j < num_gridpr_m; ++ j) {
            gsl_spline_free(Bz_interpolants_m[j]);
            gsl_spline_free(Br_interpolants_m[j]);
         
            gsl_interp_accel_free(Bz_accel_m[j]);
            gsl_interp_accel_free(Br_accel_m[j]);
        }
        delete[] Bz_interpolants_m;
        delete[] Br_interpolants_m;

        delete[] Bz_accel_m;
        delete[] Br_accel_m;
     
        delete[] Bz_values_m;
        delete[] Br_values_m;

        delete[] zvals_m;
        delete[] rvals_m;
    }
}

void FM2DMagnetoStatic_cspline::readMap()
{
    if (Bz_interpolants_m == NULL) {
        Inform msg("FM2DMS ");
        ifstream in(Filename_m.c_str());
        int tmpInt;
        string tmpString;
        double tmpDouble, Bzmax = 0.0;

        interpreteLine<string, string>(in, tmpString, tmpString);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        double *FieldstrengthBz = new double[num_gridpz_m * num_gridpr_m];
        double *FieldstrengthBr = new double[num_gridpz_m * num_gridpr_m];
      
        if (swap_m) {
            for (int i = 0; i < num_gridpz_m; i++) {
                for (int j = 0; j < num_gridpr_m; j++) {
                    interpreteLine<double, double>(in,
                                                   FieldstrengthBr[i + j * num_gridpz_m],
                                                   FieldstrengthBz[i + j * num_gridpz_m]);
                }

                if (fabs(FieldstrengthBz[i]) > Bzmax) Bzmax = fabs(FieldstrengthBz[i]);
            } 
        } else {
            for (int j = 0; j < num_gridpr_m; j++) {
                for (int i = 0; i < num_gridpz_m; i++) {
                    interpreteLine<double, double>(in, 
                                                   FieldstrengthBz[i + j * num_gridpz_m],
                                                   FieldstrengthBr[i + j * num_gridpz_m]);
                }
                for (int i = 0; i < num_gridpz_m; i++) {
                    if (fabs(FieldstrengthBz[i]) > Bzmax) Bzmax = fabs(FieldstrengthBz[i]);
                }
            }
        }

        for (int i = 0; i < num_gridpr_m * num_gridpz_m; i++) {
            FieldstrengthBz[i] /= Bzmax;
            FieldstrengthBr[i] /= Bzmax;
        }
        in.close();

        zvals_m = new double[num_gridpz_m];
        for (int i = 0; i < num_gridpz_m; ++ i) {
            zvals_m[i] = (zend_m - zbegin_m) * i / (num_gridpz_m - 1); 
        }
        rvals_m = new double[num_gridpr_m + 3];
        for (int j = -3; j < num_gridpr_m; ++ j) {
            rvals_m[j+3] = (rend_m - rbegin_m) * j / (num_gridpr_m - 1);
        }

        Bz_interpolants_m = new gsl_spline*[num_gridpr_m];
        Br_interpolants_m = new gsl_spline*[num_gridpr_m];

        Bz_accel_m = new gsl_interp_accel*[num_gridpr_m];
        Br_accel_m = new gsl_interp_accel*[num_gridpr_m];

        Bz_values_m = new double[num_gridpr_m + 3];
        Br_values_m = new double[num_gridpr_m + 3];

        for (int j = 0; j < num_gridpr_m; ++ j) {
            Bz_interpolants_m[j] = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
            Br_interpolants_m[j] = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);

            gsl_spline_init(Bz_interpolants_m[j], zvals_m, &FieldstrengthBz[j * num_gridpz_m], num_gridpz_m);
            gsl_spline_init(Br_interpolants_m[j], zvals_m, &FieldstrengthBr[j * num_gridpz_m], num_gridpz_m);
         
            Bz_accel_m[j] = gsl_interp_accel_alloc();
            Br_accel_m[j] = gsl_interp_accel_alloc();
        }

        delete[] FieldstrengthBz;
        delete[] FieldstrengthBr;

        msg << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;
    } 
}

void FM2DMagnetoStatic_cspline::freeMap()
{
    if (Bz_interpolants_m != NULL) {
        Inform msg("FM2DMS ");
        
        for (int j = 0; j < num_gridpr_m; ++ j) {
            gsl_spline_free(Bz_interpolants_m[j]);
            gsl_spline_free(Br_interpolants_m[j]);
            
            gsl_interp_accel_free(Bz_accel_m[j]);
            gsl_interp_accel_free(Br_accel_m[j]);
        }
        delete[] Bz_interpolants_m;
        delete[] Br_interpolants_m;
        
        delete[] Bz_accel_m;
        delete[] Br_accel_m;
        
        delete[] Bz_values_m;
        delete[] Br_values_m;
        
        delete[] zvals_m;
        delete[] rvals_m;

        msg << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;
    }
}

bool FM2DMagnetoStatic_cspline::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
    const double RR = sqrt(R(0)*R(0) + R(1)*R(1));

    for (int j = 0; j < num_gridpr_m; ++ j) {
        Bz_values_m[j+3] = gsl_spline_eval(Bz_interpolants_m[j], R(2), Bz_accel_m[j]);
        Br_values_m[j+3] = gsl_spline_eval(Br_interpolants_m[j], R(2), Br_accel_m[j]);
    }
    for (int j = 1; j < 4; ++ j) {
        Bz_values_m[3-j] = Bz_values_m[3+j];
        Br_values_m[3-j] = -Br_values_m[3+j];
    }

    gsl_interp_accel *Bz_acc = gsl_interp_accel_alloc();
    gsl_interp_accel *Br_acc = gsl_interp_accel_alloc();

    gsl_spline *Bz_spline = gsl_spline_alloc(gsl_interp_cspline, num_gridpr_m + 3);
    gsl_spline *Br_spline = gsl_spline_alloc(gsl_interp_cspline, num_gridpr_m + 3);

    gsl_spline_init(Bz_spline, rvals_m, Bz_values_m, num_gridpr_m + 3);
    gsl_spline_init(Br_spline, rvals_m, Br_values_m, num_gridpr_m + 3);

    double BfieldZ = gsl_spline_eval(Bz_spline, RR, Bz_acc);
    double BfieldR = gsl_spline_eval(Br_spline, RR, Br_acc);

    gsl_spline_free(Bz_spline);
    gsl_spline_free(Br_spline);

    gsl_interp_accel_free(Bz_acc);
    gsl_interp_accel_free(Br_acc);

    if( RR != 0 ) {
        B(0) += BfieldR * R(0)/RR;
        B(1) += BfieldR * R(1)/RR;
    }
    B(2) += BfieldZ;
    return false;
}

bool FM2DMagnetoStatic_cspline::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const 
{
    return false;
}

void FM2DMagnetoStatic_cspline::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM2DMagnetoStatic_cspline::swap()
{
    if (swap_m) swap_m = false;
    else swap_m = true;
}

void FM2DMagnetoStatic_cspline::getInfo(Inform *msg)
{
    (*msg) << Filename_m << " (2D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" <<endl;
}

double FM2DMagnetoStatic_cspline::getFrequency() const
{
    return 0.0;
}

void FM2DMagnetoStatic_cspline::setFrequency(double freq)
{ }
