#include <fstream>
#include <ios>

#include "Fields/FM1DElectroStatic_fast.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"
#include "fftw3.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;
using Physics::pi;

FM1DElectroStatic_fast::FM1DElectroStatic_fast(string aFilename)
    :Fieldmap(aFilename),
     FieldstrengthEz_m(NULL),
     FieldstrengthEr_m(NULL)
{
    Inform msg("FM1DES ");
    int tmpInt;
    string tmpString;
    double tmpDouble;

    Type = T1DElectroStatic;
    ifstream file(Filename_m.c_str());

    if (file.good()) {
        bool parsing_passed = interpreteLine<string, int>(file, tmpString, tmpInt);
        parsing_passed = parsing_passed &&
            interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        parsing_passed = parsing_passed &&
            interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        for (int i = 0; (i <= num_gridpz_m) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed && 
                interpreteLine<double>(file, tmpDouble);
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


FM1DElectroStatic_fast::~FM1DElectroStatic_fast()
{
    if (FieldstrengthEz_m != NULL) {
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEr_m;
    }
}

void FM1DElectroStatic_fast::readMap()
{
    if (FieldstrengthEz_m == NULL) {
        Inform msg("FM1DES ");
        ifstream in(Filename_m.c_str());

        int tmpInt;
        string tmpString;
        double tmpDouble;

        int accuracy;

        double *RealValues;
        fftw_complex* FourCoefs;
        fftw_plan p;
        double ez, ezp, ezpp, ezppp, kz;
        double somefactor, somefactor_base;
        double R;
        double Ez_max = 0.0;

        interpreteLine<string, int>(in, tmpString, accuracy);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        RealValues = (double*) fftw_malloc(sizeof(double) * 2*num_gridpz_m);
        FourCoefs  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (num_gridpz_m + 1));

        for (int i = 0; i < num_gridpz_m; i++) {
            interpreteLine<double>(in, RealValues[num_gridpz_m + i]);
            if (fabs(RealValues[num_gridpz_m + i]) > Ez_max ) {
                Ez_max = fabs(RealValues[num_gridpz_m + i]);
            }
            RealValues[num_gridpz_m - 1 - i] = RealValues[num_gridpz_m + i];
        }
        in.close();

        for (int i = 0; i < 2*num_gridpz_m; ++ i) {
            RealValues[i] /= Ez_max;
        }
        
        p = fftw_plan_dft_r2c_1d(2*num_gridpz_m, RealValues, FourCoefs, FFTW_ESTIMATE);
        fftw_execute(p);

        FourCoefs[0][0] /= (2*num_gridpz_m);
        for (int i = 1; i < accuracy; i++) {
            FourCoefs[i][0] /= num_gridpz_m;
            FourCoefs[i][1] /= num_gridpz_m;
        }

        FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];

        tmpDouble = pi / num_gridpz_m;
        for (int i = 0; i < num_gridpz_m; i++) {
            kz = tmpDouble * i + pi;
            ez = FourCoefs[0][0];
            ezp = 0.0;
            ezpp = 0.0;
            ezppp = 0.0;
            for (int l = 1; l < accuracy ; l++) {
                somefactor_base = tmpDouble / hz_m * l;
                somefactor = 1.0;
                ez    +=              ( FourCoefs[l][0] * cos(kz * l) - FourCoefs[l][1] * sin(kz * l));somefactor *= somefactor_base;
                ezp   += somefactor * (-FourCoefs[l][0] * sin(kz * l) - FourCoefs[l][1] * cos(kz * l));somefactor *= somefactor_base;
                ezpp  += somefactor * (-FourCoefs[l][0] * cos(kz * l) + FourCoefs[l][1] * sin(kz * l));somefactor *= somefactor_base;
                ezppp += somefactor * ( FourCoefs[l][0] * sin(kz * l) + FourCoefs[l][1] * cos(kz * l));
            }
            R = 0.0;
            for (int j = 0; j < num_gridpr_m; j++) {
                FieldstrengthEz_m[i + j * num_gridpz_m] =  ez - ezpp * R * R / 4.;
                FieldstrengthEr_m[i + j * num_gridpz_m] = -ezp/2. + ezppp * R * R / 16.;
                R = hr_m * (j + 1);
            }
        }

        fftw_destroy_plan(p);
        fftw_free(FourCoefs);
        fftw_free(RealValues);

        msg << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;
    }
}

void FM1DElectroStatic_fast::freeMap()
{
    if (FieldstrengthEz_m != NULL) {
        Inform msg("FM1DD ");

        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEr_m;

        msg << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;

    }
}

bool FM1DElectroStatic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
    const double RR = sqrt(R(0)*R(0) + R(1)*R(1));

    const int indexr = (int)floor(RR / hr_m);
    const double leverr = (RR / hr_m) - indexr;

    const int indexz = (int)floor((R(2)) / hz_m);
    const double leverz = (R(2) / hz_m) - indexz;

    if ((indexz < 0) || (indexz + 2 > num_gridpz_m) || (indexr < 0) || (indexr + 2 > num_gridpr_m)){
        return true;
    }

    const int index1 = indexz + indexr * num_gridpz_m;
    const int index2 = index1 + num_gridpz_m;
    const double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1]
        + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index1 + 1]
        + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index2]
        + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];

    const double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1]
        + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index1 + 1]
        + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index2]
        + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

    E(0) += EfieldR * R(0);
    E(1) += EfieldR * R(1);
    E(2) += EfieldZ;
    return false;
}

bool FM1DElectroStatic_fast::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{ }


void FM1DElectroStatic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM1DElectroStatic_fast::swap()
{ }

void FM1DElectroStatic_fast::getInfo(Inform *msg)
{
    (*msg) << Filename_m << " (1D electrostatic fast); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM1DElectroStatic_fast::getFrequency() const
{
    return 0.0;
}

void FM1DElectroStatic_fast::setFrequency(double freq)
{ }
