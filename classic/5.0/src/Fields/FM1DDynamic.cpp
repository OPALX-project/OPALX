#include <fstream>
#include <ios>

#include "Fields/FM1DDynamic.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"
#include "gsl/gsl_fft_real.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

FM1DDynamic::FM1DDynamic(string aFilename):
    Fieldmap(aFilename),
    FourCoefs_m(NULL)
{
    Inform msg("FM1DD ");

    Type = T1DDynamic;

    ifstream file(Filename_m.c_str());

    if (file.good()) {
        int tmpInt;
        string tmpString;
        double tmpDouble;

        bool parsing_passed = interpreteLine<string, int>(file, tmpString, accuracy_m);
        parsing_passed = parsing_passed &&
            interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        parsing_passed = parsing_passed &&
            interpreteLine<double>(file, frequency_m);
        parsing_passed = parsing_passed &&
            interpreteLine<double, double, int>(file, rbegin_m, rend_m, tmpInt);
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

            frequency_m *= two_pi * 1e6;
            xlrep_m = frequency_m / c;
            
            rbegin_m /= 100.;
            rend_m /= 100.;
            zbegin_m /= 100.;
            zend_m /= 100.;
            
            num_gridpz_m++;
            
            trafo_phase_m = 0.0;
        }
        length_m = num_gridpz_m * (zend_m - zbegin_m) / (num_gridpz_m - 1);
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

FM1DDynamic::~FM1DDynamic()
{
    if (FourCoefs_m != NULL) {
        delete[] FourCoefs_m;
    }
}

void FM1DDynamic::readMap()
{
    if (FourCoefs_m == NULL) {
        Inform msg("FM1DD ");
        ifstream in(Filename_m.c_str());

        string tmpString;

        int tmpInt;
        int num_gridpzp = static_cast<int>(floor(num_gridpz_m/2.) + 1.0);

        double tmpDouble;
        double Ez_max = 0.0;
        
        double *RealValues = new double[num_gridpz_m];

        gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(2*num_gridpz_m);
        gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(2*num_gridpz_m);

        FourCoefs_m = new double[2*accuracy_m - 1];

        interpreteLine<string, int>(in, tmpString, accuracy_m);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double>(in, tmpDouble);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        for (int i = 0; i < num_gridpz_m; ++ i) {
            interpreteLine<double>(in, RealValues[i]);
            if (fabs(RealValues[i]) > Ez_max) {
                Ez_max = fabs(RealValues[i]);
            }
        }
        in.close();

        for (int i = 0; i < num_gridpz_m; ++ i) {
            RealValues[i] /= Ez_max;
        }

        if (fabs(RealValues[0] - RealValues[num_gridpz_m - 1]) > 1.0e-6) {
            double* tmp_values = new double[num_gridpz_m];
            for (int i = 0; i < num_gridpz_m; ++ i) {
                tmp_values[i] = RealValues[i];
            }
            delete RealValues;
            RealValues = new double[2*num_gridpz_m];
            
            for (int i = 0; i < num_gridpz_m; ++ i) {
                RealValues[num_gridpz_m - 1 -i] = tmp_values[i];
                RealValues[num_gridpz_m + i] = tmp_values[i];
            }
            delete tmp_values;

            num_gridpz_m *= 2;
            trafo_phase_m = Physics::pi;
            length_m *= 2.0;
        }

        gsl_fft_real_transform(RealValues, 1, num_gridpz_m, real, work);

        FourCoefs_m[0] = RealValues[0] / num_gridpz_m;
        for (int i = 1; i < 2*accuracy_m-1; i++) {
            FourCoefs_m[i] = 2. * RealValues[i] / num_gridpz_m;
        }

        gsl_fft_real_workspace_free(work);
        gsl_fft_real_wavetable_free(real);

        delete[] RealValues;

        msg << typeset_msg("read in fieldmap '" + Filename_m + "'", "info") << "\n"
            << endl;

#ifdef _DEBUG
        Vector_t R = Vector_t(0.0, 0.0, zbegin_m);
        Vector_t E(0.0);
        Vector_t B(0.0);
        size_t pos_dot = Filename_m.find_last_of('.');
        string fnm = Filename_m.substr(0,pos_dot) + ".dbg";
        ofstream dbg(fnm.c_str());
        double dz = length_m / (num_gridpz_m - 1);
        msg << "* zbegin_m = " << zbegin_m << "; \n"
            << "* zend_m = " << zend_m << "; dz = " << dz << "\n"
            << "* num_gridpz_m = " << num_gridpz_m << endl;
        for (R(2) = zbegin_m; R(2) < zend_m; R(2) += dz) {
            E = 0.0;
            B = 0.0;
            getFieldstrength(R, E, B);
            dbg << R(2) << "\t" << E(2) << endl;
        }
#endif  
    }
}

void FM1DDynamic::freeMap()
{
    if (FourCoefs_m != NULL) {
        Inform msg("FM1DD ");

        delete[] FourCoefs_m;

        msg << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;
    }
}

bool FM1DDynamic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
    const double RR2 = R(0)*R(0) + R(1)*R(1);

    const double kz = two_pi * R(2)/length_m + trafo_phase_m;

    double ez = FourCoefs_m[0];
    double ezp = 0.0;
    double ezpp = 0.0;
    double ezppp = 0.0;
    double somefactor_base, somefactor;
    double coskzl;
    double sinkzl;
    
    int n = 1;
    for (int l = 1; l < accuracy_m ; l++, n += 2) {
        somefactor_base = l * two_pi / length_m;       // = \frac{dkz}{dz}
        somefactor = 1.0;
        coskzl = cos(kz * l);
        sinkzl = sin(kz * l);
        ez    +=              ( FourCoefs_m[n] * coskzl + FourCoefs_m[n+1] * sinkzl); somefactor *= somefactor_base;
        ezp   += somefactor * (-FourCoefs_m[n] * sinkzl + FourCoefs_m[n+1] * coskzl); somefactor *= somefactor_base;
        ezpp  += somefactor * (-FourCoefs_m[n] * coskzl - FourCoefs_m[n+1] * sinkzl); somefactor *= somefactor_base;
        ezppp += somefactor * ( FourCoefs_m[n] * sinkzl - FourCoefs_m[n+1] * coskzl);
    }
    const double f  = -(ezpp  + ez *  xlrep_m * xlrep_m)/16.;
    const double fp = -(ezppp + ezp * xlrep_m * xlrep_m)/16.;

    const double EfieldR = -(ezp/2. + fp * RR2);
    const double BfieldT =  (ez/2. + f * RR2) * xlrep_m / c;

    E(0) +=  EfieldR * R(0);
    E(1) +=  EfieldR * R(1);
    E(2) +=  ez + 4 * f * RR2;
    B(0) += -BfieldT * R(1);
    B(1) +=  BfieldT * R(0);

    return false;
}

bool FM1DDynamic::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{ }

void FM1DDynamic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM1DDynamic::swap()
{ }

void FM1DDynamic::getInfo(Inform *msg)
{
    (*msg) << Filename_m << " (1D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM1DDynamic::getFrequency() const
{
    return frequency_m;
}

void FM1DDynamic::setFrequency(double freq)
{
    frequency_m = freq;
}
