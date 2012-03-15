#include <fstream>
#include <ios>

#include "Fields/FM1DMagnetoStatic.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"
#include "fftw3.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

FM1DMagnetoStatic::FM1DMagnetoStatic(string aFilename)
    :Fieldmap(aFilename),
     realFourCoefs_m(NULL),
     imagFourCoefs_m(NULL)
{
    int tmpInt;
    string tmpString;
    double tmpDouble;
    Inform msg("FM1DMS ");
    ifstream file(Filename_m.c_str());

    Type = T1DMagnetoStatic;

    if (file.good()) {
        bool parsing_passed = interpreteLine<string, int>(file, tmpString, accuracy_m);
        parsing_passed = parsing_passed &&
            interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        parsing_passed = parsing_passed &&
            interpreteLine<double, double, int>(file, rbegin_m, rend_m, tmpInt);
        for (int i = 0; (i <= num_gridpz_m) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed && 
                interpreteLine<double>(file, tmpDouble);
        }

        parsing_passed = parsing_passed && 
            interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if (!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
        } else {
            rbegin_m /= 100.;
            rend_m /= 100.;
            zbegin_m /= 100.;
            zend_m /= 100.;

            num_gridpz_m++;
        }
        length_m = 2.0 * num_gridpz_m * (zend_m - zbegin_m) / (num_gridpz_m - 1);
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

FM1DMagnetoStatic::~FM1DMagnetoStatic()
{
    if (realFourCoefs_m != NULL) {
        delete[] realFourCoefs_m;
        delete[] imagFourCoefs_m;
    }
}

void FM1DMagnetoStatic::readMap()
{

#ifdef _DEBUG
    double rend = 0.0;
    int num_gridpr = 0;
#endif
    if (realFourCoefs_m == NULL) {
        Inform msg("FM1DMS ");
        ifstream in(Filename_m.c_str());
        
        int tmpInt;
        string tmpString;
        double tmpDouble;

        double *RealValues;
        fftw_complex* FourCoefs;
        fftw_plan p;
        double ez, ezp, ezpp, ezppp, kz;
        double somefactor, somefactor_base;
        double Bz_max = 0.0;

        interpreteLine<string, int>(in, tmpString, accuracy_m);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
#ifdef _DEBUG
        interpreteLine<double, double, int>(in, tmpDouble, rend, num_gridpr);
        rend /= 100.0;
#else
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
#endif

        realFourCoefs_m = new double[accuracy_m];
        imagFourCoefs_m = new double[accuracy_m - 1];

        RealValues = (double*) fftw_malloc(sizeof(double) * 2*num_gridpz_m);
        FourCoefs  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (num_gridpz_m + 1));

        for (int i = 0; i < num_gridpz_m; i++) {
            interpreteLine<double>(in, RealValues[num_gridpz_m + i]);
            if (fabs(RealValues[num_gridpz_m + i]) > Bz_max) {
                Bz_max = fabs(RealValues[num_gridpz_m + i]);
            }
            RealValues[num_gridpz_m - 1 - i] = RealValues[num_gridpz_m + i];
        }
        in.close();

        num_gridpz_m *= 2;
        for (int i = 0; i < num_gridpz_m; ++ i) {
            RealValues[i] /= Bz_max;
        }

        p = fftw_plan_dft_r2c_1d(num_gridpz_m, RealValues, FourCoefs, FFTW_ESTIMATE);
        fftw_execute(p);

        realFourCoefs_m[0] = FourCoefs[0][0] / num_gridpz_m ;
        for (int i = 1; i < accuracy_m; i++) {
            realFourCoefs_m[i] = 2.* FourCoefs[i][0] /  num_gridpz_m;
            imagFourCoefs_m[i - 1] = 2.* FourCoefs[i][1] / num_gridpz_m;
        }

        fftw_destroy_plan(p);
        fftw_free(FourCoefs);
        fftw_free(RealValues);

        msg << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;
    }

#ifdef _DEBUG    
    double dr = rend / num_gridpr;
    double dz = length_m / (num_gridpz_m - 1);
    Vector_t R = Vector_t(0.0, 0.0, zbegin_m);
    Vector_t Ef(0.0), Bf(0.0);
    size_t pos_dot = Filename_m.find_last_of('.');
    ofstream dbg((Filename_m.substr(0, pos_dot) + ".dbg").c_str());
    for (int j = 0; j <= num_gridpr; ++ j) {
        R(2) = zbegin_m;
        for (int i = 0; i < num_gridpz_m/2; ++ i) {
            Ef = 0.0;
            Bf = 0.0;
            getFieldstrength(R, Ef, Bf);
            dbg << Bf(2) << "\t"
                << Bf(0) << endl;
            R(2) += dz;
        }
        R(0) += dr;
    }
    dbg.close();
#endif

//         Vector_t Ef(0.0), Bf(0.0), R(0.01, 0.0, zbegin_m);
//         double DZ = length_m / 100;
//         for (int i = 0; i < 100; ++ i) {
//             Ef = 0.0;
//             Bf = 0.0;
//             getFieldstrength(R, Ef, Bf);
//             R(2) += DZ;
//             cout << Bf(0) << "\t"
//                  << Bf(1) << "\t"
//                  << Bf(2) << endl;
//         }

}

void FM1DMagnetoStatic::freeMap()
{
    if (realFourCoefs_m != NULL) {

        Inform msg("FM1DMS ");

        delete[] realFourCoefs_m;
        delete[] imagFourCoefs_m;

        msg << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;

    }
}

bool FM1DMagnetoStatic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
    const double RR2 = R(0)*R(0) + R(1)*R(1);

    const double kz = two_pi * R(2)/length_m + Physics::pi;

    double ez = realFourCoefs_m[0];
    double ezp = 0.0;
    double ezpp = 0.0;
    double ezppp = 0.0;
    double somefactor_base, somefactor;
    double coskzl;
    double sinkzl;

    for (int l = 1; l < accuracy_m ; l++) {
        somefactor_base =  two_pi / length_m * l;
        somefactor = 1.0;
        coskzl = cos(kz * l);
        sinkzl = sin(kz * l);
        ez    +=              ( realFourCoefs_m[l] * coskzl - imagFourCoefs_m[l-1] * sinkzl); 
        somefactor *= somefactor_base;
        ezp   += somefactor * (-realFourCoefs_m[l] * sinkzl - imagFourCoefs_m[l-1] * coskzl); 
        somefactor *= somefactor_base;
        ezpp  += somefactor * (-realFourCoefs_m[l] * coskzl + imagFourCoefs_m[l-1] * sinkzl); 
        somefactor *= somefactor_base;
        ezppp += somefactor * ( realFourCoefs_m[l] * sinkzl + imagFourCoefs_m[l-1] * coskzl);
    }

    const double BfieldR = -ezp/2. + ezppp / 16. * RR2;

    B(0) += BfieldR * R(0);
    B(1) += BfieldR * R(1);
    B(2) += ez - ezpp * RR2 / 4.;
    return false;
}

bool FM1DMagnetoStatic::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{ }

void FM1DMagnetoStatic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM1DMagnetoStatic::swap()
{ }

void FM1DMagnetoStatic::getInfo(Inform *msg)
{
    (*msg) << Filename_m << " (1D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM1DMagnetoStatic::getFrequency() const
{
    return 0.0;
}

void FM1DMagnetoStatic::setFrequency(double freq)
{ }
