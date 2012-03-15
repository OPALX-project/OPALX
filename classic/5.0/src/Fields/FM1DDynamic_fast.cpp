#include <fstream>
#include <ios>

#include "Fields/FM1DDynamic_fast.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"
#include "fftw3.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

FM1DDynamic_fast::FM1DDynamic_fast(string aFilename)
    :Fieldmap(aFilename),
     FieldstrengthEz_m(NULL),
     FieldstrengthEr_m(NULL),
     FieldstrengthBt_m(NULL)
{
    Inform msg("FM1DD ");

    Type = T1DDynamic;
    ifstream file(Filename_m.c_str());

    if (file.good()) {
        int tmpInt;
        string tmpString;
        double tmpDouble;

        char buffer[30];

        bool parsing_passed = interpreteLine<string, int>(file, tmpString, tmpInt);
        parsing_passed = parsing_passed &&
            interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        parsing_passed = parsing_passed &&
            interpreteLine<double>(file, frequency_m);
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
            frequency_m *= Physics::two_pi * 1e6;
            
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
        num_gridpz_m = num_gridpr_m = 0;
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

FM1DDynamic_fast::~FM1DDynamic_fast()
{
    if (FieldstrengthEz_m != NULL) {
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEr_m;
        delete[] FieldstrengthBt_m;
    }
}

void FM1DDynamic_fast::readMap()
{
    if (FieldstrengthEz_m == NULL) {
        Inform msg("FM1DD ");
        ifstream in(Filename_m.c_str());
        if (in.good() && num_gridpz_m > 0 && num_gridpr_m > 0) {
            int tmpInt;
            string tmpString;
            double tmpDouble;

            int accuracy;

            int num_gridpzp = (int)floor(num_gridpz_m/2.) + 1;
            double *RealValues;
            fftw_complex* FourCoefs;
            fftw_plan p;
            double ez, ezp, ezpp, ezppp, kz, kz_low, f, fp;
            double xlrep = frequency_m / c;
            double somefactor, somefactor_base;
            double R;
            double Ez_max = 0.0;
            int tmp_num_gridpz = num_gridpz_m;

            FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
            FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];
            FieldstrengthBt_m = new double[num_gridpz_m * num_gridpr_m];

            interpreteLine<string, int>(in, tmpString, accuracy);
            interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
            interpreteLine<double>(in, tmpDouble);
            interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

            RealValues = (double*) fftw_malloc(sizeof(double) * num_gridpz_m);

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
          
            kz_low = 0.0;
            if (fabs(RealValues[0] - RealValues[num_gridpz_m - 1]) > 1.0e-6) {
                double* tmp_values = new double[num_gridpz_m];
                for (int i = 0; i < num_gridpz_m; ++ i) {
                    tmp_values[i] = RealValues[i];
                }
                delete RealValues;
                delete FourCoefs;
                RealValues = (double*) fftw_malloc(2 * sizeof(double) * num_gridpz_m);
                FourCoefs  = (fftw_complex*) fftw_malloc(2 * sizeof(fftw_complex) * num_gridpzp);
              
                for (int i = 0; i < num_gridpz_m; ++ i) {
                    RealValues[num_gridpz_m - 1 -i] = tmp_values[i];
                    RealValues[num_gridpz_m + i] = tmp_values[i];
                }
                delete tmp_values;
                kz_low = Physics::pi;
                tmp_num_gridpz = 2 * num_gridpz_m;
            } else {
                FourCoefs  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_gridpzp);
            }
          

            p = fftw_plan_dft_r2c_1d(tmp_num_gridpz, RealValues, FourCoefs, FFTW_ESTIMATE);
            fftw_execute(p);

            tmpDouble = two_pi / tmp_num_gridpz;
            FourCoefs[0][0] /= tmp_num_gridpz;
            for (int i = 1; i <= accuracy; i++) {
                FourCoefs[i][0] *=  2./tmp_num_gridpz;
                FourCoefs[i][1] *= -2./tmp_num_gridpz;
            }

            for (int i = 0; i < num_gridpz_m; i++) {
                kz = tmpDouble * i + kz_low;
                ez = FourCoefs[0][0];
                ezp = 0.0;
                ezpp = 0.0;
                ezppp = 0.0;
                for (int l = 1; l <= accuracy ; l++) {
                    somefactor_base = l * tmpDouble / hz_m;
                    somefactor = 1.0;
                    const double coskzl = cos(kz * l);
                    const double sinkzl = sin(kz * l);
                    ez    +=              ( FourCoefs[l][0] * coskzl + FourCoefs[l][1] * sinkzl); somefactor *= somefactor_base;
                    ezp   += somefactor * (-FourCoefs[l][0] * sinkzl + FourCoefs[l][1] * coskzl); somefactor *= somefactor_base;
                    ezpp  += somefactor * (-FourCoefs[l][0] * coskzl - FourCoefs[l][1] * sinkzl); somefactor *= somefactor_base;
                    ezppp += somefactor * ( FourCoefs[l][0] * sinkzl - FourCoefs[l][1] * coskzl);
                }

                f  = -(ezpp  + ez *  xlrep * xlrep)/16.;
                fp = -(ezppp + ezp * xlrep * xlrep)/16.;
                for (int j = 0; j < num_gridpr_m; j++) {
                    R = j * hr_m;
                    FieldstrengthEz_m[i + j * num_gridpz_m] =   ez + 4 * f * R * R;
                    FieldstrengthEr_m[i + j * num_gridpz_m] = -(ezp/2. + fp * R * R);
                    FieldstrengthBt_m[i + j * num_gridpz_m] =  (ez/2. + f * R * R) * xlrep / c;
                }
            }

#ifdef _DEBUG
            {
                Vector_t R = Vector_t(0.0, 0.0, zbegin_m);
                Vector_t E(0.0);
                Vector_t B(0.0);
                size_t pos_dot = Filename_m.find_last_of('.');
                string fnm = Filename_m.substr(0,pos_dot) + ".dbg";
                ofstream dbg(fnm.c_str());
                double dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
                for (R(2) = zbegin_m; R(2) < zend_m; R(2) += dz) {
                    E = 0.0;
                    B = 0.0;
                    getFieldstrength(R, E, B);
                    dbg << R(2) << "\t" << E(2) << endl;
                }
            }
#endif

            fftw_destroy_plan(p);
            fftw_free(FourCoefs);
            fftw_free(RealValues);

            msg << typeset_msg("read in fieldmap '" + Filename_m + "'", "info") << "\n"
                << endl;

            in.close();
        }
    }
}

void FM1DDynamic_fast::freeMap()
{
    if (FieldstrengthEz_m != NULL) {
        Inform msg("FM1DD ");
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEr_m;
        delete[] FieldstrengthBt_m;
        
        msg << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;
    }
}


bool FM1DDynamic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
    //   if (FieldstrengthEz_m == NULL) readMap();

    const double RR = sqrt(R(0)*R(0) + R(1)*R(1));

    const int indexr = (int)floor(RR / hr_m);
    const double leverr = (RR / hr_m) - indexr;

    const int indexz = (int)floor((R(2)) / hz_m);
    const double leverz = (R(2) / hz_m) - indexz;

    if ((indexz < 0) || (indexz + 2 > num_gridpz_m) || (indexr < 0) || (indexr + 2 > num_gridpr_m)) {
        //cerr << getName() << ".getFieldstrength(): out of boudaries (z,r) = (" << R(2) << "," << RR << ")" << endl ;
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

    const double BfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBt_m[index1]
        + leverz         * (1.0 - leverr) * FieldstrengthBt_m[index1 + 1]
        + (1.0 - leverz) * leverr         * FieldstrengthBt_m[index2]
        + leverz         * leverr         * FieldstrengthBt_m[index2 + 1];

    E(0) += EfieldR * R(0);
    E(1) += EfieldR * R(1);
    E(2) += EfieldZ;
    B(0) += -BfieldT * R(1);
    B(1) +=  BfieldT * R(0);

    return false;
}

bool FM1DDynamic_fast::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{ }


void FM1DDynamic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM1DDynamic_fast::swap()
{ }

void FM1DDynamic_fast::getInfo(Inform *msg)
{
    (*msg) << Filename_m << " (1D dynamic fast); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM1DDynamic_fast::getFrequency() const
{
    return frequency_m;
}

void FM1DDynamic_fast::setFrequency(double freq)
{
    frequency_m = freq;
}
