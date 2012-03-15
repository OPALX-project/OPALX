#include <fstream>
#include <ios>

#include "Fields/FM2DDynamic.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"

using namespace std;
using Physics::mu_0;

FM2DDynamic::FM2DDynamic(string aFilename)
    :Fieldmap(aFilename),
     FieldstrengthEz_m(NULL),
     FieldstrengthEr_m(NULL),
     FieldstrengthHt_m(NULL)
{
    Inform msg("FM2DD ");
    int tmpInt;
    string tmpString;
    double tmpDouble;

    Type = T2DDynamic;
    ifstream file(Filename_m.c_str());

    if (file.good()) {
        bool parsing_passed = interpreteLine<string, string>(file, tmpString, tmpString);
        if (tmpString == "ZX") {
            swap_m = true;
            parsing_passed = parsing_passed && 
                interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
            parsing_passed = parsing_passed &&
                interpreteLine<double>(file, frequency_m);
            parsing_passed = parsing_passed && 
                interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        } else if (tmpString == "XZ") {
            swap_m = false;
            parsing_passed = parsing_passed && 
                interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
            parsing_passed = parsing_passed && 
                interpreteLine<double>(file, frequency_m);
            parsing_passed = parsing_passed && 
                interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        }
        else {
            cerr << "unknown orientation of 2D dynamic fieldmap" << endl;
            parsing_passed = false;
            zbegin_m = 0.0;
            zend_m = -1e-3;
        }

        for (long i = 0; (i < (num_gridpz_m + 1) * (num_gridpr_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed && interpreteLine<double, double, double, double>(file, tmpDouble, tmpDouble, tmpDouble, tmpDouble);
        }

        parsing_passed = parsing_passed && 
            interpreteEOF(file);

        file.close();

        if (!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
        } else {
          
            frequency_m *= Physics::two_pi * 1e6;
          
            rbegin_m /= 100.0;
            rend_m /= 100.0;
            zbegin_m /= 100.0;
            zend_m /= 100.0;
          
            hr_m = (rend_m - rbegin_m) / num_gridpr_m;
            hz_m = (zend_m - zbegin_m) / num_gridpz_m;
          
            num_gridpr_m++;
            num_gridpz_m++;
        }
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}


FM2DDynamic::~FM2DDynamic()
{
    if (FieldstrengthEz_m != NULL) {
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEr_m;
        delete[] FieldstrengthHt_m;
    }
}

void FM2DDynamic::readMap()
{
    if (FieldstrengthEz_m == NULL) {
        Inform msg("FM2DD ");
        ifstream in(Filename_m.c_str());
        int tmpInt;
        string tmpString;
        double tmpDouble, Ezmax = 0.0;

        interpreteLine<string, string>(in, tmpString, tmpString);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double>(in, tmpDouble);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthHt_m = new double[num_gridpz_m * num_gridpr_m];

        if (swap_m) {
            for (int i = 0; i < num_gridpz_m; i++) {
                for (int j = 0; j < num_gridpr_m; j++) {
                    interpreteLine<double, double, double, double>(in,
                                                                   FieldstrengthEr_m[i + j * num_gridpz_m],
                                                                   FieldstrengthEz_m[i + j * num_gridpz_m],
                                                                   FieldstrengthHt_m[i + j * num_gridpz_m],
                                                                   tmpDouble);
                }
                if (fabs(FieldstrengthEz_m[i]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[i]);
            }
        } else {
            for (int j = 0; j < num_gridpr_m; j++) {
                for (int i = 0; i < num_gridpz_m; i++) {
                    interpreteLine<double, double, double, double>(in,
                                                                   FieldstrengthEz_m[i + j * num_gridpz_m],
                                                                   FieldstrengthEr_m[i + j * num_gridpz_m],
                                                                   tmpDouble,
                                                                   FieldstrengthHt_m[i + j * num_gridpz_m]);
                }
            }
            for (int i = 0; i < num_gridpz_m; i++)
                if (fabs(FieldstrengthEz_m[i]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[i]);
        }
        
        for (int i = 0; i < num_gridpr_m * num_gridpz_m; i++) {
            FieldstrengthEz_m[i] /= Ezmax;
            FieldstrengthEr_m[i] /= Ezmax;
            FieldstrengthHt_m[i] /= Ezmax;
        }
        in.close();

        msg << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;

    }
}

void FM2DDynamic::freeMap()
{
    if (FieldstrengthEz_m != NULL)
    {
        Inform msg("FM2DD ");
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEr_m;
        delete[] FieldstrengthHt_m;

        msg << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << "\n"
            << endl;
    }
}

bool FM2DDynamic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const
{
    const double RR = sqrt(R(0)*R(0) + R(1)*R(1));

    const int indexr = (int)floor(RR / hr_m);
    const double leverr = RR / hr_m - indexr;

    const int indexz = (int)floor(R(2) / hz_m);
    const double leverz = R(2) / hz_m - indexz;

    if ((indexz < 0) || (indexz + 2 > num_gridpz_m))
        return false;
    if (indexr + 2 > num_gridpr_m)
        return true;

    const int index1 = indexz + indexr * num_gridpz_m;
    const int index2 = index1 + num_gridpz_m;

    double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1]
        + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index1 + 1]
        + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index2]
        + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];

    double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1]
        + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index1 + 1]
        + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index2]
        + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

    double HfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthHt_m[index1]
        + leverz         * (1.0 - leverr) * FieldstrengthHt_m[index1 + 1]
        + (1.0 - leverz) * leverr         * FieldstrengthHt_m[index2]
        + leverz         * leverr         * FieldstrengthHt_m[index2 + 1];

    if( RR > 1e-10 ) {
        E(0) += EfieldR * R(0)/RR;
        E(1) += EfieldR * R(1)/RR;
        B(0) += -mu_0 * HfieldT * R(1) / RR * 1e-6;
        B(1) +=  mu_0 * HfieldT * R(0) / RR * 1e-6;
    }
    E(2) += EfieldZ;
    
    return false;
}

bool FM2DDynamic::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{

}

void FM2DDynamic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const
{
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM2DDynamic::swap()
{
    if (swap_m) swap_m = false;
    else swap_m = true;
}

void FM2DDynamic::getInfo(Inform *msg)
{
    (*msg) << Filename_m << " (2D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM2DDynamic::getFrequency() const
{
    return frequency_m;
}

void FM2DDynamic::setFrequency(double freq)
{
    frequency_m = freq;
}
