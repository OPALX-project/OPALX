#ifndef OPAL_LASERPROFILE_HH
#define OPAL_LaserProfile_HH
// ------------------------------------------------------------------------
// $RCSfile: LaserProfile.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LaserProfile
//
// ------------------------------------------------------------------------
// #include <Ippl.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram2d.h>
#include <string>
// #include <hdf5.h>


class LaserProfile {

public:
    LaserProfile() {
    }

    LaserProfile(std::string fn, std::string image, double cut) {
        ReadFile(fn, image, cut);
    }

    ~LaserProfile() {
        gsl_histogram2d_pdf_free(pdf_m);
        gsl_histogram2d_free(hist2d_m);
        gsl_rng_free(r_m);
    }
    void ReadFile(std::string fn, std::string image, double cut);

    void SaveDist();
    void SampleDist();
    void GetXY(double *s_x, double *s_y);

#ifdef LASERPROFILE_TEST
    void GetX();
#endif

    // void BackGroundCut(double cut );
    void GetProfileMax(unsigned short int *profileMax_m, unsigned short int   *image);

private:
    double X, Y;
    //unsigned short int profileMax_m;
    int sizeX_m, sizeY_m;
    gsl_histogram2d *hist2d_m;
    const gsl_rng_type *rngT_m;
    gsl_rng *r_m;
    gsl_histogram2d_pdf *pdf_m;

};
#endif

