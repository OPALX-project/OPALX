// ------------------------------------------------------------------------
// $RCSfile: LaserProfile.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LaserProfile
//
// ------------------------------------------------------------------------
#include "Distribution/LaserProfile.h"
#include "Utility/Inform.h"

using namespace std;

void LaserProfile::ReadFile(string fn, string imagestr, double cut) {

    Inform m("LaserProfile::ReadFile ");
    hid_t h5 = H5Fopen(fn.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gr = H5Gopen2 (h5, imagestr.c_str(), H5P_DEFAULT);

    hid_t dsetBit = H5Dopen2 (gr, "CameraBits", H5P_DEFAULT);
    dsetBit = H5Dopen2 (gr, "CameraGain", H5P_DEFAULT);

    double bitData[] = { -1 };
    H5Dread(dsetBit, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bitData);
    int  bits = static_cast<int>(bitData[0]);


    hsize_t dim[2];
    hid_t   dset = H5Dopen2 (gr, "CameraData", H5P_DEFAULT);

    hid_t filespace = H5Dget_space(dset);
    H5Sget_simple_extent_dims(filespace, dim, NULL);
    filespace = H5Screate_simple(2, dim, NULL);
    /*
    filespace = H5Screate_simple(3, dim, NULL);
    */

    sizeX_m        = dim[0];
    sizeY_m        = dim[1];
    /*
    int nImage     = dim[0];
    sizeX_m        = dim[1];
    sizeY_m        = dim[2];
    */

    // selection of hyperslab in hdf5 file:
    hsize_t start[]  = {0, 0}; // Start of hyperslab
    hsize_t count[]  = {sizeX_m, sizeY_m};   // Block count
    /*
    hsize_t start[]  = {imageNumber,0,0};  // Start of hyperslab
    hsize_t count[]  = {imageNumber,sizeX_m,sizeY_m};    // Block count
    */

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

    unsigned short int  *image = new unsigned short int  [sizeX_m*sizeY_m];
    hid_t mem = H5Screate_simple(2, count, NULL);
    /*
    hid_t mem = H5Screate_simple(3, count, NULL);
    */

    H5Dread(dset, H5T_NATIVE_USHORT, mem, filespace, H5P_DEFAULT, image);

    m << "# Read image done .... sizeX= " << sizeX_m << " sizeY= " << sizeY_m << " Bits= " << bits << endl;

    int pixel = 0;
    hsize_t col = 1;
    hsize_t row = 1;

    hist2d_m = gsl_histogram2d_alloc(sizeX_m, sizeY_m);
    gsl_histogram2d_set_ranges_uniform(hist2d_m, 0.0, 1.0, 0.0, 1.0);   // all bins are set to 0 too

    cout << "sizeX_m=" << sizeX_m << "; sizeY_m=" << sizeY_m << endl;

    bool doBackGroundCut_m = 1;
    unsigned short int profileMax_m;
    if(doBackGroundCut_m)
        GetProfileMax(&profileMax_m, image);

    //    exit ();

    /*

    for (col=1; col<=sizeX_m; col++) {
    for (row=1; row<=sizeY_m; row++) {
      //        int val = (image[pixel++] >> shift); // lowest bits are insignificant
      int val = (image[pixel++]);
        //      if ((200 < col < 300) && (200 < row < 300))
        //      gsl_histogram2d_accumulate (hist2d_m, col/sizeX_m, row/sizeY_m, val);
        if ( (700<col)&&(col<1000)&&(1<row)&&(row<300) ) {
        gsl_histogram2d_accumulate (hist2d_m, double(col)/sizeX_m, double(row)/sizeY_m, val);
        if( ((col==5)||(col==250)||(col==500)||(col==800))&&((row==100)||(row==300)||(row==500)||(row==700)) ) {
            cout << col << " " << row << " val=" << val << endl;
            cout << "col/sizeX_m=" << double(col)/sizeX_m << "; row/sizeY_m=" << double(row)/sizeY_m << endl;
        }
        }
    }
    }
    cout << "xmax=" << gsl_histogram2d_xmax(hist2d_m) << "; xmin=" << gsl_histogram2d_xmin(hist2d_m) << endl;
    cout << "ymax=" << gsl_histogram2d_ymax(hist2d_m) << "; ymin=" << gsl_histogram2d_ymin(hist2d_m) << endl;
    cout << "max_val=" << gsl_histogram2d_max_val(hist2d_m) << "; min_val=" << gsl_histogram2d_min_val(hist2d_m) << endl;
    cout << "xmean=" << gsl_histogram2d_xmean(hist2d_m) << "; xsigma=" << gsl_histogram2d_xsigma(hist2d_m) << endl;
    cout << "ymean=" << gsl_histogram2d_ymean(hist2d_m) << "; ysigma=" << gsl_histogram2d_ysigma(hist2d_m) << endl;
    cout << "sum=" << gsl_histogram2d_sum(hist2d_m) << endl;

    */

    /*

    if (doBackGroundCut_m) {
      GetProfileMax();
      double wkjdfkjdf=0.14;
      for (col=1; col<=sizeX_m; col++) {
        for (row=1; row<=sizeY_m; row++) {

            val =

        // Background cut starts here
        if ( (val*1.0)/(MAX*1.0) < wkjdfkjdf )
            val = 0;
        if ( (val*1.0)/(profileMax_m*1.0) >= wkjdfkjdf )
            val = ((val*1.0)/profileMax_m - wkjdfkjdf)/(1-wkjdfkjdf) * profileMax_m;
        // end of cut

        if ( (700<col)&&(col<1000)&&(1<row)&&(row<300) ) {
          gsl_histogram2d_accumulate (hist2d_m, double(col)/sizeX_m, double(row)/sizeY_m, val);
          if( ((col==5)||(col==250)||(col==500)||(col==800))&&((row==100)||(row==300)||(row==500)||(row==700)) ) {
            cout << col << " " << row << " val=" << val << endl;
            cout << "col/sizeX_m=" << double(col)/sizeX_m << "; row/sizeY_m=" << double(row)/sizeY_m << endl;
          }
        }
        }
      }
      cout << "profileMax = " << profileMax << endl;
    }

    */

    /*

    else {

    */

    double insCut = 0.45;

    //  double MAX = gsl_histogram2d_max_val(hist2d_m); // MAX =Maximum im Laserprofil
    //  cout << "profileMax_m = " << profileMax_m << endl;
    for(col = sizeX_m; col >= 1; col--) {
        for(row = 1; row <= sizeY_m; row++) {
            //
            //        if (pixel == 0)
            //          cout << "pixel == 0; image[pixel] = " << image[pixel] << "; image[pixel+1] = " << image[pixel+1] << endl;

            // int val = (image[pixel++] >> shift); // lowest bits are insignificant
            double val;
            /*
                    if (pixel == 0)
                      { val = (image[pixel++]);
                        cout << "val = " << val << "; image[pixel] = " << image[pixel] << endl; }
                    else
                      val = (image[pixel++]);
            */
            val = (image[pixel++]);
            // Background cut starts here
            if((val * 1.0) / (profileMax_m * 1.0) < insCut)
                val = 0;
            if((val * 1.0) / (profileMax_m * 1.0) >= insCut)
                val = ((val * 1.0) / profileMax_m - insCut) / (1 - insCut) * profileMax_m;
            // end of cut



            if((700 < col) && (col < 1000) && (1 < row) && (row < 300)) {
                gsl_histogram2d_accumulate(hist2d_m, double(col) / sizeX_m, double(row) / sizeY_m, val);
                //          if( ((col==5)||(col==250)||(col==500)||(col==800))&&((row==100)||(row==300)||(row==500)||(row==700)) ) {
                //            cout << col << " " << row << " val=" << val << endl;
                //            cout << "col/sizeX_m=" << double(col)/sizeX_m << "; row/sizeY_m=" << double(row)/sizeY_m << endl;
                //          }
            }



            //      if ((200 < col < 300) && (200 < row < 300))
            //          gsl_histogram2d_accumulate (hist2d_m, col/sizeX_m, row/sizeY_m, val);
        }
    }

    /*

    }

    */


    SaveDist();

    gsl_rng_env_setup();

    const gsl_rng_type *T = gsl_rng_default;
    r_m = gsl_rng_alloc(T);

    pdf_m = gsl_histogram2d_pdf_alloc(sizeX_m, sizeY_m);
    gsl_histogram2d_pdf_init(pdf_m, hist2d_m);

    SampleDist();

    double aa, bb;
    GetXY(&aa, &bb);
    cout << "x aa=" << aa << "; y bb=" << bb << endl;

}



/*-----*/
void LaserProfile::SaveDist() {
    FILE  *fh = fopen("dist2dhist.dat", "w");
    gsl_histogram2d_fprintf(fh, hist2d_m, "%g", "%g");
    fclose(fh);
}


void LaserProfile::SampleDist() {
    FILE  *fh = fopen("distsampled.dat", "w");
    double x, y;

    for(int i = 0; i < 10000; i++) {
        GetXY(&x, &y);
        fprintf(fh, "%g %g \n", x, y);
    }


    fclose(fh);
}

void LaserProfile::GetXY(double *s_x, double *s_y) {
    double x, y;
    double u = gsl_rng_uniform(r_m);
    double v = gsl_rng_uniform(r_m);
    gsl_histogram2d_pdf_sample(pdf_m, u, v, &x, &y);

    *s_x = x;
    *s_y = y;

    //  *s_x = x*sizeX_m;
    //  *s_y = y*sizeY_m;
}

void LaserProfile::GetProfileMax(unsigned short int *profileMax_m, unsigned short int image[]) {

    int image_len = sizeX_m * sizeY_m;
    //*profileMax_m = image[i];       // start with max = first element

    *profileMax_m = 0;

    cout << " profileMax_m = " << profileMax_m << " ; *profileMax_m = " << *profileMax_m << endl;
    for(int i = 0; i < image_len; i++) {
        if(image[i] > *profileMax_m)
            *profileMax_m = image[i];
    }
}


#ifdef LASERPROFILE_TEST
/*---------------------------------------------*/
/*---------------------------------------------*/
/*-----*/
double LaserProfile::GetX() {
    double x, y;
    double u = gsl_rng_uniform(r_m);
    double v = gsl_rng_uniform(r_m);
    gsl_histogram2d_pdf_sample(pdf_m, u, v, &x, &y);
    printf("%g %g\n", x * sizeX_m, y * sizeY);

    return x * sizeX_m;
}

/*-----*/
double LaserProfile::GetY() {
    double x, y;
    double u = gsl_rng_uniform(r_m);
    double v = gsl_rng_uniform(r_m);
    gsl_histogram2d_pdf_sample(pdf_m, u, v, &x, &y);
    printf("%g %g\n", x * sizeX_m, y * sizeY_m);

    return y * sizeY_m;
}


int main(int argc, const char *argv[]) {

    LaserProfile test(argc, argv);

    test.GetX();
    test.GetX();
    test.GetX();
    test.GetX();
    test.GetX();
    test.GetX();
    test.GetX();

    cout << "\n";

    test.GetY();
    test.GetY();
    test.GetY();
    test.GetY();
    test.GetY();
    test.GetY();

    SaveDist();

    double a, b;
    test.GetXY(&a, &b);
    cout << "\n";
    cout << a << "\t" << b << "\n";

    return 0;
}
#endif

/*---------------------------------------------*/
/*---------------------------------------------*/
