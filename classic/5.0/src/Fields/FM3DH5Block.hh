#ifndef CLASSIC_FIELDMAP3DH5BLOCK_HH
#define CLASSIC_FIELDMAP3DH5BLOCK_HH

#include "Fields/Fieldmap.hh"
#include "hdf5.h"
#include "H5hut.h"
#include <vector>

class FM3DH5Block: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *msg);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
    virtual void getOnaxisEz(std::vector<std::pair<double, double> > & F);

private:
    FM3DH5Block(std::string aFilename);
    ~FM3DH5Block();

    virtual void readMap();
    virtual void freeMap();

    std::vector<h5_float64_t> FieldstrengthEz_m;    /**< 3D array with Ez */
    std::vector<h5_float64_t> FieldstrengthEx_m;    /**< 3D array with Ex */
    std::vector<h5_float64_t> FieldstrengthEy_m;    /**< 3D array with Ey */
    std::vector<h5_float64_t> FieldstrengthHz_m;    /**< 3D array with Hz */
    std::vector<h5_float64_t> FieldstrengthHx_m;    /**< 3D array with Hx */
    std::vector<h5_float64_t> FieldstrengthHy_m;    /**< 3D array with Hy */

    h5_float64_t frequency_m;

    h5_float64_t xbegin_m;
    h5_float64_t xend_m;
    //     int xcentral_idx_m;

    h5_float64_t ybegin_m;
    h5_float64_t yend_m;
    //     int ycentral_idx_m;

    h5_float64_t zbegin_m;
    h5_float64_t zend_m;

    h5_float64_t hx_m;                   /**< length between points in grid, x-direction */
    h5_float64_t hy_m;                   /**< length between points in grid, y-direction */
    h5_float64_t hz_m;                   /**< length between points in grid, z-direction */
    h5_int64_t num_gridpx_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    h5_int64_t num_gridpy_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    h5_int64_t num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    bool swap_m;
    friend class Fieldmap;
};

#endif
