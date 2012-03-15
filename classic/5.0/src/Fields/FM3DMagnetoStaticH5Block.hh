#ifndef CLASSIC_FIELDMAP3DMAGNETOSTATICH5BLOCK_HH
#define CLASSIC_FIELDMAP3DMAGNETOSTATICH5BLOCK_HH

#include "Fields/Fieldmap.hh"
#include "hdf5.h"
#include "H5PartTypes.h"

using namespace std;

class FM3DMagnetoStaticH5Block: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *msg);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
private:
    FM3DMagnetoStaticH5Block(string aFilename);
    ~FM3DMagnetoStaticH5Block();

    virtual void readMap();
    virtual void freeMap();

    h5part_float64_t *FieldstrengthEz_m;    /**< 3D array with Ez */
    h5part_float64_t *FieldstrengthEx_m;    /**< 3D array with Ex */
    h5part_float64_t *FieldstrengthEy_m;    /**< 3D array with Ey */
    h5part_float64_t *FieldstrengthBz_m;    /**< 3D array with Bz */
    h5part_float64_t *FieldstrengthBx_m;    /**< 3D array with Bx */
    h5part_float64_t *FieldstrengthBy_m;    /**< 3D array with By */

    h5part_float64_t frequency_m;

    h5part_float64_t xbegin_m;
    h5part_float64_t xend_m;
    //     int xcentral_idx_m;

    h5part_float64_t ybegin_m;
    h5part_float64_t yend_m;
    //     int ycentral_idx_m;

    h5part_float64_t zbegin_m;
    h5part_float64_t zend_m;

    h5part_float64_t hx_m;                   /**< length between points in grid, x-direction */
    h5part_float64_t hy_m;                   /**< length between points in grid, y-direction */
    h5part_float64_t hz_m;                   /**< length between points in grid, z-direction */
    h5part_int64_t num_gridpx_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    h5part_int64_t num_gridpy_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    h5part_int64_t num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    bool swap_m;
    friend class Fieldmap;
};

#endif
