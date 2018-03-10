#ifndef AMR_SLICE_WRITER_H
#define AMR_SLICE_WRITER_H

#include "Amr/AbstractAmrWriter.h"

#include <boost/filesystem.hpp>

class AmrSliceWriter : public AbstractAmrWriter {
    
public:
    
    AmrSliceWriter();
    
    /*!
     * This method writes only a slice of the potential along the
     * x-axis centered in y- and z-direction. It works only for
     * single-core. Additionally the maximum grid size has to be
     * equal to the number of grid points per dimension.
     */
    void writeFields(const amr::AmrFieldContainer_t& rho,
                     const amr::AmrFieldContainer_t& phi,
                     const amr::AmrFieldContainer_t& efield,
                     const amr::AmrIntArray_t& refRatio,
                     const amr::AmrGeomContainer_t& geom,
                     const double& time,
                     const double& scale);
    
    /*!
     * @param bunch_p
     */
    void writeBunch(const AmrPartBunch* bunch_p,
                    const double& time,
                    const double& scale);
    
private:
    boost::filesystem::path dir_m;      ///< directory where to write files
    std::string filename_m;             ///< of input file
};

#endif
