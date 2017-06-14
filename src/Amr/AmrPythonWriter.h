#ifndef AMR_PYTHON_WRITER_H
#define AMR_PYTHON_WRITER_H

#include "Amr/AbstractAmrWriter.h"

#include <boost/filesystem.hpp>

/*!
 * This concrete class writes AMR grids
 * and particles that can be visualized by
 * amrPlots/part_boxes.py that is part of
 * the pyOPALTools repository.
 */
class AmrPythonWriter : public AbstractAmrWriter {
    
public:
    
    AmrPythonWriter();
    
    void writeGrids(const std::string& dir,
                    const amr::AmrFieldContainer_t& rho,
                    const amr::AmrFieldContainer_t& phi,
                    const amr::AmrFieldContainer_t& efield,
                    const amr::AmrIntArray_t& refRatio,
                    const amr::AmrGeomContainer_t& geom,
                    const double& time);
    
    void writeBunch(const AmrPartBunch* bunch_p);
    
private:
    boost::filesystem::path dir_m;  ///< directory where to write files
};

#endif
