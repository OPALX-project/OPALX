#ifndef AMR_PYTHON_WRITER_H
#define AMR_PYTHON_WRITER_H

#include "Amr/AbstractAmrWriter.h"

#include <boost/filesystem.hpp>

/*!
 * This concrete class writes AMR grids
 * and particles that can be visualized by
 * amrPlots/bunch_and_grids.py that is part of
 * the pyOPALTools repository.
 */
class AmrPythonWriter : public AbstractAmrWriter {
    
public:
    
    AmrPythonWriter();
    
    void writeFields(const amr::AmrFieldContainer_t& rho,
                     const amr::AmrFieldContainer_t& phi,
                     const amr::AmrFieldContainer_t& efield,
                     const amr::AmrIntArray_t& refRatio,
                     const amr::AmrGeomContainer_t& geom,
                     const double& time,
                     const double& scale);
    
    /*!
     * Write the particle coordinates and momenta to a file
     * bunch_**********.dat where * represents a digit for
     * specifying the step number. The grids of all levels
     * are stored in the appropriate file grids_**********.dat.
     * All files are stored in the simulation subdirectory
     * data/amr/python
     * 
     * @param bunch_p
     */
    void writeBunch(const AmrPartBunch* bunch_p,
                    const double& time,
                    const double& scale);
    
private:
    boost::filesystem::path dir_m;  ///< directory where to write files
};

#endif
