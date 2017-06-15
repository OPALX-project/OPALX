#ifndef AMR_YT_WRITER_H
#define AMR_YT_WRITER_H

#include "Amr/AbstractAmrWriter.h"

#include <boost/filesystem.hpp>

/*!
 * This concrete class writes output files
 * that are readable by
 * <a href="http://yt-project.org/">yt</a>.
 * The format can be read by the fork
 * <a href="https://Yavin_4@bitbucket.org/Yavin_4/opal-yt">opal-yt</a>
 * that uses accelerator units instead of astrophysics units.
 * The functions are copied from AMReX and modified to fit
 * our needs.
 */
class AmrYtWriter : public AbstractAmrWriter {
    
public:
    
    AmrYtWriter(int step);
    
    void writeGrids(const amr::AmrFieldContainer_t& rho,
                    const amr::AmrFieldContainer_t& phi,
                    const amr::AmrFieldContainer_t& efield,
                    const amr::AmrIntArray_t& refRatio,
                    const amr::AmrGeomContainer_t& geom,
                    const double& time);
    
    
    void writeBunch(const AmrPartBunch* bunch_p);
    
private:
    boost::filesystem::path dir_m;  ///< directory where to write files
    int step_m;                     ///< that we write
};

#endif
