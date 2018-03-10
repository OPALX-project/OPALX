#ifndef ABSTRACT_AMR_WRITER_H
#define ABSTRACT_AMR_WRITER_H

#include "Amr/AmrObject.h"
#include "Amr/AmrDefs.h"
#include "Algorithms/AmrPartBunch.h"

/*!
 * Abstract base class for writing AMR data to files in
 * order to plot them.
 */
class AbstractAmrWriter {
    
public:
    /*!
     * @param rho is the charge density on all levels
     * @param phi is the electrostatic potential on all levels
     * @param efield are the electric field components on all levels
     * @param refRatio are the refinement ratios among the levels
     * @param geom are the geometries of all levels
     * @param nLevel available
     * @param time specifies the step.
     * @param scale used for mapping
     */
    virtual void writeFields(const amr::AmrFieldContainer_t& rho,
                             const amr::AmrFieldContainer_t& phi,
                             const amr::AmrFieldContainer_t& efield,
                             const amr::AmrIntArray_t& refRatio,
                             const amr::AmrGeomContainer_t& geom,
                             const int& nLevel,
                             const double& time,
                             const double& scale = 1.0) = 0;
    
    /*!
     * @param bunch_p
     * @param time
     * @param scale used for mapping
     */
    virtual void writeBunch(const AmrPartBunch* bunch_p,
                            const double& time,
                            const double& scale = 1.0) = 0;
    
    virtual ~AbstractAmrWriter() { }
    
};

#endif
