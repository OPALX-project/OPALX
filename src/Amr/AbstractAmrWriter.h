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
     * @param time specifies the step.
     */
    virtual void writeGrids(const amr::AmrFieldContainer_t& rho,
                            const amr::AmrFieldContainer_t& phi,
                            const amr::AmrFieldContainer_t& efield,
                            const amr::AmrIntArray_t& refRatio,
                            const amr::AmrGeomContainer_t& geom,
                            const double& time) = 0;
    
    /*!
     * @param bunch_p
     */
    virtual void writeBunch(const AmrPartBunch* bunch_p) = 0;
    
    virtual ~AbstractAmrWriter() { }
    
};

#endif
