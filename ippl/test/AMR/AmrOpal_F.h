#ifndef _Adv_F_H_
#define _Adv_F_H_
#include <BLFort.H>

/*!
 * @file AmrOpal_F.h
 * @author Weiqun Zhang
 * @date October 2016, LBNL
 * @details Fortran - C++ interface for tagging
 * @brief Fortran - C++ interface for calling the tagging function
 */

extern "C" 
{
//     void initdata(const int& level, const Real& time, 
// 		  const int* lo, const int* hi,
// 		  BL_FORT_FAB_ARG_3D(state),
// 		  const Real* dx, const Real* problo);

//     void get_face_velocity(const int& level, const Real& time, 
// 			   D_DECL(BL_FORT_FAB_ARG(xvel),
// 				  BL_FORT_FAB_ARG(yvel),
// 				  BL_FORT_FAB_ARG(zvel)),
// 			   const Real* dx, const Real* problo);

    /*!
     * Calls the Fortran routine for tagging.
     * @param tag is the integer tag array
     * @param tag_lo is the lower index extent of tag array
     * @param tag_hi is the upper index extent of tag array
     * @param state is the state array
     * @param tagval is the integer value to tag cell for refinement
     * @param clearval is the integer value to untag cell
     * @param lo is the lower left corner of the work region
     * @param hi is the upper right corner of the work region
     * @param dx is the cell size
     * @param problo is the physical location of lower left corner of prob domain
     * @param time is the problem evolution time (not used)
     * @param phierr is the number of particles per cell
     */
    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
		     const BL_FORT_FAB_ARG_3D(state),
		     const int* tagval, const int* clearval,
		     const int* lo, const int* hi,
		     const Real* dx, const Real* problo,
		     const Real* time, const Real* phierr);

//     void advect(const Real& time, const int* lo, const int*hi,
// 		const BL_FORT_FAB_ARG_3D(statein),
// 		BL_FORT_FAB_ARG_3D(stateout),
// 		D_DECL(const BL_FORT_FAB_ARG_3D(xvel),
// 		       const BL_FORT_FAB_ARG_3D(yvel),
// 		       const BL_FORT_FAB_ARG_3D(zvel)),
// 		D_DECL(BL_FORT_FAB_ARG_3D(fx),
// 		       BL_FORT_FAB_ARG_3D(fy),
// 		       BL_FORT_FAB_ARG_3D(fz)),
// 		const Real* dx, const Real& dt);
}

#endif
