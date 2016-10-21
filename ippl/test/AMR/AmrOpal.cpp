#include "AmrOpal.h"


AmrOpal::AmrOpal() { }

AmrOpal::~AmrOpal() { }


void AmrOpal::ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ngrow*/) { }

void AmrOpal::ErrorEst(int lev, MultiFab& mf, TagBoxArray& tags, Real time) {

//     const int clearval = TagBox::CLEAR;
//     const int   tagval = TagBox::SET;
// 
//     const Real* dx      = geom[lev].CellSize();
//     const Real* prob_lo = geom[lev].ProbLo();
//     int nPart = 2;
// 
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//     {
//         Array<int>  itags;
// 	
// 	for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
// 	{
// 	    const Box&  tilebx  = mfi.tilebox();
// 
//             TagBox&     tagfab  = tags[mfi];
// 	    
// 	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
// 	    // So we are going to get a temporary integer array.
// 	    tagfab.get_itags(itags, tilebx);
// 	    
//             // data pointer and index space
// 	    int*        tptr    = itags.dataPtr();
// 	    const int*  tlo     = tilebx.loVect();
// 	    const int*  thi     = tilebx.hiVect();
// 
// 	    state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
// 			BL_TO_FORTRAN_3D(mf[mfi]),
// 			&tagval, &clearval, 
// 			ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
// 			ZFILL(dx), ZFILL(prob_lo), &time, nPart);
// 	    //
// 	    // Now update the tags in the TagBox.
// 	    //
// 	    tagfab.tags_and_untags(itags, tilebx);
// 	}
//     }
}
