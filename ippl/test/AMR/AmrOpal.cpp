#include "AmrOpal.h"
#include "AmrOpal_F.h"

// AmrOpal::AmrOpal() { }
AmrOpal::AmrOpal(const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord, PartBunchBase* bunch)
    : AmrCore(rb, max_level_in, n_cell_in, coord), bunch_m(dynamic_cast<AmrPartBunch*>(bunch))
{
//     nPartPerCell_m.resize(max_level_in, PArrayManage);
}

AmrOpal::~AmrOpal() { }


void AmrOpal::initBaseLevel() {
    const Real time = 0.0;

    // define coarse level BoxArray and DistributionMap
    {
        finest_level = 0;
        const BoxArray& ba = MakeBaseGrids();
        DistributionMapping dm(ba, ParallelDescriptor::NProcs());

        SetBoxArray(0 /*level*/, ba);
        SetDistributionMap(0 /*level*/, dm);
    }

//     if (max_level > 0) // build fine levels
//     {
//         Array<BoxArray> new_grids(max_level+1);
//         new_grids[0] = grids[0];
//         do
//         {
//             int new_finest;
//             
//             // Add (at most) one level at a time.
//             MakeNewGrids(finest_level,time,new_finest,new_grids);
// 
//             if (new_finest <= finest_level) break;
//             finest_level = new_finest;
//             
//             DistributionMapping dm(new_grids[new_finest], ParallelDescriptor::NProcs());
// 
//             SetBoxArray(new_finest /*level*/, new_grids[new_finest]);
//             SetDistributionMap(new_finest, dm);
// //             MakeNewLevel(new_finest, time, new_grids[new_finest], dm);
// 
// //             InitLevelData(new_finest);
// 	}
// 	while (finest_level < max_level);
//     }
    
    
}

void AmrOpal::initFineLevels() {
    
    
}


void AmrOpal::doTagging(int lev, MultiFab& mf, TagBoxArray& tags) {
    
    
}


void AmrOpal::ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ngrow*/) {
    
    
    BoxArray ba = this->boxArray(lev);
    MultiFab nPartPerCell(ba, 1, 1);
    bunch_m->AssignDensitySingleLevel(lev, nPartPerCell, 0, 1, 0);
    
//     nPartPerCell_m.set(lev, new MultiFab(ba, 1 , 1));
//     bunch_m->AssignDensitySingleLevel(lev, nPartPerCell_m[lev], 0, 1, 0);
    
    
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real nPart = 2;
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        
        for (MFIter mfi(nPartPerCell,true); mfi.isValid(); ++mfi) {
            const Box&  tilebx  = mfi.tilebox();
            
            TagBox&     tagfab  = tags[mfi];
            
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
            tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D(nPartPerCell[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &nPart);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
    
}


void
AmrOpal::regrid (int lbase, Real time)
{
    int new_finest;
    Array<BoxArray> new_grids(finest_level+2);
    MakeNewGrids(lbase, time, new_finest, new_grids);

    BL_ASSERT(new_finest <= finest_level+1);

    DistributionMapping::FlushCache();

    for (int lev = lbase+1; lev <= new_finest; ++lev)
    {
	if (lev <= finest_level) // an old level
	{
	    if (new_grids[lev] != grids[lev]) // otherwise nothing
	    {
		DistributionMapping new_dmap(new_grids[lev], ParallelDescriptor::NProcs());
		RemakeLevel(lev, time, new_grids[lev], new_dmap);
	    }
	}
	else  // a new level
	{
	    DistributionMapping new_dmap(new_grids[lev], ParallelDescriptor::NProcs());
	    MakeNewLevel(lev, time, new_grids[lev], new_dmap);
	}
    }
}


void
AmrOpal::RemakeLevel (int lev, Real time,
		     const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
/*    const int ncomp = phi_new[lev]->nComp();
    const int nghost = phi_new[lev]->nGrow();

    auto new_state = std::unique_ptr<MultiFab>(new MultiFab(new_grids, ncomp, nghost, new_dmap));
    auto old_state = std::unique_ptr<MultiFab>(new MultiFab(new_grids, ncomp, nghost, new_dmap));

    FillPatch(lev, time, *new_state, 0, ncomp);

    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev] = std::unique_ptr<FluxRegister>
	    (new FluxRegister(grids[lev], refRatio(lev-1), lev, ncomp, dmap[lev]));
    } */   
}

void
AmrOpal::MakeNewLevel (int lev, Real time,
		      const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
//     const int ncomp = 1;
//     const int nghost = 0;
// 
//     SetBoxArray(lev, new_grids);
//     SetDistributionMap(lev, new_dmap);
// 
//     phi_new[lev] = std::unique_ptr<MultiFab>(new MultiFab(grids[lev], ncomp, nghost, dmap[lev]));
//     phi_old[lev] = std::unique_ptr<MultiFab>(new MultiFab(grids[lev], ncomp, nghost, dmap[lev]));
// 
//     t_new[lev] = time;
//     t_old[lev] = time - 1.e200;
// 
//     if (lev > 0 && do_reflux) {
// 	flux_reg[lev] = std::unique_ptr<FluxRegister>
// 	    (new FluxRegister(grids[lev], refRatio(lev-1), lev, ncomp, dmap[lev]));
//     }
}
