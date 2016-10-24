#include "AmrOpal.h"
#include "AmrOpal_F.h"
#include <PlotFileUtil.H>


// AmrOpal::AmrOpal() { }
AmrOpal::AmrOpal(const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord, PartBunchBase* bunch)
    : AmrCore(rb, max_level_in, n_cell_in, coord),
    bunch_m(dynamic_cast<AmrPartBunch*>(bunch))
{
    initBaseLevel();
    nPartPerCell_m.resize(max_level_in + 1, PArrayManage);
    
    nPartPerCell_m.set(0, new MultiFab(this->boxArray(0), 1, 1));
}

AmrOpal::~AmrOpal() { }



/* init the base level by using the distributionmap and BoxArray of the
 * particles (called in constructor of AmrOpal
 */
void AmrOpal::initBaseLevel() {
    finest_level = 0; // AmrCore protected member variable
    const BoxArray& ba = bunch_m->ParticleBoxArray(0 /*level*/);
    const DistributionMapping& dm = bunch_m->ParticleDistributionMap(0 /*level*/);
    
    SetBoxArray(0 /*level*/, ba);
    SetDistributionMap(0 /*level*/, dm);
}

// void AmrOpal::initFineLevels() { }


void AmrOpal::writePlotFile(std::string filename, int step) {
    
    std::vector<std::string> varnames(1, "rho");
    
    Array<const MultiFab*> tmp(nPartPerCell_m.size());
    for (int i = 0; i < nPartPerCell_m.size(); ++i)
        tmp[i] = &nPartPerCell_m[i];
    
    const auto& mf = tmp;
    
//     std::cout << "Size: " << nPartPerCell_m.size() << std::endl;
//     std::cout << "Finest level: " << finest_level << std::endl;
    
    Array<int> istep(this->maxLevel(), step);
    
    BoxLib::WriteMultiLevelPlotifle(filename, finest_level, mf, varnames,
                                    Geom(), 0.0, istep, refRatio());
}


void AmrOpal::ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ngrow*/) {
    
    bunch_m->AssignDensitySingleLevel(0, nPartPerCell_m[lev], lev);
    
//     BoxArray ba = this->boxArray(lev);
//     
//     std::cerr << "DEFINE level = " << lev << std::endl;
//     nPartPerCell_m.set(lev, new MultiFab(ba, 1, 1));
    
    // rho_index, sub_cycle, PArray<MultiFab>& mf, lev_min = 0, ncomp, finest_level = -1
//     bunch_m->AssignDensity(0, false, nPartPerCell_m, 0, 1, lev);
    
    std::cout << nPartPerCell_m[lev].min(0) << " " << nPartPerCell_m[lev].max(0) << std::endl;
//     std::cout << nPartPerCell_m[lev].minIndex(0) << " " << nPartPerCell_m[lev].maxIndex(0) << std::endl;
    
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real nPart = 1;
    
//     std::cout << "dx = " << *dx << " prob_lo = " << *prob_lo << std::endl;
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        
        for (MFIter mfi(nPartPerCell_m[lev],true); mfi.isValid(); ++mfi) {
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
                        BL_TO_FORTRAN_3D(nPartPerCell_m[lev][mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &nPart);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
    std::cerr << "FINEST: " << finest_level << std::endl;
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
                
                /*
                 * particles need to know the BoxArray
                 * and DistributionMapping
                 */
                bunch_m->SetParticleBoxArray(lev, new_grids[lev]);
                bunch_m->SetParticleDistributionMap(lev, new_dmap);
	    }
	}
	else  // a new level
	{
	    DistributionMapping new_dmap(new_grids[lev], ParallelDescriptor::NProcs());
	    MakeNewLevel(lev, time, new_grids[lev], new_dmap);
            
            /*
             * particles need to know the BoxArray
             * and DistributionMapping
             */
            bunch_m->SetParticleBoxArray(lev, new_grids[lev]);
            bunch_m->SetParticleDistributionMap(lev, new_dmap);
	}
    }
    
    std::cout << "NEW: " << new_finest << " OLD: " << finest_level << std::endl;
    
    if (new_finest > finest_level)
        finest_level = new_finest;
    
    // update to multilevel
    bunch_m->myUpdate();
}


void
AmrOpal::RemakeLevel (int lev, Real time,
		     const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
//     const int ncomp = phi_new[lev]->nComp();
//     const int nghost = phi_new[lev]->nGrow();

//     auto new_state = std::unique_ptr<MultiFab>(new MultiFab(new_grids, ncomp, nghost, new_dmap));
//     auto old_state = std::unique_ptr<MultiFab>(new MultiFab(new_grids, ncomp, nghost, new_dmap));

//     FillPatch(lev, time, *new_state, 0, ncomp);
    
//     if ( new_grids.empty() )
//         return;
    
    
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nPartPerCell_m.set(lev, new MultiFab(new_grids, 1, 1));

//     std::cerr << nPartPerCell_m[lev].is_nodal() << std::endl;
    
//     std::swap(new_state, phi_new[lev]);
//     std::swap(old_state, phi_old[lev]);

//     t_new[lev] = time;
//     t_old[lev] = time - 1.e200;

//     if (lev > 0 && do_reflux) {
// 	flux_reg[lev] = std::unique_ptr<FluxRegister>
// 	    (new FluxRegister(grids[lev], refRatio(lev-1), lev, ncomp, dmap[lev]));
//     }
}

void
AmrOpal::MakeNewLevel (int lev, Real time,
		      const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
//     const int ncomp = 1;
//     const int nghost = 0;

    
    
//     if ( new_grids.empty() )
//         return;
    
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nPartPerCell_m.set(lev, new MultiFab(new_grids, 1, 1));
    
//     std::cerr << nPartPerCell_m[lev].is_nodal() << std::endl;

//     phi_new[lev] = std::unique_ptr<MultiFab>(new MultiFab(grids[lev], ncomp, nghost, dmap[lev]));
//     phi_old[lev] = std::unique_ptr<MultiFab>(new MultiFab(grids[lev], ncomp, nghost, dmap[lev]));

//     t_new[lev] = time;
//     t_old[lev] = time - 1.e200;

//     if (lev > 0 && do_reflux) {
// 	flux_reg[lev] = std::unique_ptr<FluxRegister>
// 	    (new FluxRegister(grids[lev], refRatio(lev-1), lev, ncomp, dmap[lev]));
//     }
}
