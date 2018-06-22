#include "AmrOpal.h"
#include "AmrOpal_F.h"
#include <AMReX_PlotFileUtil.H>
// #include <MultiFabUtil.H>

#include "Solver.h"

#include <map>


// AmrOpal::AmrOpal() { }
AmrOpal::AmrOpal(const amrex::RealBox* rb, int max_level_in, const amrex::Vector<int>& n_cell_in, int coord,
#ifdef IPPL_AMR
                 PartBunchAmr<amrplayout_t>* bunch)
#else
                 PartBunchBase* bunch)
#endif
    : AmrMesh(rb, max_level_in, n_cell_in, coord),
#ifdef IPPL_AMR
      bunch_m(bunch),
#else
      bunch_m(dynamic_cast<AmrPartBunch*>(bunch)),
#endif
      tagging_m(kChargeDensity),
      scaling_m(0.75),
      nCharge_m(1.0e-15),
      minNumPart_m(1),
      maxNumPart_m(1),
      meshScaling_m(Vector_t(1.0, 1.0, 1.0))
{
    initBaseLevel();
    
    nChargePerCell_m.resize(max_level_in + 1);
    nChargePerCell_m[0] = std::unique_ptr<amrex::MultiFab>(new amrex::MultiFab(this->boxArray(0), this->DistributionMap(0), 1, 1));
    nChargePerCell_m[0]->setVal(0.0, 1);
}

AmrOpal::AmrOpal(const amrex::RealBox* rb, int max_level_in, const amrex::Vector<int>& n_cell_in, int coord)
    : AmrMesh(rb, max_level_in, n_cell_in, coord),
      tagging_m(kChargeDensity),
      scaling_m(0.75),
      nCharge_m(1.0e-15),
      minNumPart_m(1),
      maxNumPart_m(1),
      meshScaling_m(Vector_t(1.0, 1.0, 1.0))
{
    finest_level = 0;
    
    const amrex::BoxArray& ba = MakeBaseGrids();
    
    amrex::DistributionMapping dm(ba, amrex::ParallelDescriptor::NProcs());
    
    nChargePerCell_m.resize(max_level_in + 1);
    
    MakeNewLevel(0, 0.0, ba, dm);
}

AmrOpal::AmrOpal(const amrex::RealBox* rb, int max_level_in, const amrex::Vector<int>& n_cell_in, int coord,
                 const std::vector<int>& refratio)
    : AmrMesh(rb, max_level_in, n_cell_in, coord, refratio),
      tagging_m(kChargeDensity),
      scaling_m(0.75),
      nCharge_m(1.0e-15),
      minNumPart_m(1),
      maxNumPart_m(1),
      meshScaling_m(Vector_t(1.0, 1.0, 1.0))
{
    finest_level = 0;
    
    const amrex::BoxArray& ba = MakeBaseGrids();
    
    amrex::DistributionMapping dm(ba, amrex::ParallelDescriptor::NProcs());
    
    nChargePerCell_m.resize(max_level_in + 1);
    
    MakeNewLevel(0, 0.0, ba, dm);
}


AmrOpal::~AmrOpal() { }



/* init the base level by using the distributionmap and BoxArray of the
 * particles (called in constructor of AmrOpal
 */
void AmrOpal::initBaseLevel() {
    finest_level = 0; // AmrMesh protected member variable
#ifdef IPPL_AMR
    amrplayout_t* PLayout = &bunch_m->getLayout();
    const amrex::BoxArray& ba = PLayout->ParticleBoxArray(0);
    const amrex::DistributionMapping& dm = PLayout->ParticleDistributionMap(0 /*level*/);
#else
    const amrex::BoxArray& ba = bunch_m->ParticleBoxArray(0 /*level*/);
    const amrex::DistributionMapping& dm = bunch_m->ParticleDistributionMap(0 /*level*/);
#endif
    
    SetBoxArray(0 /*level*/, ba);
    SetDistributionMap(0 /*level*/, dm);
}

// void AmrOpal::initFineLevels() { }


void AmrOpal::writePlotFileYt(std::string filename, int step) {
    mfs_mt chargeOnGrid;
    chargeOnGrid.resize(finest_level + 1);
    
    for (int i = 0; i < finest_level; ++i) {
        chargeOnGrid[i].reset(new amrex::MultiFab(this->boxArray(i), this->DistributionMap(i), 1, 0));
        chargeOnGrid[i]->setVal(0.0);
    }
        
#if defined IPPL_AMR
        bunch_m->AssignDensity(bunch_m->qm, false, chargeOnGrid, 0, finest_level);
#elif !defined IPPL_AMR
        bunch_m->AssignDensity(0, false, chargeOnGrid, 0, 1, finest_level);
#endif
    
    std::string dir = filename;
    int nLevels = chargeOnGrid.size();
    //
    // Only let 64 CPUs be writing at any one time.
    //
    amrex::VisMF::SetNOutFiles(64);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (amrex::ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(dir, 0755))
            amrex::CreateDirectoryFailed(dir);
    //
    // Force other processors to wait till directory is built.
    //
    amrex::ParallelDescriptor::Barrier();

    std::string HeaderFileName = dir + "/Header";

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int nData = chargeOnGrid[0]->nComp();
    
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);
        HeaderFile << "HyperCLaw-V1.1\n"; //"NavierStokes-V1.1\n"

        HeaderFile << nData << '\n';

	// variable names
        for (int ivar = 1; ivar <= chargeOnGrid[0]->nComp(); ivar++) {
          HeaderFile << "rho\n";
        }
        
        
        // dimensionality
        HeaderFile << AMREX_SPACEDIM << '\n';
        
        // time
        HeaderFile << 0.0 << '\n';
        HeaderFile << nLevels - 1 << std::endl; // maximum level number (0=single level)
        
        // physical domain
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            HeaderFile << geom[0].ProbLo(i) << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            HeaderFile << geom[0].ProbHi(i) << ' ';
        HeaderFile << std::endl;
        
        // reference ratio
        for (int i = 0; i < nLevels - 1; ++i)
            HeaderFile << 2 << ' ';
        HeaderFile << std::endl;
        
        // geometry domain for all levels
        for (int i = 0; i < nLevels; ++i)
            HeaderFile << geom[i].Domain() << ' ';
        HeaderFile << std::endl;
        
        // number of time steps
        for (int i = 0; i < nLevels - 1; ++i)
            HeaderFile << 0 << " ";
        HeaderFile << std::endl;
        
        // cell sizes for all level
        for (int i = 0; i < nLevels; ++i) {
            for (int k = 0; k < AMREX_SPACEDIM; k++)
                HeaderFile << geom[i].CellSize()[k] << ' ';
            HeaderFile << '\n';
        }
        
        // coordinate system
        HeaderFile << geom[0].Coord() << '\n';
        HeaderFile << "0\n"; // write boundary data
    }
    
    for (int lev = 0; lev < nLevels; ++lev) {
        // Build the directory to hold the MultiFab at this level.
        // The name is relative to the directory containing the Header file.
        //
        static const std::string BaseName = "/Cell";
    
        std::string Level = amrex::Concatenate("Level_", lev, 1);
        //
        // Now for the full pathname of that directory.
        //
        std::string FullPath = dir;
        if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
            FullPath += '/';
        FullPath += Level;
        //
        // Only the I/O processor makes the directory if it doesn't already exist.
        //
        if (amrex::ParallelDescriptor::IOProcessor())
            if (!amrex::UtilCreateDirectory(FullPath, 0755))
                amrex::CreateDirectoryFailed(FullPath);
        //
        // Force other processors to wait till directory is built.
        //
        amrex::ParallelDescriptor::Barrier();
        
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            HeaderFile << lev << ' ' << chargeOnGrid[lev]->boxArray().size() << ' ' << 0 << '\n';
            HeaderFile << 0 << '\n';    // # time steps at this level
    
            for (int i = 0; i < chargeOnGrid[lev]->boxArray().size(); ++i)
            {
                amrex::RealBox loc = amrex::RealBox(chargeOnGrid[lev]->boxArray()[i],geom[lev].CellSize(),geom[lev].ProbLo());
                for (int n = 0; n < AMREX_SPACEDIM; n++)
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
            }
    
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            HeaderFile << PathNameInHeader << '\n';
        }
        
        amrex::MultiFab data(chargeOnGrid[lev]->boxArray(), chargeOnGrid[lev]->DistributionMap(), nData, 0);
        // dst, src, srccomp, dstcomp, numcomp, nghost
        /*
        * srccomp: the component to copy
        * dstcmop: the component where to copy
        * numcomp: how many components to copy
        */
        amrex::MultiFab::Copy(data, *chargeOnGrid[lev],    0, 0, 1, 0);
        //
        // Use the Full pathname when naming the MultiFab.
        //
        std::string TheFullPath = FullPath;
        TheFullPath += BaseName;
    
        amrex::VisMF::Write(data,TheFullPath);
    }
    
    for (int i = 0; i <= finest_level; ++i) {
        chargeOnGrid[i].reset(nullptr);
    }
}

void AmrOpal::writePlotFile(std::string filename, int step) {
    mfs_mt chargeOnGrid;
    chargeOnGrid.resize(finest_level + 1);
    
    for (int i = 0; i < finest_level; ++i) {
        chargeOnGrid[i].reset(new amrex::MultiFab(this->boxArray(i), this->DistributionMap(i), 1, 0));
        chargeOnGrid[i]->setVal(0.0);
    }
    
    amrex::Vector<std::string> varnames(1, "rho");
    
    amrex::Vector<const amrex::MultiFab*> tmp(finest_level + 1);
    for (/*unsigned*/ int i = 0; i < finest_level + 1; ++i) {
        tmp[i] = chargeOnGrid[i].get();
    }
    
    const auto& mf = tmp;
    
    amrex::Vector<int> istep(finest_level+1, step);
    
    amrex::WriteMultiLevelPlotfile(filename, finest_level + 1, mf, varnames,
                                   Geom(), 0.0, istep, refRatio());
    
    for (int i = 0; i <= finest_level; ++i) {
        chargeOnGrid[i].reset(nullptr);
    }
}


void AmrOpal::ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) {
    switch ( tagging_m ) {
        case kChargeDensity:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
        case kPotentialStrength:
            tagForPotentialStrength_m(lev, tags, time, ngrow);
            break;
        case kEfieldStrength:
            tagForEfieldStrength_m(lev, tags, time, ngrow);
            break;
        case kMomentum:
            tagForMomentum_m(lev, tags, time, ngrow);
            break;
        case kMaxNumParticles:
            tagForMaxNumParticles_m(lev, tags, time, ngrow);
            break;
        case kMinNumParticles:
            tagForMinNumParticles_m(lev, tags, time, ngrow);
            break;
        case kCenteredRegion:
            tagForCenteredRegion_m(lev, tags, time, ngrow);
            break;
        default:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
    }
}


void AmrOpal::MakeNewLevelFromScratch(int lev, amrex::Real time, const amrex::BoxArray& ba,
                                      const amrex::DistributionMapping& dm)
{
    amrex::Abort("How did we get her!");
}


void AmrOpal::MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray& ba,
                                     const amrex::DistributionMapping& dm)
{
    amrex::Abort("How did we get her!");
}


void
AmrOpal::regrid (int lbase, amrex::Real time)
{
    int new_finest = 0;
    amrex::Vector<amrex::BoxArray> new_grids(finest_level+2);
    
    MakeNewGrids(lbase, time, new_finest, new_grids);

    BL_ASSERT(new_finest <= finest_level+1);
    
    for (int lev = lbase+1; lev <= new_finest; ++lev)
    {
        if (lev <= finest_level) // an old level
        {
            if (new_grids[lev] != grids[lev]) // otherwise nothing
            {
                amrex::DistributionMapping new_dmap(new_grids[lev], amrex::ParallelDescriptor::NProcs());
                RemakeLevel(lev, time, new_grids[lev], new_dmap);
                
                /*
                 * particles need to know the BoxArray
                 * and DistributionMapping
                 */
#ifdef IPPL_AMR
                amrplayout_t* PLayout = &bunch_m->getLayout();
                PLayout->SetParticleBoxArray(lev, new_grids[lev]);
                PLayout->SetParticleDistributionMap(lev, new_dmap);
#else
                bunch_m->SetParticleBoxArray(lev, new_grids[lev]);
                bunch_m->SetParticleDistributionMap(lev, new_dmap);
#endif
	    }
	}
	else  // a new level
	{
	    amrex::DistributionMapping new_dmap(new_grids[lev], amrex::ParallelDescriptor::NProcs());
	    MakeNewLevel(lev, time, new_grids[lev], new_dmap);
            
            /*
             * particles need to know the BoxArray
             * and DistributionMapping
             */
#ifdef IPPL_AMR
            amrplayout_t* PLayout = &bunch_m->getLayout();
            PLayout->SetParticleBoxArray(lev, new_grids[lev]);
            PLayout->SetParticleDistributionMap(lev, new_dmap);
#else
            bunch_m->SetParticleBoxArray(lev, new_grids[lev]);
            bunch_m->SetParticleDistributionMap(lev, new_dmap);
#endif
        }
    }
    
    for (int lev = new_finest+1; lev <= finest_level; ++lev) {
        ClearLevel(lev);
        ClearBoxArray(lev);
        ClearDistributionMap(lev);
    }
    
    finest_level = new_finest;
    
    // update to multilevel
#ifdef IPPL_AMR
    bunch_m->update(lbase, finest_level);
#else
    bunch_m->myUpdate();
#endif
}


void
AmrOpal::RemakeLevel (int lev, amrex::Real time,
                      const amrex::BoxArray& new_grids, const amrex::DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m[lev].reset(new amrex::MultiFab(new_grids, new_dmap, 1, 1));
}

void
AmrOpal::MakeNewLevel (int lev, amrex::Real time,
                       const amrex::BoxArray& new_grids, const amrex::DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m[lev] = std::unique_ptr<amrex::MultiFab>(new amrex::MultiFab(new_grids, dmap[lev], 1, 1));
}

void AmrOpal::ClearLevel(int lev) {
    
    nChargePerCell_m[lev].reset(nullptr);
    
    ClearBoxArray(lev);
    ClearDistributionMap(lev);
}


void AmrOpal::scatter_m(int lev) {
    for (int i = lev; i <= finest_level; ++i) {
        nChargePerCell_m[i]->setVal(0.0, 1);
    }
    
    mfs_mt partMF(finest_level + 1);
    for (int i = lev; i <= finest_level; i++) {
        const amrex::BoxArray& ba = boxArray()[i];
        const amrex::DistributionMapping& dmap = DistributionMap(i);
        partMF[i].reset(new amrex::MultiFab(ba, dmap, 1, 2));
        partMF[i]->setVal(0.0, 2);
   }
    
    bunch_m->AssignDensityFort(bunch_m->qm, partMF, lev, 1, finest_level);
    
    for (int i = lev; i <= finest_level; ++i) {
        amrex::MultiFab::Copy(*nChargePerCell_m[i], *partMF[i], 0, 0, 1, 0);
    }
}


// use scale instead of time to scale charge due to mapping of particles
void AmrOpal::tagForChargeDensity_m(int lev, amrex::TagBoxArray& tags, amrex::Real scale/*time*/, int ngrow) {
    
    scatter_m(lev);
    
    for (int i = lev; i <= finest_level; ++i) {
        nChargePerCell_m[i]->mult(scale, 1);
    }
    
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;

    const amrex::Real* dx      = geom[lev].CellSize();
    const amrex::Real* prob_lo = geom[lev].ProbLo();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        amrex::Vector<int>  itags;
        for (amrex::MFIter mfi(*nChargePerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
            const amrex::Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            amrex::TagBox&     tagfab  = tags[mfi];
            
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
            tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D((*nChargePerCell_m[lev])[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &scale, &nCharge_m);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
}

void AmrOpal::tagForPotentialStrength_m(int lev, amrex::TagBoxArray& tags, amrex::Real scale/*time*/, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the potential is higher than 75 percent of the maximum potential
     * value of this level.
     */
    
    // 1. Assign charge onto grid of level lev
    int lbase   = 0;
    int lfinest = finest_level;
    
    int nlevels = lfinest - lbase + 1;
    
    mfs_mt nPartPerCell(nlevels);
    
    for ( int i = 0; i < nlevels; ++ i) {
        nPartPerCell[i] = std::unique_ptr<amrex::MultiFab>(new amrex::MultiFab(this->boxArray(i),
                                                             this->DistributionMap(i),
                                                             1, 2)
                                                 );
        nPartPerCell[i]->setVal(0.0, 2);
    }
    
    #ifdef IPPL_AMR
        bunch_m->AssignDensityFort(bunch_m->qm, nPartPerCell, lbase, 1, lfinest);
//         bunch_m->AssignCellDensitySingleLevelFort(bunch_m->qm, *nPartPerCell[lbase], lev);
//         bunch_m->AssignDensitySingleLevel(bunch_m->qm, *nPartPerCell[lbase], lev);
    #else
        bunch_m->AssignDensitySingleLevel(0, *nPartPerCell[lbase], lev);
    #endif
    
    for ( int i = 0; i < nlevels; ++ i)
        nPartPerCell[i]->mult(scale, 0, 1);
    
    // 2. Solve Poisson's equation on level lev
    
    amrex::Real offset = 0.0;
    
    if ( geom[0].isAllPeriodic() ) {
        for (std::size_t i = 0; i < bunch_m->getLocalNum(); ++i) {
             #ifdef IPPL_AMR
                offset += bunch_m->qm[i];
            #else
                offset += ((PartBunchBase*)bunch_m)->getQM(i);
            #endif
        }
        offset /= geom[0].ProbSize();
    }
    
    Solver::container_t phi(nlevels);
    Solver::container_t grad_phi(nlevels);
    for ( int i = 0; i < nlevels; ++i ) {
        //                                                                # component # ghost cells                                                                                                                                          
        phi[i] = std::unique_ptr<amrex::MultiFab>(
            new amrex::MultiFab(this->boxArray(i), this->DistributionMap(i), 1          , 1));
        grad_phi[i] = std::unique_ptr<amrex::MultiFab>(
            new amrex::MultiFab(this->boxArray(i), this->DistributionMap(i), AMREX_SPACEDIM, 1));
        phi[i]->setVal(0.0, 1);
        grad_phi[i]->setVal(0.0, 1);
    }
    
    Solver sol;
    sol.solve_for_accel(nPartPerCell,
                        phi,
                        grad_phi,
                        geom,
                        lbase,
                        lfinest,
                        offset,
                        false,
                        false); // we need no gradient
    
    // 3. Tag cells for refinement
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;

    const amrex::Real* dx      = geom[lev].CellSize();
    const amrex::Real* prob_lo = geom[lev].ProbLo();
    
    amrex::Real minp = 0.0;
    amrex::Real maxp = 0.0;
    maxp = scaling_m * phi[lev]->max(0);
    minp = scaling_m * phi[lev]->min(0);
    amrex::Real value = std::max( std::abs(minp), std::abs(maxp) );
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        amrex::Vector<int>  itags;
        for (amrex::MFIter mfi(*phi[lev],false/*true*/); mfi.isValid(); ++mfi) {
            const amrex::Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            amrex::TagBox&     tagfab  = tags[mfi];
            
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
            tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            tag_potential_strength(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D((*phi[lev])[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &scale, &value);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
    
    // releases memory itself
    for (int i = 0; i < nlevels; ++i)
        nPartPerCell[i].reset(nullptr);
}

void AmrOpal::tagForEfieldStrength_m(int lev, amrex::TagBoxArray& tags, amrex::Real scale/*time*/, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the potential is higher than 75 percent of the maximum potential
     * value of this level.
     */
    
    // 1. Assign charge onto grid of level lev
    int lbase   = 0;
    int lfinest = finest_level;
    
    int nlevels = lfinest - lbase + 1;
    
    mfs_mt nPartPerCell(nlevels);
    
    for ( int i = 0; i < nlevels; ++ i) {
        nPartPerCell[i] = std::unique_ptr<amrex::MultiFab>(new amrex::MultiFab(this->boxArray(i),
                                                             this->DistributionMap(i),
                                                             1, 2)
                                                 );
        nPartPerCell[i]->setVal(0.0, 2);
    }
    
    #ifdef IPPL_AMR
        bunch_m->AssignDensityFort(bunch_m->qm, nPartPerCell, lbase, 1, lfinest);
//         bunch_m->AssignCellDensitySingleLevelFort(bunch_m->qm, *nPartPerCell[lbase], lev);
//         bunch_m->AssignDensitySingleLevel(bunch_m->qm, *nPartPerCell[lbase], lev);
    #else
        bunch_m->AssignDensitySingleLevel(0, *nPartPerCell[lbase], lev);
    #endif
    
    for ( int i = 0; i < nlevels; ++ i)
        nPartPerCell[i]->mult(scale, 0, 1);
    
    // 2. Solve Poisson's equation on level lev
    amrex::Real offset = 0.0;
    
    if ( geom[0].isAllPeriodic() ) {
        for (std::size_t i = 0; i < bunch_m->getLocalNum(); ++i) {
             #ifdef IPPL_AMR
                offset += bunch_m->qm[i];
            #else
                offset += ((PartBunchBase*)bunch_m)->getQM(i);
            #endif
        }
        offset /= geom[0].ProbSize();
    }
    
    Solver::container_t phi(nlevels);
    Solver::container_t grad_phi(nlevels);
    for ( int i = 0; i < nlevels; ++i ) {
        //                                                                # component # ghost cells                                                                                                                                          
        phi[i] = std::unique_ptr<amrex::MultiFab>(
            new amrex::MultiFab(this->boxArray(i), this->DistributionMap(i), 1          , 1));
        grad_phi[i] = std::unique_ptr<amrex::MultiFab>(
            new amrex::MultiFab(this->boxArray(i), this->DistributionMap(i), AMREX_SPACEDIM, 1));
        phi[i]->setVal(0.0, 1);
        grad_phi[i]->setVal(0.0, 1);
    }
    
    Solver sol;
    sol.solve_for_accel(nPartPerCell,
                        phi,
                        grad_phi,
                        geom,
                        lbase,
                        lfinest,
                        offset,
                        false);
    
    // 3. Tag cells for refinement
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;
    
    amrex::Real min[3] = {0.0, 0.0, 0.0};
    amrex::Real max[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        max[i] = scaling_m * grad_phi[lev]->max(i);
        min[i] = scaling_m * grad_phi[lev]->min(i);
    }
    amrex::Real efield[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)
        efield[i] = std::max( std::abs(min[i]), std::abs(max[i]) );
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
//         Array<int>  itags;
        for (amrex::MFIter mfi(*grad_phi[lev],false/*true*/); mfi.isValid(); ++mfi) {
            const amrex::Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            amrex::TagBox&     tagfab  = tags[mfi];
            amrex::FArrayBox&  fab     = (*grad_phi[lev])[mfi];
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
//             tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
//             int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        amrex::IntVect iv(i,j,k);
                        if (std::abs(fab(iv, 0)) >= efield[0])
                            tagfab(iv) = tagval;
                        else if (std::abs(fab(iv, 1)) >= efield[1])
                            tagfab(iv) = tagval;
                        else if (std::abs(fab(iv, 2)) >= efield[2])
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
                            
                        
            //
            // Now update the tags in the TagBox.
            //
//             tagfab.tags_and_untags(itags, tilebx);
        }
    }
    
    
    // releases memory itself
    for (int i = 0; i < nlevels; ++i)
        nPartPerCell[i].reset(nullptr);
}


void AmrOpal::tagForMomentum_m(int lev, amrex::TagBoxArray& tags, amrex::Real scale/*time*/, int ngrow) {
#ifdef IPPL_AMR
    amrplayout_t* PLayout = &bunch_m->getLayout();
    const amrex::BoxArray& ba = PLayout->ParticleBoxArray(lev);
    const amrex::DistributionMapping& dm = PLayout->ParticleDistributionMap(lev);
    
    
    const auto& LocalNumPerLevel = bunch_m->getLocalNumPerLevel();
    
    size_t lBegin = LocalNumPerLevel.begin(lev);
    size_t lEnd   = LocalNumPerLevel.end(lev);
    
    Vector_t pmax = Vector_t(0.0, 0.0, 0.0);
    for (size_t i = lBegin; i < lEnd; ++i) {
        const Vector_t& tmp = bunch_m->P[i];
        pmax = ( dot(tmp, tmp) > dot(pmax, pmax) ) ? tmp : pmax;
    }
    
    double momentum2 = scaling_m * dot(pmax, pmax);
    
    std::map<amrex::IntVect, bool> cells;
    for (size_t i = lBegin; i < lEnd; ++i) {
        if ( dot(bunch_m->P[i], bunch_m->P[i]) >= momentum2 ) {
            amrex::IntVect iv = PLayout->Index(bunch_m->R[i], lev);
            cells[iv] = true;
        }
    }
        
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
//         Array<int>  itags;
        for (amrex::MFIter mfi(*nChargePerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
            const amrex::Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            amrex::TagBox&     tagfab  = tags[mfi];
            
            // data pointer and index space
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        amrex::IntVect iv(i, j, k);
                        if ( cells.find(iv) != cells.end() )
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
            //
            // Now update the tags in the TagBox.
            //
//             tagfab.tags_and_untags(itags, tilebx);
        }
    }
#endif
}


void AmrOpal::tagForMaxNumParticles_m(int lev, amrex::TagBoxArray& tags, amrex::Real scale/*time*/, int ngrow) {
#ifdef IPPL_AMR
    amrplayout_t* PLayout = &bunch_m->getLayout();
    const amrex::BoxArray& ba = PLayout->ParticleBoxArray(lev);
    const amrex::DistributionMapping& dm = PLayout->ParticleDistributionMap(lev);
    
    
    const auto& LocalNumPerLevel = bunch_m->getLocalNumPerLevel();
    
    size_t lBegin = LocalNumPerLevel.begin(lev);
    size_t lEnd   = LocalNumPerLevel.end(lev);
    
    // count particles per cell
    std::map<amrex::IntVect, size_t> cells;
    for (size_t i = lBegin; i < lEnd; ++i) {
        amrex::IntVect iv = PLayout->Index(bunch_m->R[i], lev);
        ++cells[iv];
    }
        
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
//         Array<int>  itags;
        for (amrex::MFIter mfi(*nChargePerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
            const amrex::Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            amrex::TagBox&     tagfab  = tags[mfi];
//             tagfab.get_itags(itags, tilebx);

            // data pointer and index space
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        amrex::IntVect iv(i, j, k);
                        if ( cells.find(iv) != cells.end() && cells[iv] <= maxNumPart_m )
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
            //
            // Now update the tags in the TagBox.
            //
//             tagfab.tags_and_untags(itags, tilebx);
        }
    }
#endif
}


void AmrOpal::tagForMinNumParticles_m(int lev, amrex::TagBoxArray& tags, amrex::Real scale/*time*/, int ngrow) {
#ifdef IPPL_AMR
    amrplayout_t* PLayout = &bunch_m->getLayout();
    const amrex::BoxArray& ba = PLayout->ParticleBoxArray(lev);
    const amrex::DistributionMapping& dm = PLayout->ParticleDistributionMap(lev);
    
    // count particles per cell
    std::map<amrex::IntVect, size_t> cells;
    for (size_t i = 0; i < bunch_m->getLocalNum(); ++i) {
        amrex::IntVect iv = PLayout->Index(bunch_m->R[i], lev);
        ++cells[iv];
    }
        
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
//         Array<int>  itags;
        for (amrex::MFIter mfi(*nChargePerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
            const amrex::Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            amrex::TagBox&     tagfab  = tags[mfi];
            
            // data pointer and index space
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        amrex::IntVect iv(i, j, k);
                        if ( cells.find(iv) != cells.end() && cells[iv] >= minNumPart_m )
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
            //
            // Now update the tags in the TagBox.
            //
//             tagfab.tags_and_untags(itags, tilebx);
        }
    }
#endif
}


void AmrOpal::tagForCenteredRegion_m(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) {
    
    for (int i = lev; i <= finest_level; ++i) {
        nChargePerCell_m[i]->setVal(0.0, 1);
        #ifdef IPPL_AMR
            bunch_m->AssignDensitySingleLevel(bunch_m->qm, *nChargePerCell_m[i], i);
        #else
            bunch_m->AssignDensitySingleLevel(0, *nChargePerCell_m[i], i);
        #endif
    }
    
    for (int i = finest_level-1; i >= lev; --i) {
        amrex::MultiFab tmp(nChargePerCell_m[i]->boxArray(), nChargePerCell_m[i]->DistributionMap(), 1, 0);
        tmp.setVal(0.0);
        amrex::average_down(*nChargePerCell_m[i+1], tmp, 0, 1, refRatio(i));
        amrex::MultiFab::Add(*nChargePerCell_m[i], tmp, 0, 0, 1, 0);
    }
    
    
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;

    const amrex::Real* dx      = geom[lev].CellSize();
    const amrex::Real* prob_lo = geom[lev].ProbLo();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        amrex::Vector<int>  itags;
        for (amrex::MFIter mfi(*nChargePerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
            const amrex::Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            amrex::TagBox&     tagfab  = tags[mfi];
            
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
            tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            centered_region(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D((*nChargePerCell_m[lev])[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &nCharge_m);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
}
