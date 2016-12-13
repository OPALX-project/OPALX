#include "AmrOpal.h"
#include "AmrOpal_F.h"
#include <PlotFileUtil.H>
// #include <MultiFabUtil.H>


// AmrOpal::AmrOpal() { }
AmrOpal::AmrOpal(const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord, PartBunchBase* bunch)
    : AmrCore(rb, max_level_in, n_cell_in, coord),
    bunch_m(dynamic_cast<AmrPartBunch*>(bunch))
{
    initBaseLevel();
    nPartPerCell_m.resize(max_level_in + 1);//, PArrayManage);
//     nPartPerCell_m.set(0, new MultiFab(this->boxArray(0), 1, 1));
    
    nPartPerCell_m[0] = std::unique_ptr<MultiFab>(new MultiFab(this->boxArray(0), 1, 1, this->DistributionMap(0)));
    nPartPerCell_m[0]->setVal(0.0);
    
    chargeOnGrid_m.resize(max_level_in + 1);
    chargeOnGrid_m[0] = std::unique_ptr<MultiFab>(new MultiFab(this->boxArray(0), 1, 0, this->DistributionMap(0)));
    
//     bunch_m->AssignDensitySingleLevel(0, *nPartPerCell_m[0], 0);
//     bunch_m->myUpdate();
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


void AmrOpal::writePlotFileYt(std::string filename, int step) {
    
    std::string dir = filename;
    int nLevels = chargeOnGrid_m.size();
    //
    // Only let 64 CPUs be writing at any one time.
    //
    VisMF::SetNOutFiles(64);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(dir, 0755))
            BoxLib::CreateDirectoryFailed(dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    std::string HeaderFileName = dir + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int nData = chargeOnGrid_m[0]->nComp();
    
    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);
        HeaderFile << "HyperCLaw-V1.1\n"; //"NavierStokes-V1.1\n"

        HeaderFile << nData << '\n';

	// variable names
        for (int ivar = 1; ivar <= chargeOnGrid_m[0]->nComp(); ivar++) {
          HeaderFile << "rho\n";
        }
        
        
        // dimensionality
        HeaderFile << BL_SPACEDIM << '\n';
        
        // time
        HeaderFile << 0.0 << '\n';
        HeaderFile << nLevels - 1 << std::endl; // maximum level number (0=single level)
        
        // physical domain
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom[0].ProbLo(i) << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
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
            for (int k = 0; k < BL_SPACEDIM; k++)
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
    
        std::string Level = BoxLib::Concatenate("Level_", lev, 1);
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
        if (ParallelDescriptor::IOProcessor())
            if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
                BoxLib::CreateDirectoryFailed(FullPath);
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        
        if (ParallelDescriptor::IOProcessor())
        {
            HeaderFile << lev << ' ' << chargeOnGrid_m[lev]->boxArray().size() << ' ' << 0 << '\n';
            HeaderFile << 0 << '\n';    // # time steps at this level
    
            for (int i = 0; i < chargeOnGrid_m[lev]->boxArray().size(); ++i)
            {
                RealBox loc = RealBox(chargeOnGrid_m[lev]->boxArray()[i],geom[lev].CellSize(),geom[lev].ProbLo());
                for (int n = 0; n < BL_SPACEDIM; n++)
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
            }
    
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            HeaderFile << PathNameInHeader << '\n';
        }
        
        MultiFab data(chargeOnGrid_m[lev]->boxArray(), nData, 0);
        
        // dst, src, srccomp, dstcomp, numcomp, nghost
        /*
        * srccomp: the component to copy
        * dstcmop: the component where to copy
        * numcomp: how many components to copy
        */
//         data.copy(*nPartPerCell_m[lev],0,0,0);
        MultiFab::Copy(data, *chargeOnGrid_m[lev],    0, 0, 1, 0);
        
        //
        // Use the Full pathname when naming the MultiFab.
        //
        std::string TheFullPath = FullPath;
        TheFullPath += BaseName;
    
        VisMF::Write(data,TheFullPath);
    }
    
}

void AmrOpal::writePlotFile(std::string filename, int step) {
    Array<std::string> varnames(1, "rho");
    
    Array<const MultiFab*> tmp(finest_level + 1/*nPartPerCell_m.size()*/);
    for (/*unsigned*/ int i = 0; i < finest_level + 1/*nPartPerCell_m.size()*/; ++i) {
        tmp[i] = chargeOnGrid_m[i].get();
    }
    
    const auto& mf = tmp;
    
    Array<int> istep(finest_level+1, step);
    
    BoxLib::WriteMultiLevelPlotfile(filename, finest_level + 1, mf, varnames,
                                    Geom(), 0.0, istep, refRatio());
}


void AmrOpal::ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ngrow*/) {
    
    for (int i = lev; i <= finest_level; ++i) {
        nPartPerCell_m[i]->setVal(0.0);
        bunch_m->AssignDensitySingleLevel(0, *nPartPerCell_m[i], i);
    }

    for (int i = finest_level-1; i >= lev; --i) {
        MultiFab tmp(nPartPerCell_m[i]->boxArray(), 1, 0, nPartPerCell_m[i]->DistributionMap());
        tmp.setVal(0.0);
        BoxLib::average_down(*nPartPerCell_m[i+1], tmp, 0, 1, refRatio(i));
        MultiFab::Add(*nPartPerCell_m[i], tmp, 0, 0, 1, 0);
    }
    
    
//     std::cout << lev << " " << nPartPerCell_m[lev]->min(0) << " " << nPartPerCell_m[lev]->max(0) << std::endl;
    
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real nPart = 1.0e-15;
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        
        for (MFIter mfi(*nPartPerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
            const Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            TagBox&     tagfab  = tags[mfi];
            
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
            tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D((*nPartPerCell_m[lev])[mfi]),
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
    
//     if (new_finest > finest_level)
        finest_level = new_finest;
    
    // update to multilevel
    bunch_m->myUpdate();
    
    
//     for (int i = 0; i <= finest_level; ++i) {
//         nPartPerCell_m[i]->setVal(0.0);
//         bunch_m->AssignDensitySingleLevel(0, *nPartPerCell_m[i], i);
//     }
    
    
//     bunch_m->AssignDensity(0, false, nPartPerCell_m, 0, 1, finest_level);
    
//     for (int i = finest_level-1; i >= 0; --i) {
//         MultiFab tmp(nPartPerCell_m[i]->boxArray(), 1, 0, nPartPerCell_m[i]->DistributionMap());
//         tmp.setVal(0.0);
//         BoxLib::average_down(*nPartPerCell_m[i+1], tmp, 0, 1, refRatio(i));
//         MultiFab::Add(*nPartPerCell_m[i], tmp, 0, 0, 1, 0);
//     }
}


void
AmrOpal::RemakeLevel (int lev, Real time,
		     const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
    
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nPartPerCell_m[lev].reset(new MultiFab(new_grids, 1, 1, new_dmap));
    
    chargeOnGrid_m[lev].reset(new MultiFab(new_grids, 1, 0, new_dmap));
}

void
AmrOpal::MakeNewLevel (int lev, Real time,
		      const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nPartPerCell_m[lev] = std::unique_ptr<MultiFab>(new MultiFab(new_grids, 1, 1, dmap[lev]));
    
    chargeOnGrid_m[lev] = std::unique_ptr<MultiFab>(new MultiFab(new_grids, 1, 0, dmap[lev]));
}

void AmrOpal::ClearLevel(int lev) {
    
    nPartPerCell_m[lev].reset(nullptr);
    chargeOnGrid_m[lev].reset(nullptr);
    ClearBoxArray(lev);
    ClearDistributionMap(lev);
}
