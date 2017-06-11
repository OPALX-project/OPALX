#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>

#define UNIQUE_PTR 1

#ifndef UNIQUE_PTR
    #include <PArray.H>
#endif

/*!
 * @file writePlotFile.H
 * @author Matthias Frey
 * @date October, 2016, LBNL
 * @details This functions writes the density, potential
 * and the electric field components to a AMReX plotfile
 * that can be visualized using yt.
 * @brief Write a AMReX plotfile
 */

/*!
 * @param dir where to store
 * @param rho is the density on all levels
 * @param phi is the potential on all levels
 * @param efield are the electric field components on all levels
 * @param refRatio are the refinement ratios among the levels
 * @param geom are the geometries of all levels
 * @param time specifies the step.
 */
inline void
#ifdef UNIQUE_PTR
writePlotFile (const std::string& dir,
               const std::vector<std::unique_ptr<MultiFab> >& rho,
               const std::vector<std::unique_ptr<MultiFab> >& phi,
               const std::vector<std::unique_ptr<MultiFab> >& efield,
               const Array<int>& refRatio,
               const Array<Geometry>& geom,
	       const Real&        time)
#else
writePlotFile (const std::string& dir,
               const PArray<MultiFab>& rho,
               const PArray<MultiFab>& phi,
               const PArray<MultiFab>& efield,
               const Array<int>& refRatio,
               const Array<Geometry>& geom,
	       const Real&        time)
#endif
{
    int nLevels = rho.size();
    //
    // Only let 64 CPUs be writing at any one time.
    //
    VisMF::SetNOutFiles(64);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(dir, 0755))
            amrex::CreateDirectoryFailed(dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    std::string HeaderFileName = dir + "/OpalAmrHeader";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

#ifdef UNIQUE_PTR
    int nData = rho[0]->nComp() + phi[0]->nComp() + efield[0]->nComp();
#else
    int nData = rho[0].nComp() + phi[0].nComp() + efield[0].nComp();
#endif
    
    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);
        HeaderFile << "OpalAmr-V1.0\n";

        HeaderFile << nData << '\n';

	// variable names
#ifdef UNIQUE_PTR
        for (int ivar = 1; ivar <= rho[0]->nComp(); ivar++) {
#else
        for (int ivar = 1; ivar <= rho[0].nComp(); ivar++) {
#endif
          HeaderFile << "rho\n";
        }
        
#ifdef UNIQUE_PTR
        for (int ivar = 1; ivar <= phi[0]->nComp(); ivar++) {
#else
        for (int ivar = 1; ivar <= phi[0].nComp(); ivar++) {
#endif
          HeaderFile << "potential\n";
        }
        
//         for (int ivar = 1; ivar <= efield.nComp(); ivar++) {
          HeaderFile << "Ex\nEy\nEz\n";
//         }
        
        
        // dimensionality
        HeaderFile << BL_SPACEDIM << '\n';
        
        // time
        HeaderFile << time << '\n';
        HeaderFile << nLevels - 1 << std::endl; // maximum level number (0=single level)
        
        // physical domain
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom[0].ProbLo(i) << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom[0].ProbHi(i) << ' ';
        HeaderFile << std::endl;
        
        // reference ratio
        for (std::size_t i = 0; i < refRatio.size(); ++i)
            HeaderFile << refRatio[i] << ' ';
        HeaderFile << std::endl;
        
        // geometry domain for all levels
        for (int i = 0; i < nLevels; ++i)
            HeaderFile << geom[i].Domain() << ' ';
        HeaderFile << std::endl;
        
        // number of time steps
        HeaderFile << 0 << " " << std::endl;
        
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
        if (ParallelDescriptor::IOProcessor())
            if (!amrex::UtilCreateDirectory(FullPath, 0755))
                amrex::CreateDirectoryFailed(FullPath);
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        
        if (ParallelDescriptor::IOProcessor())
        {
#ifdef UNIQUE_PTR
            HeaderFile << lev << ' ' << rho[lev]->boxArray().size() << ' ' << 0 << '\n';
#else
            HeaderFile << lev << ' ' << rho[lev].boxArray().size() << ' ' << 0 << '\n';
#endif
            HeaderFile << 0 << '\n';    // # time steps at this level
    
#ifdef UNIQUE_PTR
            for (int i = 0; i < rho[lev]->boxArray().size(); ++i)
#else
            for (int i = 0; i < rho[lev].boxArray().size(); ++i)
#endif
            {
#ifdef UNIQUE_PTR
                RealBox loc = RealBox(rho[lev]->boxArray()[i],geom[lev].CellSize(),geom[lev].ProbLo());
#else
                RealBox loc = RealBox(rho[lev].boxArray()[i],geom[lev].CellSize(),geom[lev].ProbLo());
#endif
                for (int n = 0; n < BL_SPACEDIM; n++)
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
            }
    
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            HeaderFile << PathNameInHeader << '\n';
        }

#ifdef UNIQUE_PTR
        MultiFab data(rho[lev]->boxArray(), rho[lev]->DistributionMap(), nData, 0);
#else
        MultiFab data(rho[lev].boxArray(), rho[lev]->DistributionMap(), nData, 0);
#endif
        
        // dst, src, srccomp, dstcomp, numcomp, nghost
        /*
        * srccomp: the component to copy
        * dstcmop: the component where to copy
        * numcomp: how many components to copy
        */
#ifdef UNIQUE_PTR
        MultiFab::Copy(data, *rho[lev],    0, 0, 1, 0);
        MultiFab::Copy(data, *phi[lev],    0, 1, 1, 0);
        MultiFab::Copy(data, *efield[lev], 0, 2, 3, 0); // (Ex, Ey, Ez)
#else
        MultiFab::Copy(data, rho[lev],    0, 0, 1, 0);
        MultiFab::Copy(data, phi[lev],    0, 1, 1, 0);
        MultiFab::Copy(data, efield[lev], 0, 2, 3, 0); // (Ex, Ey, Ez)
#endif
//         MultiFab::Copy(data, efield[lev], 0, 2, 1, 0);   // Ex
//         MultiFab::Copy(data, efield[lev], 1, 3, 1, 0);   // Ey
//         MultiFab::Copy(data, efield[lev], 2, 4, 1, 0);   // Ez
        
        
        //
        // Use the Full pathname when naming the MultiFab.
        //
        std::string TheFullPath = FullPath;
        TheFullPath += BaseName;
    
        VisMF::Write(data,TheFullPath);
    }
}
