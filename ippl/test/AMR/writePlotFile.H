#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>

using amrex::MultiFab;
using amrex::Array;
using amrex::Geometry;
using amrex::VisMF;
using amrex::Real;
using amrex::RealBox;

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
writePlotFile (const std::string& dir,
               const std::vector<std::unique_ptr<MultiFab> >& rho,
               const std::vector<std::unique_ptr<MultiFab> >& phi,
               const std::vector<std::unique_ptr<MultiFab> >& efield,
               const Array<int>& refRatio,
               const Array<Geometry>& geom,
	       const Real&        time,
               double scale = 1.0)
{
    int nLevels = rho.size();
    //
    // Only let 64 CPUs be writing at any one time.
    //
    VisMF::SetNOutFiles(64);
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

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int nData = rho[0]->nComp() + phi[0]->nComp() + efield[0]->nComp();
    
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);
//         HeaderFile << "OpalAmr-V1.0\n";
        HeaderFile << "HyperClaw-V1.1\n";

        HeaderFile << nData << '\n';

	// variable names
        for (int ivar = 1; ivar <= rho[0]->nComp(); ivar++)
          HeaderFile << "rho\n";
        
        for (int ivar = 1; ivar <= phi[0]->nComp(); ivar++)
          HeaderFile << "phi\n";
        
//         for (int ivar = 1; ivar <= efield.nComp(); ivar++) {
#if BL_SPACEDIM == 3
          HeaderFile << "Ex\nEy\nEz\n";
#else
          HeaderFile << "Ex\nEy\n";
#endif
//         }
        
        
        // dimensionality
        HeaderFile << BL_SPACEDIM << '\n';
        
        // time
        HeaderFile << time << '\n';
        HeaderFile << nLevels - 1 << std::endl; // maximum level number (0=single level)
        
        // physical domain
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom[0].ProbLo(i) / scale << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom[0].ProbHi(i) / scale << ' ';
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
                HeaderFile << geom[i].CellSize()[k] / scale << ' ';
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
            HeaderFile << lev << ' ' << rho[lev]->boxArray().size() << ' ' << 0 << '\n';
            HeaderFile << 0 << '\n';    // # time steps at this level
    
            for (int i = 0; i < rho[lev]->boxArray().size(); ++i)
            {
                amrex::Real dx[BL_SPACEDIM] = {
                    D_DECL(
                        geom[lev].CellSize(0),
                        geom[lev].CellSize(1),
                        geom[lev].CellSize(2)
                    )
                };
                dx[0] /= scale;
                dx[1] /= scale;
#if BL_SPACEDIM == 3
                dx[2] /= scale;
#endif
                amrex::Real lo[BL_SPACEDIM] = {
                    D_DECL(
                        geom[lev].ProbLo(0),
                        geom[lev].ProbLo(1),
                        geom[lev].ProbLo(2)
                    )
                };
                lo[0] /= scale;
                lo[1] /= scale;
#if BL_SPACEDIM == 3
                lo[2] /= scale;
#endif
                
                
                RealBox loc = RealBox(rho[lev]->boxArray()[i],
                                      &dx[0], &lo[0]);
                
                for (int n = 0; n < BL_SPACEDIM; n++)
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
            }
    
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            HeaderFile << PathNameInHeader << '\n';
        }

        MultiFab data(rho[lev]->boxArray(), rho[lev]->DistributionMap(), nData, 0);
        
        // dst, src, srccomp, dstcomp, numcomp, nghost
        /*
        * srccomp: the component to copy
        * dstcmop: the component where to copy
        * numcomp: how many components to copy
        */
        MultiFab::Copy(data, *rho[lev],    0, 0, 1, 0);
        MultiFab::Copy(data, *phi[lev],    0, 1, 1, 0);
        MultiFab::Copy(data, *efield[lev], 0, 2, BL_SPACEDIM, 0); // (Ex, Ey, Ez)
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
