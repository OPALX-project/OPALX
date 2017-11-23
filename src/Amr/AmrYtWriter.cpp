#include "AmrYtWriter.h"

#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"


AmrYtWriter::AmrYtWriter(int step) : step_m(step)
{
    namespace fs = boost::filesystem;
    
    fs::path dir = OpalData::getInstance()->getInputBasename();
    dir_m = dir.parent_path() / "data" / "amr" / "yt";
    
    if ( Ippl::myNode() == 0 && !fs::exists(dir_m) ) {
        try {
            fs::create_directories( dir_m );
        } catch(const fs::filesystem_error& ex) {
            throw OpalException("AmrPythonWriter::AmrPythonWriter()",
                                ex.what());
        }
    }
    
    Ippl::Comm->barrier();
}


void AmrYtWriter::writeFields(const amr::AmrFieldContainer_t& rho,
                              const amr::AmrFieldContainer_t& phi,
                              const amr::AmrFieldContainer_t& efield,
                              const amr::AmrIntArray_t& refRatio,
                              const amr::AmrGeomContainer_t& geom,
                              const double& time,
                              const double& scale)
{
    /* We need to scale the geometry and cell sizes according to the
     * particle mapping
     */
    std::string dir = amrex::Concatenate((dir_m / "plt").string(), step_m, 10);
    
    int nLevels = rho.size();
    //
    // Only let 64 CPUs be writing at any one time.
    //
    amrex::VisMF::SetNOutFiles(64);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if ( Ippl::myNode() == 0 )
        if (!amrex::UtilCreateDirectory(dir, 0755))
            amrex::CreateDirectoryFailed(dir);
    //
    // Force other processors to wait till directory is built.
    //
    Ippl::Comm->barrier();

    std::string HeaderFileName = dir + "/OpalAmrHeader";

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int nData = rho[0]->nComp() + phi[0]->nComp() + efield[0]->nComp();
    
    if ( Ippl::myNode() == 0 )
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
        for (int ivar = 1; ivar <= rho[0]->nComp(); ivar++)
          HeaderFile << "rho\n";
        
        for (int ivar = 1; ivar <= phi[0]->nComp(); ivar++)
          HeaderFile << "potential\n";
        
          HeaderFile << "Ex\nEy\nEz\n";
        
        // dimensionality
        HeaderFile << AMREX_SPACEDIM << '\n';
        
        // time
        HeaderFile << time << '\n';
        HeaderFile << nLevels - 1 << std::endl; // maximum level number (0=single level)
        
        // physical domain
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            HeaderFile << geom[0].ProbLo(i) / scale << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            HeaderFile << geom[0].ProbHi(i) / scale << ' ';
        HeaderFile << std::endl;
        
        // reference ratio
        for (int i = 0; i < refRatio.size(); ++i)
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
            for (int k = 0; k < AMREX_SPACEDIM; k++)
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
        if ( Ippl::myNode() == 0 )
            if (!amrex::UtilCreateDirectory(FullPath, 0755))
                amrex::CreateDirectoryFailed(FullPath);
        //
        // Force other processors to wait till directory is built.
        //
        Ippl::Comm->barrier();
        
        if ( Ippl::myNode() == 0 )
        {
            HeaderFile << lev << ' ' << rho[lev]->boxArray().size() << ' ' << 0 << '\n';
            HeaderFile << 0 << '\n';    // # time steps at this level
    
            for (int i = 0; i < rho[lev]->boxArray().size(); ++i)
            {
                amrex::Real dx[3] = {
                    geom[lev].CellSize(0),
                    geom[lev].CellSize(1),
                    geom[lev].CellSize(2)
                };
                dx[0] /= scale;
                dx[1] /= scale;
                dx[2] /= scale;
                amrex::Real lo[3] = {
                    geom[lev].ProbLo(0),
                    geom[lev].ProbLo(1),
                    geom[lev].ProbLo(2)
                };
                lo[0] /= scale;
                lo[1] /= scale;
                lo[2] /= scale;
                amr::AmrDomain_t loc = amr::AmrDomain_t(rho[lev]->boxArray()[i],
                                                        &dx[0], &lo[0]);
//                                                         geom[lev].CellSize(),
//                                                         geom[lev].ProbLo());
                for (int n = 0; n < AMREX_SPACEDIM; n++)
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
            }
    
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            HeaderFile << PathNameInHeader << '\n';
        }
        
        //
        // We combine all of the multifabs 
        //
        amr::AmrField_t data(rho[lev]->boxArray(), rho[lev]->DistributionMap(), nData, 0);
        
        //
        // Cull data -- use no ghost cells.
        //
        // dst, src, srccomp, dstcomp, numcomp, nghost
        /*
        * srccomp: the component to copy
        * dstcmop: the component where to copy
        * numcomp: how many components to copy
        */
        amr::AmrField_t::Copy(data, *rho[lev],    0, 0, 1, 0);
        amr::AmrField_t::Copy(data, *phi[lev],    0, 1, 1, 0);
        amr::AmrField_t::Copy(data, *efield[lev], 0, 2, 3, 0); // (Ex, Ey, Ez)
        
        //
        // Use the Full pathname when naming the MultiFab.
        //
        std::string TheFullPath = FullPath;
        TheFullPath += BaseName;
    
        amrex::VisMF::Write(data,TheFullPath);
    }
}


void AmrYtWriter::writeBunch(const AmrPartBunch* bunch_p,
                             const double& scale)
{
    
    throw OpalException("AmrYtWriter::writeBunch",
                        "Not yet implemented.");
}
