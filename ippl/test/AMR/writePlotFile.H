#include <fstream>
#include <iomanip>

#include <Geometry.H>
#include <VisMF.H>

void
writePlotFile (const std::string& dir,
               const MultiFab&    rho,
               const MultiFab&    phi,
               const MultiFab&    efield,
               const Geometry&    geom,
	       const Real&        time)
{
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

    int nData = rho.nComp() + phi.nComp() + efield.nComp();
    
    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);
        HeaderFile << "HyperCLaw-V1.1\n";

        HeaderFile << nData << '\n';

	// variable names
        for (int ivar = 1; ivar <= rho.nComp(); ivar++) {
          HeaderFile << "rho\n";
        }
        
        for (int ivar = 1; ivar <= phi.nComp(); ivar++) {
          HeaderFile << "phi\n";
        }
        
//         for (int ivar = 1; ivar <= efield.nComp(); ivar++) {
          HeaderFile << "Ex\nEy\nEz\n";
//         }
        
        
	// dimensionality
        HeaderFile << BL_SPACEDIM << '\n';
	// time
        HeaderFile << time << '\n';
	// maximum level number (0=single level)
        HeaderFile << 0 << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom.ProbLo(i) << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom.ProbHi(i) << ' ';
        HeaderFile << '\n';
        HeaderFile << '\n';
        HeaderFile << geom.Domain() << ' ';
        HeaderFile << '\n';
        HeaderFile << 0 << ' ';
        HeaderFile << '\n';
        for (int k = 0; k < BL_SPACEDIM; k++)
            HeaderFile << geom.CellSize()[k] << ' ';
        HeaderFile << '\n';
        HeaderFile << geom.Coord() << '\n';
        HeaderFile << "0\n";
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";

    std::string Level = BoxLib::Concatenate("Level_", 0, 1);
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

    
    std::cerr << "Hallo 1" << std::endl;
    
    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile << 0 << ' ' << rho.boxArray().size() << ' ' << 0 << '\n';
        HeaderFile << 0 << '\n';

        for (int i = 0; i < rho.boxArray().size(); ++i)
        {
            RealBox loc = RealBox(rho.boxArray()[i],geom.CellSize(),geom.ProbLo());
            for (int n = 0; n < BL_SPACEDIM; n++)
                HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
        }

        std::string PathNameInHeader = Level;
        PathNameInHeader += BaseName;
        HeaderFile << PathNameInHeader << '\n';
    }
    
    std::cerr << "Hallo 2 nData = " << nData << std::endl;
    
    MultiFab data(rho.boxArray(), nData, 0);
    
    std::cerr << data.nComp() << " " << nData << std::endl;
    
    std::cerr << "Hallo 3" << std::endl;
    // dst, src, srccomp, dstcomp, numcomp, nghost
    /*
     * srccomp: the component to copy
     * dstcmop: the component where to copy
     * numcomp: how many compoents to copy
     */
    MultiFab::Copy(data, rho,    0, 0, 1, 0);
    MultiFab::Copy(data, phi,    0, 1, 1, 0);
    MultiFab::Copy(data, efield, 0, 2, 1, 0);   // Ex
    MultiFab::Copy(data, efield, 1, 3, 1, 0);   // Ey
    MultiFab::Copy(data, efield, 2, 4, 1, 0);   // Ez
    
    
    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;

    VisMF::Write(data,TheFullPath);
}
