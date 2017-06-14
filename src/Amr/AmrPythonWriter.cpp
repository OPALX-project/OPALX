#include "AmrPythonWriter.h"

#include "Utilities/OpalException.h"
#include "AbstractObjects/OpalData.h"

#include <iomanip>
#include <sstream>


AmrPythonWriter::AmrPythonWriter()
{
    namespace fs = boost::filesystem;
    
    fs::path dir = OpalData::getInstance()->getInputBasename();
    dir_m = dir.parent_path() / "data" / "amr";
    
    if ( Ippl::myNode() == 0 && !fs::exists(dir_m) )
        fs::create_directory( dir_m );
    
    Ippl::Comm->barrier();
}


void AmrPythonWriter::writeGrids(const std::string& dir,
                                 const amr::AmrFieldContainer_t& rho,
                                 const amr::AmrFieldContainer_t& phi,
                                 const amr::AmrFieldContainer_t& efield,
                                 const amr::AmrIntArray_t& refRatio,
                                 const amr::AmrGeomContainer_t& geom,
                                 const double& time)
{
    throw OpalException("AmrPythonWriter::writeGrids",
                        "Not yet implemented.");
}

void AmrPythonWriter::writeBunch(const AmrPartBunch* bunch_p) {
    namespace fs = boost::filesystem;
    
    std::stringstream sstep;
    
    long long step = bunch_p->getLocalTrackStep();
    
    sstep << std::setfill('0') << std::setw(6) << std::to_string(step);
    
    fs::path gDir = dir_m / ("grids_" + sstep.str() + ".dat");
    fs::path pDir = dir_m / ("bunch_" + sstep.str() + ".dat");
    
    std::string fGrid     = gDir.string();
    std::string fParticle = pDir.string();
    
    AmrLayout_t* playout = (AmrLayout_t*)(&bunch_p->getLayout());
    
    // only one core writes at a time
    Inform::WriteMode wm = Inform::OVERWRITE;
    for (int i = 0; i < Ippl::getNodes(); ++i) {
        if ( i == Ippl::myNode() ) {
            wm = (i == 0) ? wm : Inform::APPEND;
            
            Inform msg_grid("", fGrid.c_str(), wm, i);
            for (int lev = 0; lev < playout->finestLevel() + 1; ++lev) {
                const amr::AmrGeometry_t& geom = playout->Geom(lev);
                for (int grid = 0; grid < playout->ParticleBoxArray(lev).size(); ++grid) {
                    msg_grid << lev << ' ';
                    amr::AmrDomain_t loc = amr::AmrDomain_t(playout->boxArray(lev)[grid],
                                                            geom.CellSize(),
                                                            geom.ProbLo());
                    for (int n = 0; n < BL_SPACEDIM; n++)
                        msg_grid << loc.lo(n) << ' ' << loc.hi(n) << ' ';
                    msg_grid << endl;
                }
            }
            
            Inform msg_particles("", fParticle.c_str(), wm, i);
            for (size_t i = 0; i < bunch_p->getLocalNum(); ++i) {
                msg_particles << bunch_p->R[i](0) << " "
                              << bunch_p->R[i](1) << " "
                              << bunch_p->R[i](2) << " "
                              << bunch_p->P[i](0) << " "
                              << bunch_p->P[i](1) << " "
                              << bunch_p->P[i](2)
                              << endl;
            }
        }
        Ippl::Comm->barrier();
    }
}