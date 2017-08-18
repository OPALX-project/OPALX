#include "AmrPythonWriter.h"

#include "Utilities/OpalException.h"
#include "AbstractObjects/OpalData.h"

#include <iomanip>
#include <sstream>


AmrPythonWriter::AmrPythonWriter()
{
    namespace fs = boost::filesystem;
    
    fs::path dir = OpalData::getInstance()->getInputBasename();
    dir_m = dir.parent_path() / "data" / "amr" / "python";
    
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


void AmrPythonWriter::writeFields(const amr::AmrFieldContainer_t& rho,
                                  const amr::AmrFieldContainer_t& phi,
                                  const amr::AmrFieldContainer_t& efield,
                                  const amr::AmrIntArray_t& refRatio,
                                  const amr::AmrGeomContainer_t& geom,
                                  const double& time,
                                  const double& scale)
{
    throw OpalException("AmrPythonWriter::writeFields",
                        "Not yet implemented.");
}

void AmrPythonWriter::writeBunch(const AmrPartBunch* bunch_p,
                                 const double& scale)
{
    /* We need to scale the geometry and cell sizes according to the
     * particle mapping
     */
    namespace fs = boost::filesystem;
    
    std::stringstream sstep;
    
    long long step = bunch_p->getLocalTrackStep();
    
    sstep << std::setfill('0') << std::setw(10) << std::to_string(step);
    
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
                amrex::Real dx[3] = {
                    geom.CellSize(0),
                    geom.CellSize(1),
                    geom.CellSize(2)
                };
                dx[0] /= scale;
                dx[1] /= scale;
                dx[2] /= scale;
                amrex::Real lo[3] = {
                    geom.ProbLo(0),
                    geom.ProbLo(1),
                    geom.ProbLo(2)
                };
                lo[0] /= scale;
                lo[1] /= scale;
                lo[2] /= scale;
                for (int grid = 0; grid < playout->ParticleBoxArray(lev).size(); ++grid) {
                    msg_grid << lev << ' ';
                    amr::AmrDomain_t loc = amr::AmrDomain_t(playout->boxArray(lev)[grid],
                                                            &dx[0], &lo[0]);
//                                                             geom.CellSize(),
//                                                             geom.ProbLo());
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