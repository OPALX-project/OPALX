#include "AmrSliceWriter.h"

#include "AbstractObjects/OpalData.h"

AmrSliceWriter::AmrSliceWriter()
{
    namespace fs = boost::filesystem;
    
    fs::path dir = OpalData::getInstance()->getInputBasename();
    
    filename_m = dir.filename().string();
    
    dir_m = dir.parent_path() / "data";
    
    if ( Ippl::myNode() == 0 && !fs::exists(dir_m) ) {
        try {
            fs::create_directories( dir_m );
        } catch(const fs::filesystem_error& ex) {
            throw OpalException("AmrSliceWriter::AmrSliceWriter()",
                                ex.what());
        }
    }
    
    Ippl::Comm->barrier();
}


void AmrSliceWriter::writeFields(const amr::AmrFieldContainer_t& rho,
                                 const amr::AmrFieldContainer_t& phi,
                                 const amr::AmrFieldContainer_t& efield,
                                 const amr::AmrIntArray_t& refRatio,
                                 const amr::AmrGeomContainer_t& geom,
                                 const double& time)
{
    namespace fs = boost::filesystem;
    
    std::string step = std::to_string(int(time));
    std::ofstream fstr;
    fstr.precision(9);
    
    /*
     * dump the charge density field
     */
    INFOMSG("*** START DUMPING SCALAR FIELD ***" << endl);
    std::string rho_fn = (dir_m / ( filename_m +
                                    std::string("-rho_scalar-") +
                                    step
                                  )
                         ).string();
    
    fstr.open(rho_fn.c_str(), std::ios::out);
    
    int level = 0;
    for (amrex::MFIter mfi(*rho[level]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& fab = (*rho[level])[mfi];
        
        for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
            for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                    amr::AmrIntVect_t iv(x, y, z);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    fstr << x + 1 << " " << y + 1 << " " << z + 1 << " "
                         << fab(iv, 0)  << std::endl;
                }
            }
        }
    }
    fstr.close();
    INFOMSG("*** FINISHED DUMPING SCALAR FIELD ***" << endl);
    
    /*
     * dump the electrostatic potential
     */
    INFOMSG("*** START DUMPING SCALAR FIELD ***" << endl);
    std::string phi_fn = (dir_m / ( filename_m +
                                    std::string("-phi_scalar-") +
                                    step
                                    )
                         ).string();
    
    fstr.open(phi_fn.c_str(), std::ios::out);
    
    for (amrex::MFIter mfi(*phi[level]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& fab = (*phi[level])[mfi];
        
        for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
            for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                    amr::AmrIntVect_t iv(x, y, z);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    fstr << x + 1 << " " << y + 1 << " " << z + 1 << " "
                         << fab(iv, 0)  << std::endl;
                }
            }
        }
    }
    fstr.close();
    INFOMSG("*** FINISHED DUMPING SCALAR FIELD ***" << endl);
    
    /*
     * dump the electric field
     */
    INFOMSG("*** START DUMPING E FIELD ***" << endl);
    std::string e_field = (dir_m / ( filename_m +
                                     std::string("-e_field-") +
                                     step
                                    )
                          ).string();
    
    fstr.open(e_field.c_str(), std::ios::out);
        
    for (amrex::MFIter mfi(*efield[level]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& fab = (*efield[level])[mfi];
        
        for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
            for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                    amr::AmrIntVect_t iv(x, y, z);
                        // add one in order to have same convention as PartBunch::computeSelfField()
                    fstr << x + 1 << " " << y + 1 << " " << z + 1 << " ( "
                         << fab(iv, 0) << " , " << fab(iv, 1) << " , " << fab(iv, 2) << " )" << std::endl;
                }
            }
        }
    }
    fstr.close();
    INFOMSG("*** FINISHED DUMPING E FIELD ***" << endl);
}

void AmrSliceWriter::writeBunch(const AmrPartBunch* bunch_p) {
    throw OpalException("AmrSliceWriter::writeBunch()",
                        "Writing a slice of a bunch is not supported.");
}