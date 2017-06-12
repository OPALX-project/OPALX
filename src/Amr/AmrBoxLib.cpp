#include "AmrBoxLib.h"

#include "Algorithms/AmrPartBunch.h"
// #include "Solvers/AmrPoissonSolver.h"
#include "Solvers/PoissonSolver.h"

#include "AmrBoxLib_F.h"
#include <AMReX_MultiFabUtil.H>

#include <AMReX_ParmParse.H> // used in initialize function

AmrBoxLib::AmrBoxLib() : AmrObject(),
                         amrex::AmrMesh(),
                         nChargePerCell_m(0),
                         bunch_mp(nullptr),
                         layout_mp(nullptr),
                         rho_m(0),
                         phi_m(0),
                         eg_m(0),
                         fieldDBGStep_m(0)
{}


AmrBoxLib::AmrBoxLib(TaggingCriteria tagging,
                     double scaling,
                     double nCharge)
    : AmrObject(tagging, scaling, nCharge),
      amrex::AmrMesh(),
      nChargePerCell_m(0),
      bunch_mp(nullptr),
      layout_mp(nullptr),
      rho_m(0),
      phi_m(0),
      eg_m(0),
      fieldDBGStep_m(0)
{}


// AmrBoxLib::AmrBoxLib(const DomainBoundary_t& realbox,
//               const NDIndex<3>& nGridPts,
//               short maxLevel,
//               const RefineRatios_t& refRatio)
//     : AmrObject(realbox, nGridPts, maxLevel, refRatio)
// {
//     
//     
// }



AmrBoxLib::AmrBoxLib(const AmrDomain_t& domain,
                     const AmrIntArray_t& nGridPts,
                     short maxLevel)
    : AmrObject(),
      amrex::AmrMesh(&domain, maxLevel, nGridPts, 0 /* cartesian */),
      nChargePerCell_m(maxLevel + 1),
      bunch_mp(nullptr),
      layout_mp(nullptr),
      rho_m(maxLevel + 1),
      phi_m(maxLevel + 1),
      eg_m(maxLevel + 1),
      fieldDBGStep_m(0)
{}


AmrBoxLib::AmrBoxLib(const AmrDomain_t& domain,
                     const AmrIntArray_t& nGridPts,
                     short maxLevel,
                     AmrPartBunch* bunch)
    : AmrObject(),
      amrex::AmrMesh(&domain, maxLevel, nGridPts, 0 /* cartesian */),
      nChargePerCell_m(maxLevel + 1),
      bunch_mp(bunch),
      layout_mp(static_cast<AmrLayout_t*>(&bunch->getLayout())),
      rho_m(maxLevel + 1),
      phi_m(maxLevel + 1),
      eg_m(maxLevel + 1),
      fieldDBGStep_m(0)
{
    /*
     * The layout needs to know how many levels we can make.
     */
    layout_mp->resize(maxLevel);
    
    initBaseLevel_m(nGridPts);
    
//     layout_p->define(this->Geom(),
//                      this->boxArray(),
//                      this->DistributionMap());
    
//     bunch->update();
    
//     std::cout << "After resize" << std::endl;
                  
//     AmrPartBunch::pbase_t* amrpbase_p = bunch->getAmrParticleBase();
//     amrpbase_p->resizeContainerGDB(maxLevel);
//     amrpbase_p->setGDB(this->Geom(),
//                        this->boxArray(),
//                        this->DistributionMap(),
//                        this->
//                       );
//     for (int i = 0; i < maxLevel + 1; ++i) {
//         std::cout << "level " << i << ": " << this->Geom(i) << std::endl;
//     }
    
    // set mesh spacing of bunch
    updateMesh();
    
//     std::cout << "After updateMesh" << std::endl;
}


void AmrBoxLib::setBunch(AmrPartBunch* bunch) {
    bunch_mp = bunch;
    layout_mp = static_cast<AmrLayout_t*>(&bunch->getLayout());
}


// void AmrBoxLib::initialize(const DomainBoundary_t& lower,
//                            const DomainBoundary_t& upper,
//                            const NDIndex<3>& nGridPts,
//                            short maxLevel,
//                            const RefineRatios_t& refRatio)
// {
//     // setup the physical domain
//     RealBox domain;
//     for (int i = 0; i < 3; ++i) {
//         domain.setLo(i, lower[i]); // m
//         domain.setHi(i, upper[i]); // m
//     }
//     
//     Array<int> nCells(3);
//     for (int i = 0; i < 3; ++i)
//         nCells = nGridPts[i];
//     
//     amrex::AmrMesh::Initialize();
//     Geometry::Setup(&domain, 0 /* cartesian coordinates */);
//     amrex::AmrMesh::InitAmrCore(maxLevel, nCells);
// }


void AmrBoxLib::regrid(int lbase, int lfine, double time) {
    int new_finest = 0;
    AmrGridContainer_t new_grids(finest_level+2);
    
    MakeNewGrids(lbase, time, new_finest, new_grids);

    BL_ASSERT(new_finest <= finest_level+1);
    
    for (int lev = lbase+1; lev <= new_finest; ++lev)
    {
        if (lev <= finest_level) // an old level
        {
            if (new_grids[lev] != grids[lev]) // otherwise nothing
            {
                AmrProcMap_t new_dmap(new_grids[lev]);
                RemakeLevel(lev, time, new_grids[lev], new_dmap);
	    }
	}
	else  // a new level
	{
	    AmrProcMap_t new_dmap(new_grids[lev]);
	    MakeNewLevel(lev, time, new_grids[lev], new_dmap);
        }
    }
    
    for (int lev = new_finest+1; lev <= finest_level; ++lev) {
        ClearLevel(lev);
        ClearBoxArray(lev);
        ClearDistributionMap(lev);
    }
    
    finest_level = new_finest;
    
    layout_mp->setFinestLevel(finest_level);
}


AmrBoxLib::VectorPair_t AmrBoxLib::getEExtrema() {
    Vector_t maxE = Vector_t(0, 0, 0), minE = Vector_t(0, 0, 0);
    /* check for maximum / minimum values over all levels
     * and all components, i.e. Ex, Ey, Ez
     */
    for (unsigned int lev = 0; lev < eg_m.size(); ++lev) {
        for (std::size_t i = 0; i < bunch_mp->Dimension; ++i) {
            /* calls:
             *  - max(comp, nghost (to search), local)
             *  - min(cmop, nghost, local)
             */
            double max = eg_m[lev]->max(i, 0, false);
            maxE[i] = (maxE[i] < max) ? max : maxE[i];
            
            double min = eg_m[lev]->min(i, 0, false);
            minE[i] = (minE[i] > min) ? min : minE[i];
        }
    }
    
    return VectorPair_t(maxE, minE);
}


double AmrBoxLib::getRho(int x, int y, int z) {
    //TODO
//     for (int lev = phi_m.size() - 2; lev >= 0; lev--)
//         BoxLib::average_down(phi_m[lev+1], phi_m[lev], 0, 1, this->refRatio(lev));
    throw OpalException("AmrBoxLib::getRho(x, y, z) ", "Not yet Implemented.");
    return 0.0;
}


void AmrBoxLib::computeSelfFields() {
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    
    // bring on Amr domain
    Vector_t rmin, rmax;
    bounds(bunch_mp->R, rmin, rmax);
    
    
    double scalefactor = layout_mp->domainMapping(*amrpbase_p);
    
    
    /// from charge (C) to charge density (C/m^3).
    bunch_mp->Q *= bunch_mp->dt;
    amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, finest_level);
    bunch_mp->Q /= bunch_mp->dt;
    
    int baseLevel = 0;
    int finestLevel = (&amrpbase_p->getAmrLayout())->finestLevel();
    
    int nLevel = finestLevel + 1;
    rho_m.resize(nLevel);
    phi_m.resize(nLevel);
    eg_m.resize(nLevel);
    
    //calculating mesh-scale factor
    double gammaz = sum(bunch_mp->P)[2] / bunch_mp->getTotalNum();
    gammaz *= gammaz;
    gammaz = std::sqrt(gammaz + 1.0);
    
    double tmp = 1.0 / bunch_mp->getdT() / gammaz * scalefactor;
    for (int i = 0; i <= finestLevel; ++i) {
        this->rho_m[i]->mult(tmp, 0, 1);
    }
    
#ifdef DBG_SCALARFIELD
    if ( Ippl::getNodes() > 1 )
        throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ", "Dumping only in serial execution.");

    INFOMSG("*** START DUMPING SCALAR FIELD ***" << endl);
    std::ofstream fstr1;
    fstr1.precision(9);

    std::string SfileName = OpalData::getInstance()->getInputBasename();

    std::string rho_fn = std::string("data/") + SfileName + std::string("-rho_scalar-") + std::to_string(fieldDBGStep_m);
    fstr1.open(rho_fn.c_str(), std::ios::out);
    
    int level = 0;
    for (amrex::MFIter mfi(*rho_m[level]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& fab = (*rho_m[level])[mfi];
        
        for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
            for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                    AmrIntVect_t iv(x, y, z);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    fstr1 << x + 1 << " " << y + 1 << " " << z + 1 << " "
                          << fab(iv, 0)  << std::endl;
                }
            }
        }
    }
    fstr1.close();
    INFOMSG("*** FINISHED DUMPING SCALAR FIELD ***" << endl);
#endif
    
    for (int i = 0; i <= finestLevel; ++i)
        this->rho_m[i]->mult(-1.0 / Physics::epsilon_0, 0, 1);
    
    // charge density is in rho_m
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    
    IpplTimings::startTimer(bunch_mp->compPotenTimer_m);
    solver->solve(rho_m, phi_m, eg_m, baseLevel, finestLevel);
    IpplTimings::stopTimer(bunch_mp->compPotenTimer_m);
    
    // apply scale of electric-field in order to undo the transformation
    for (int i = 0; i <= finestLevel; ++i)
        this->eg_m[i]->mult(scalefactor, 0, 3);
    
#ifdef DBG_SCALARFIELD
    INFOMSG("*** START DUMPING SCALAR FIELD ***" << endl);
    std::ofstream fstr2;
    fstr2.precision(9);

    std::string phi_fn = std::string("data/") + SfileName + std::string("-phi_scalar-") + std::to_string(fieldDBGStep_m);
    fstr2.open(phi_fn.c_str(), std::ios::out);
    
    for (amrex::MFIter mfi(*phi_m[level]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& fab = (*phi_m[level])[mfi];
        
        for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
            for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                    AmrIntVect_t iv(x, y, z);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    fstr2 << x + 1 << " " << y + 1 << " " << z + 1 << " "
                          << fab(iv, 0) << std::endl;
                }
            }
        }
    }
    fstr2.close();
    INFOMSG("*** FINISHED DUMPING SCALAR FIELD ***" << endl);
#endif
    
#ifdef DBG_SCALARFIELD
        INFOMSG("*** START DUMPING E FIELD ***" << endl);
        //ostringstream oss;
        //MPI_File file;
        //MPI_Status status;
        //MPI_Info fileinfo;
        //MPI_File_open(Ippl::getComm(), "rho_scalar", MPI_MODE_WRONLY | MPI_MODE_CREATE, fileinfo, &file);
        std::ofstream fstr;
        fstr.precision(9);

        std::string e_field = std::string("data/") + SfileName + std::string("-e_field-") + std::to_string(fieldDBGStep_m);
        fstr.open(e_field.c_str(), std::ios::out);
        
        for (amrex::MFIter mfi(*eg_m[level]); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            const amrex::FArrayBox& fab = (*eg_m[level])[mfi];
            
            for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
                for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                    for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                        AmrIntVect_t iv(x, y, z);
                        // add one in order to have same convention as PartBunch::computeSelfField()
                        fstr << x + 1 << " " << y + 1 << " " << z + 1 << " ( "
                             << fab(iv, 0) << " , " << fab(iv, 1) << " , " << fab(iv, 2) << " )" << std::endl;
                    }
                }
            }
        }

        fstr.close();
        fieldDBGStep_m++;

        INFOMSG("*** FINISHED DUMPING E FIELD ***" << endl);
#endif
    
    amrpbase_p->gather(bunch_mp->Ef, this->eg_m, bunch_mp->R, 0, finest_level);
    
    layout_mp->domainMapping(*amrpbase_p, true);
    
    bunch_mp->Ef *= Vector_t(gammaz,
                             gammaz,
                             1.0 / gammaz);
    
    /** Magnetic field in x and y direction induced by the eletric field
     *
     *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
     *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
     *  \f[ B_z = B_z^{'} = 0 \f]
     *
     */
    double betaC = sqrt(gammaz * gammaz - 1.0) / gammaz / Physics::c;

    bunch_mp->Bf(0) = bunch_mp->Bf(0) - betaC * bunch_mp->Ef(1);
    bunch_mp->Bf(1) = bunch_mp->Bf(1) + betaC * bunch_mp->Ef(0);
//     throw OpalException("AmrBoxLib::computeSelfFields() ", "Not yet Implemented.");
}


void AmrBoxLib::computeSelfFields(int bin) {
    //TODO
    throw OpalException("AmrBoxLib::computeSelfFields(int bin) ", "Not yet Implemented.");
}


void AmrBoxLib::computeSelfFields_cycl(double gamma) {   
    /*
     * The potential is not scaled according to domain modification.
     * 
     */
    
    
    /*
     * scatter charges onto grid
     */
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    
    // map on Amr domain
//     bunch_mp->python_format(0);
    double scalefactor = layout_mp->domainMapping(*amrpbase_p);
//     bunch_mp->python_format(1); std::cout << "Written." << std::endl; std::cin.get();
    
    /// from charge (C) to charge density (C/m^3).
    amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, finest_level);
    
    int baseLevel = 0;
    int nLevel = finest_level + 1;
    double invGamma = 1.0 / gamma;
    
    /// Lorentz transformation
    /// In particle rest frame, the longitudinal length (y for cyclotron) enlarged
    for (int i = 0; i <= finest_level; ++i) {
        this->rho_m[i]->mult(invGamma * scalefactor, 0 /*comp*/, 1 /*ncomp*/);
        
        if ( this->rho_m[i]->contains_nan(false) )
            throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ",
                                "NANs at level " + std::to_string(i) + ".");
    }
    
    
#ifdef DBG_SCALARFIELD
    if ( Ippl::getNodes() > 1 )
        throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ", "Dumping only in serial execution.");

    INFOMSG("*** START DUMPING SCALAR FIELD ***" << endl);
    std::ofstream fstr1;
    fstr1.precision(9);

    std::string SfileName = OpalData::getInstance()->getInputBasename();

    std::string rho_fn = std::string("data/") + SfileName + std::string("-rho_scalar-") + std::to_string(fieldDBGStep_m);
    fstr1.open(rho_fn.c_str(), std::ios::out);
    
    int level = 0;
    for (amrex::MFIter mfi(*rho_m[level]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& fab = (*rho_m[level])[mfi];
        
        for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
            for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                    AmrIntVect_t iv(x, y, z);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    fstr1 << x + 1 << " " << y + 1 << " " << z + 1 << " "
                          << fab(iv, 0)  << std::endl;
                }
            }
        }
    }
    fstr1.close();
    INFOMSG("*** FINISHED DUMPING SCALAR FIELD ***" << endl);
#endif
    
    
    // charge density is in rho_m
    // calculate Possion equation (with coefficient: -1/(eps))
    for (int i = 0; i <= finest_level; ++i) {
        this->rho_m[i]->mult(-1.0 / Physics::epsilon_0, 0, 1);
    }
    
    
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    IpplTimings::startTimer(bunch_mp->compPotenTimer_m);
    solver->solve(rho_m, phi_m, eg_m, baseLevel, finest_level);
    IpplTimings::stopTimer(bunch_mp->compPotenTimer_m);
    
    // apply scale of electric-field in order to undo the transformation
    for (int i = 0; i <= finestLevel; ++i)
        this->eg_m[i]->mult(scalefactor, 0, 3);
    
    for (int i = 0; i <= finest_level; ++i) {
        if ( this->eg_m[i]->contains_nan(false) )
            throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ",
                                "Ef: NANs at level " + std::to_string(i) + ".");
    }
    
//     for (int i = 0; i <= finest_level; ++i) {
//         if ( this->phi_m[i]->contains_nan(false) )
//             throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ",
//                                 "Pot: NANs at level " + std::to_string(i) + ".");
//     }
    
#ifdef DBG_SCALARFIELD
    INFOMSG("*** START DUMPING SCALAR FIELD ***" << endl);
    std::ofstream fstr2;
    fstr2.precision(9);

    std::string phi_fn = std::string("data/") + SfileName + std::string("-phi_scalar-") + std::to_string(fieldDBGStep_m);
    fstr2.open(phi_fn.c_str(), std::ios::out);
    
    for (amrex::MFIter mfi(*phi_m[level]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& fab = (*phi_m[level])[mfi];
        
        for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
            for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                    AmrIntVect_t iv(x, y, z);
                    // add one in order to have same convention as PartBunch::computeSelfField()
                    fstr2 << x + 1 << " " << y + 1 << " " << z + 1 << " "
                          << fab(iv, 0)  << std::endl;
                }
            }
        }
    }
    fstr2.close();
    INFOMSG("*** FINISHED DUMPING SCALAR FIELD ***" << endl);
#endif
    
    /// Back Lorentz transformation
//     for (int i = 0; i <= finest_level; ++i) {
//         this->eg_m[i]->mult(gamma, 0, 3, 1); // x-direction
//         this->eg_m[i]->mult(1.0 / gamma, 1, 3, 1); // y-direction
//         this->eg_m[i]->mult(gamma, 2, 3, 1); // z-direction
//     }
    
#ifdef DBG_SCALARFIELD
        INFOMSG("*** START DUMPING E FIELD ***" << endl);
        std::ofstream fstr;
        fstr.precision(9);

        std::string e_field = std::string("data/") + SfileName + std::string("-e_field-") + std::to_string(fieldDBGStep_m);
        fstr.open(e_field.c_str(), std::ios::out);
        
        for (amrex::MFIter mfi(*eg_m[level]); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            const amrex::FArrayBox& fab = (*eg_m[level])[mfi];
            
            for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
                for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                    for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                        AmrIntVect_t iv(x, y, z);
                        // add one in order to have same convention as PartBunch::computeSelfField()
                        fstr << x + 1 << " " << y + 1 << " " << z + 1 << " ( "
                             << fab(iv, 0) << " , " << fab(iv, 1) << " , " << fab(iv, 2) << " )" << std::endl;
                    }
                }
            }
        }

        fstr.close();
        fieldDBGStep_m++;

        INFOMSG("*** FINISHED DUMPING E FIELD ***" << endl);
#endif
    
    amrpbase_p->gather(bunch_mp->Ef, this->eg_m, bunch_mp->R, 0, finest_level);
    
    // undo domain change
    layout_mp->domainMapping(*amrpbase_p, true);
    
    bunch_mp->Ef *= Vector_t(gamma,
                             invGamma,
                             gamma);
    
//     std::cout << bunch_mp->Ef[0] << std::endl; std::cin.get();
    
    /// calculate coefficient
    // Relativistic E&M says gamma*v/c^2 = gamma*beta/c = sqrt(gamma*gamma-1)/c
    // but because we already transformed E_trans into the moving frame we have to
    // add 1/gamma so we are using the E_trans from the rest frame -DW
    double betaC = std::sqrt(gamma * gamma - 1.0) * invGamma / Physics::c;
    
    /// calculate B field from E field
    bunch_mp->Bf(0) =  betaC * bunch_mp->Ef(2);
    bunch_mp->Bf(2) = -betaC * bunch_mp->Ef(0);
}


void AmrBoxLib::computeSelfFields_cycl(int bin) {
    //TODO
    throw OpalException("AmrBoxLib::computeSelfFields_cycl(int bin) ", "Not yet Implemented.");
}


void AmrBoxLib::updateMesh() {
    //FIXME What about resizing mesh, i.e. geometry?
    const AmrGeomContainer_t& geom = this->Geom();
    
    const AmrReal_t* tmp = geom[0].CellSize();
    
    Vector_t hr;
    for (int i = 0; i < 3; ++i)
        hr[i] = tmp[i];
    
    bunch_mp->setBaseLevelMeshSpacing(hr);
}

Vektor<int, 3> AmrBoxLib::getBaseLevelGridPoints() {
    const AmrGeomContainer_t& geom = this->Geom();
    const amrex::Box& bx = geom[0].Domain();
    
    const AmrIntVect_t& low = bx.smallEnd();
    const AmrIntVect_t& high = bx.bigEnd();
    
    return Vektor<int, 3>(high[0] - low[0] + 1,
                          high[1] - low[1] + 1,
                          high[2] - low[2] + 1);
}


inline int AmrBoxLib::maxLevel() {
    return amrex::AmrMesh::maxLevel();
}


inline int AmrBoxLib::finestLevel() {
    return amrex::AmrMesh::finestLevel();
}


// void AmrBoxLib::updateBunch()  {
//     const Array<Geometry>& geom = this->Geom();
//     const Array<AmrProcMap_t>& dmap = this->DistributionMap();
//     const Array<BoxArray>& ba = this->boxArray();
//     const Array<IntVect>& ref_rato = this->refRatio ();
//         
//         
//     Array<int> rr( ref_rato.size() );
//     for (unsigned int i = 0; i < rr.size(); ++i) {
//         rr[i] = ref_rato[i][0];
//     }        
//     
//     layout_mp->define(geom, ba, dmap, rr);
//         
//     bunch_mp->update();
//         
// //     this->setBunch(bunch_mp);
// }


void AmrBoxLib::RemakeLevel (int lev, AmrReal_t time,
                             const AmrGrid_t& new_grids,
                             const AmrProcMap_t& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m[lev].reset(new AmrField_t(new_grids, new_dmap, 1, 1));
    nChargePerCell_m[lev]->setVal(0.0, 1);
    
    rho_m[lev].reset(new AmrField_t(new_grids, new_dmap, 1, 1));
    rho_m[lev]->setVal(0.0, 1);
    
    /*
     * particles need to know the BoxArray
     * and DistributionMapping
     */
     layout_mp->SetParticleBoxArray(lev, new_grids);
     layout_mp->SetParticleDistributionMap(lev, new_dmap);
}


void AmrBoxLib::MakeNewLevel (int lev, AmrReal_t time,
                              const AmrGrid_t& new_grids,
                              const AmrProcMap_t& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    
    nChargePerCell_m[lev].reset(new AmrField_t(new_grids, new_dmap, 1, 1));
    nChargePerCell_m[lev]->setVal(0.0, 1);
    
    
    rho_m[lev].reset(new AmrField_t(new_grids, new_dmap, 1, 1));
    rho_m[lev]->setVal(0.0, 1);
    
    
    /*
     * particles need to know the BoxArray
     * and DistributionMapping
     */
    layout_mp->SetParticleBoxArray(lev, new_grids);
    layout_mp->SetParticleDistributionMap(lev, new_dmap);
}


void AmrBoxLib::ClearLevel(int lev) {
    nChargePerCell_m[lev].reset(nullptr);
    rho_m[lev].reset(nullptr);
    ClearBoxArray(lev);
    ClearDistributionMap(lev);
}


void AmrBoxLib::ErrorEst(int lev, TagBoxArray_t& tags,
                         AmrReal_t time, int ngrow)
{
    switch ( tagging_m ) {
        case CHARGE_DENSITY:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
        case POTENTIAL:
            tagForPotentialStrength_m(lev, tags, time, ngrow);
            break;
        case EFIELD:
            tagForEfield_m(lev, tags, time, ngrow);
            break;
        default:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
    }
}


void AmrBoxLib::MakeNewLevelFromScratch(int lev, AmrReal_t time,
                                  const AmrGrid_t& ba,
                                  const AmrProcMap_t& dm)
{
    throw OpalException("AmrBoxLib::MakeNewLevelFromScratch()", "Shouldn't be called.");
}


void AmrBoxLib::MakeNewLevelFromCoarse (int lev, AmrReal_t time,
                                        const AmrGrid_t& ba,
                                        const AmrProcMap_t& dm)
{
    throw OpalException("AmrBoxLib::MakeNewLevelFromCoarse()", "Shouldn't be called.");
}


void AmrBoxLib::tagForChargeDensity_m(int lev, TagBoxArray_t& tags,
                                      AmrReal_t time, int ngrow)
{
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    for (int i = lev; i <= finest_level; ++i)
        nChargePerCell_m[i]->setVal(0.0, 1);
    
    // bring on Amr domain
    Vector_t rmin, rmax;
    bounds(bunch_mp->R, rmin, rmax);
    
    layout_mp->domainMapping(*amrpbase_p);
    
    // the new scatter function averages the value also down to the coarsest level
    amrpbase_p->scatter(bunch_mp->Q, nChargePerCell_m,
                        bunch_mp->R, lev, finest_level);
    
    // undo domain change
    layout_mp->domainMapping(*amrpbase_p, true);
    
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;

    const AmrReal_t* dx      = geom[lev].CellSize();
    const AmrReal_t* prob_lo = geom[lev].ProbLo();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        AmrIntArray_t  itags;
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
                        ZFILL(dx), ZFILL(prob_lo), &time, &nCharge_m);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
}


void AmrBoxLib::tagForPotentialStrength_m(int lev, TagBoxArray_t& tags, AmrReal_t time, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the potential is higher than 75 percent of the maximum potential
     * value of this level.
     */
    
    // 1. Assign charge onto grid of level lev
    int baseLevel   = 0;
    int finestLevel = 0;
    
    AmrFieldContainer_t nPartPerCell(1);
    nPartPerCell[baseLevel] = std::unique_ptr<AmrField_t>(new AmrField_t(this->boxArray(lev),
                                                                         this->DistributionMap(lev),
                                                                         1, 1));
    nPartPerCell[baseLevel]->setVal(0.0, 1);
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    // single level scatter
    amrpbase_p->scatter(bunch_mp->Q, *nPartPerCell[baseLevel], bunch_mp->R, lev);
    
    if ( geom[0].isAllPeriodic() ) {
        AmrReal_t offset = 0.0;
        for (std::size_t i = 0; i < bunch_mp->getLocalNum(); ++i)
            offset += bunch_mp->Q[i];
        offset /= geom[0].ProbSize();
        for (int lev = baseLevel; lev <= finestLevel; lev++)
            nPartPerCell[lev]->plus(+offset, 0, 1, 0);
    }
    
    // 2. Solve Poisson's equation on level lev
    AmrFieldContainer_t phi(1);
    AmrFieldContainer_t efield(1);
    
//     AmrPoissonSolver<AmrBoxLib> *solver = bunch_mp->template getFieldSolver<AmrBoxLib>();
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    
    solver->solve(nPartPerCell, phi, efield, baseLevel, finestLevel);
    
    // 3. Tag cells for refinement
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;

    const AmrReal_t* dx      = geom[lev].CellSize();
    const AmrReal_t* prob_lo = geom[lev].ProbLo();
    
    AmrReal_t minp = 0.0;
    AmrReal_t maxp = 0.0;
    maxp = scaling_m * phi[baseLevel]->max(0);
    minp = scaling_m * phi[baseLevel]->min(0);
    AmrReal_t threshold = std::max( std::abs(minp), std::abs(maxp) );
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        AmrIntArray_t  itags;
        for (amrex::MFIter mfi(*phi[baseLevel],false); mfi.isValid(); ++mfi) {
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
                        BL_TO_FORTRAN_3D((*phi[baseLevel])[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &threshold);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
    
    phi[0].reset(nullptr);
    efield[0].reset(nullptr);
    nPartPerCell[0].reset(nullptr);
}


void AmrBoxLib::tagForEfield_m(int lev, TagBoxArray_t& tags, AmrReal_t time, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the efield is higher than 75 percent of the maximum efield
     * value of this level.
     */
    
    if ( bunch_mp->Q[0] == 0 )
        return;
    
    // 1. Assign charge onto grid of level lev
    int baseLevel   = 0;
    int finestLevel = 0;
    
    AmrFieldContainer_t nPartPerCell(1);
    nPartPerCell[baseLevel] = std::unique_ptr<AmrField_t>(new AmrField_t(this->boxArray(lev),
                                                                         this->DistributionMap(lev),
                                                                         1, 1));
    nPartPerCell[baseLevel]->setVal(0.0, 1);
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    // single level scatter
    amrpbase_p->scatter(bunch_mp->Q, *nPartPerCell[baseLevel], bunch_mp->R, lev);
    
    std::cout << "lev = " << lev << " " << bunch_mp->Q[0] << std::endl;
    
    if ( geom[0].isAllPeriodic() ) {
        std::cout << "All periodic" << std::endl;
        AmrReal_t offset = 0.0;
        for (std::size_t i = 0; i < bunch_mp->getLocalNum(); ++i)
            offset += bunch_mp->Q[i];
        offset /= geom[0].ProbSize();
        for (int ilev = baseLevel; ilev <= finestLevel; ilev++)
            nPartPerCell[ilev]->plus(+offset, 0, 1, 0);
    }
    
    // charge density is in rho_m
    // calculate Possion equation (with coefficient: -1/(eps))
    for (int i = baseLevel; i <= finestLevel; ++i) {
        nPartPerCell[i]->mult(-1.0 / Physics::epsilon_0, 0, 1);
    }
    
    std::cout << "sum: " << nPartPerCell[baseLevel]->sum()
              << " min: " << nPartPerCell[baseLevel]->min(0)
              << " max: " << nPartPerCell[baseLevel]->max(0) << std::endl;
    
    // 2. Solve Poisson's equation on level lev
    AmrFieldContainer_t phi;
    AmrFieldContainer_t efield;
    
//     AmrPoissonSolver<AmrBoxLib> *solver = bunch_mp->template getFieldSolver<AmrBoxLib>();
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    
    solver->solve(nPartPerCell, phi, efield, baseLevel, finestLevel);
    
    std::cout << "potential: " << phi[0]->sum() << " " << phi[0]->min(0) << " " << phi[0]->max(0) << std::endl;
    
    // 3. Tag cells for refinement
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;
    
    AmrReal_t min[3] = {0.0, 0.0, 0.0};
    AmrReal_t max[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        max[i] = scaling_m * efield[baseLevel]->max(i);
        min[i] = scaling_m * efield[baseLevel]->min(i);
        std::cout << max[i] << " " << min[i] << " " << efield[baseLevel]->max(i) << " " << efield[baseLevel]->min(i) << std::endl;
    }
    AmrReal_t threshold[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)
        threshold[i] = std::max( std::abs(min[i]), std::abs(max[i]) );
     std::cin.get();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        AmrIntArray_t  itags;
        // mfi(efield[baseLevel], true)
        for (amrex::MFIter mfi(*efield[baseLevel],false); mfi.isValid(); ++mfi) {
            const amrex::Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            amrex::TagBox&     tagfab  = tags[mfi];
            amrex::FArrayBox&  fab     = (*efield[baseLevel])[mfi];
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
//             tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        AmrIntVect_t iv(i,j,k);
                        if (std::abs(fab(iv, 0)) >= threshold[0])
                            tagfab(iv) = tagval;
                        else if (std::abs(fab(iv, 1)) >= threshold[1])
                            tagfab(iv) = tagval;
                        else if (std::abs(fab(iv, 2)) >= threshold[2])
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
    
    phi[baseLevel].reset(nullptr);
    efield[baseLevel].reset(nullptr);
    nPartPerCell[baseLevel].reset(nullptr);
    
}


void AmrBoxLib::initBaseLevel_m(const AmrIntArray_t& nGridPts) {
    // we need to set the AmrCore variable
    finest_level = 0;
    
    AmrIntVect_t low(0, 0, 0); 
    AmrIntVect_t high(nGridPts[0] - 1,
                      nGridPts[1] - 1,
                      nGridPts[2] - 1);
    
    const amrex::Box bx(low, high);
    AmrGrid_t ba(bx);
    ba.maxSize( this->maxGridSize(0) );
//     this->SetBoxArray(0, ba);
    
    AmrProcMap_t dmap;
    dmap.define(ba, amrex::ParallelDescriptor::NProcs());
//     this->SetDistributionMap(0, dmap);
    
    this->RemakeLevel(0, 0.0, ba, dmap);
    
    for (int i = 0; i < maxLevel() + 1; ++i) {
        std::cout << "level " << i << ": " << this->Geom(i) << std::endl;
    }
    
    layout_mp->define(this->Geom());
    layout_mp->define(this->ref_ratio);
    
//     const BoxArray& ba = layout_mp->ParticleBoxArray(0);
//     const AmrProcMap_t& dm = layout_mp->ParticleDistributionMap(0 /*level*/);
    
//     SetBoxArray(0 /*level*/, ba);
//     SetDistributionMap(0 /*level*/, dm);
    
//     rho_m[0] = std::unique_ptr<AmrField_t>(new AmrField_t(ba, 1, 1, dm));
//     rho_m[0]->setVal(0.0);
    
//     nChargePerCell_m[0] = std::unique_ptr<AmrField_t>(new AmrField_t(ba, 1, 1, dm));
//     nChargePerCell_m[0]->setVal(0.0);
}


// void AmrBoxLib::resizeBaseLevel_m(const AmrDomain_t& domain,
//                                   const AmrIntArray_t& nGridPts) {
//     
//     Geometry::ProbDomain(domain);
//     
//     // This says we are using Cartesian coordinates
//     int coord = 0;
//     
//         // This sets the boundary conditions to be doubly or triply periodic
//         int is_per[BL_SPACEDIM];
//         for (int i = 0; i < BL_SPACEDIM; i++) 
//             is_per[i] = 0; 
//     
//     bunch_mp->
// }

