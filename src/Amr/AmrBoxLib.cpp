#include "AmrBoxLib.h"

#include "Algorithms/AmrPartBunch.h"
// #include "Solvers/AmrPoissonSolver.h"
#include "Solvers/PoissonSolver.h"

#include "AmrBoxLib_F.h"
#include <MultiFabUtil.H>

#include <ParmParse.H> // used in initialize function

AmrBoxLib::AmrBoxLib() : AmrObject(),
                         AmrCore(),
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
      AmrCore(),
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
      AmrCore(&domain, maxLevel, nGridPts, 0 /* cartesian */),
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
      AmrCore(&domain, maxLevel, nGridPts, 0 /* cartesian */),
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
//     AmrCore::Initialize();
//     Geometry::Setup(&domain, 0 /* cartesian coordinates */);
//     AmrCore::InitAmrCore(maxLevel, nCells);
// }


void AmrBoxLib::regrid(int lbase, int lfine, double time) {
    int new_finest = 0;
    AmrGridContainer_t new_grids(finest_level+2);
    
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
	    }
	}
	else  // a new level
	{
	    DistributionMapping new_dmap(new_grids[lev], ParallelDescriptor::NProcs());
	    MakeNewLevel(lev, time, new_grids[lev], new_dmap);
        }
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
    
//     //FIXME now we regrid in every selffield
//     this->regrid(0, maxLevel(), 0.0);
//     
//     //FIXME Lorentz transformation
//     //scatter charges onto grid
// //     bunch_mp->Q *= bunch_mp->dt;
//     AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
//     amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, -1);
// //     bunch_mp->Q /= bunch_mp->dt;
//     
//     int baseLevel = 0;
//     int finestLevel = (&amrpbase_p->getAmrLayout())->finestLevel();
//     
//     int nLevel = finestLevel + 1;
//     rho_m.resize(nLevel);
//     phi_m.resize(nLevel);
//     eg_m.resize(nLevel);
//     
//     double invDt = 1.0 / bunch_mp->getdT() * bunch_mp->getCouplingConstant();
//     for (int i = 0; i <= finestLevel; ++i) {
//         this->rho_m[i]->mult(invDt, 0, 1);
//     }
//     
//     //calculating mesh-scale factor
//     double gammaz = sum(bunch_mp->P)[2] / bunch_mp->getTotalNum();
//     gammaz *= gammaz;
//     gammaz = std::sqrt(gammaz + 1.0);
//     
//     // charge density is in rho_m
//     PoissonSolver *solver = bunch_mp->getFieldSolver();
//     
//     IpplTimings::startTimer(bunch_mp->compPotenTimer_m);
//     solver->solve(rho_m, phi_m, eg_m, baseLevel, finestLevel);
//     IpplTimings::stopTimer(bunch_mp->compPotenTimer_m);
//     
//     amrpbase_p->gather(bunch_mp->Ef, this->eg_m, bunch_mp->R, 0, -1);
//     
//     /** Magnetic field in x and y direction induced by the eletric field
//      *
//      *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
//      *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
//      *  \f[ B_z = B_z^{'} = 0 \f]
//      *
//      */
//     double betaC = sqrt(gammaz * gammaz - 1.0) / gammaz / Physics::c;
// 
//     bunch_mp->Bf(0) = bunch_mp->Bf(0) - betaC * bunch_mp->Ef(1);
//     bunch_mp->Bf(1) = bunch_mp->Bf(1) + betaC * bunch_mp->Ef(0);
    throw OpalException("AmrBoxLib::computeSelfFields() ", "Not yet Implemented.");
}


void AmrBoxLib::computeSelfFields(int bin) {
    //TODO
    throw OpalException("AmrBoxLib::computeSelfFields(int bin) ", "Not yet Implemented.");
}


void AmrBoxLib::computeSelfFields_cycl(double gamma) {    
//     for (int i = 0; i <= finest_level; ++i) {
//         std::cout << "Level " << i << ": " << grids[i] << std::endl;
//         std::cout << "Level " << i << ": " << layout_mp->ParticleBoxArray(i) << std::endl;
//     }
    
//     std::cout << "regrid done" << std::endl;
    
//     bunch_mp->python_format(0); std::cout << "Written." << std::endl; std::cin.get();
    
    //FIXME Lorentz transformation
    
    
    
    //scatter charges onto grid
//     bunch_mp->Q *= bunch_mp->dt;
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    
    /// from charge (C) to charge density (C/m^3).
    amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, finest_level);
//     bunch_mp->Q /= bunch_mp->dt;
    int baseLevel = 0;
    int nLevel = finest_level + 1;
    double invGamma = 1.0 / gamma /* bunch_mp->getCouplingConstant()*/;
    double scalefactor = 1.0; //bunch_mp->getdT();
    
//     std::cout << "Coupling: " << bunch_mp->getCouplingConstant() << std::endl; std::cin.get();
    
//     std::cout << "finest_level = " << finest_level << std::endl;
//     std::cout << "finest_level = " << (&amrpbase_p->getAmrLayout())->finestLevel() << std::endl; std::cin.get();
    
    /// Lorentz transformation
    /// In particle rest frame, the longitudinal length (y for cyclotron) enlarged
    for (int i = 0; i <= finest_level; ++i) {
        this->rho_m[i]->mult(invGamma / scalefactor, 0, 1);
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
    for (MFIter mfi(*rho_m[level]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const FArrayBox& fab = (*rho_m[level])[mfi];
        
        for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
            for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                    IntVect iv(x, y, z);
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
        this->rho_m[i]->mult(-1.0 / Physics::epsilon_0 / scalefactor, 0, 1);
    }
    
    
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    IpplTimings::startTimer(bunch_mp->compPotenTimer_m);
    solver->solve(rho_m, phi_m, eg_m, baseLevel, finest_level);
    IpplTimings::stopTimer(bunch_mp->compPotenTimer_m);
    
#ifdef DBG_SCALARFIELD
    INFOMSG("*** START DUMPING SCALAR FIELD ***" << endl);
    std::ofstream fstr2;
    fstr2.precision(9);

    std::string phi_fn = std::string("data/") + SfileName + std::string("-phi_scalar-") + std::to_string(fieldDBGStep_m);
    fstr2.open(phi_fn.c_str(), std::ios::out);
    
    for (MFIter mfi(*phi_m[level]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const FArrayBox& fab = (*phi_m[level])[mfi];
        
        for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
            for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                    IntVect iv(x, y, z);
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
        //ostringstream oss;
        //MPI_File file;
        //MPI_Status status;
        //MPI_Info fileinfo;
        //MPI_File_open(Ippl::getComm(), "rho_scalar", MPI_MODE_WRONLY | MPI_MODE_CREATE, fileinfo, &file);
        std::ofstream fstr;
        fstr.precision(9);

        std::string e_field = std::string("data/") + SfileName + std::string("-e_field-") + std::to_string(fieldDBGStep_m);
        fstr.open(e_field.c_str(), std::ios::out);
        
        for (MFIter mfi(*eg_m[level]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const FArrayBox& fab = (*eg_m[level])[mfi];
            
            for (int x = bx.loVect()[0]; x <= bx.hiVect()[0]; ++x) {
                for (int y = bx.loVect()[1]; y <= bx.hiVect()[1]; ++y) {
                    for (int z = bx.loVect()[2]; z <= bx.hiVect()[2]; ++z) {
                        IntVect iv(x, y, z);
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
    
    bunch_mp->Ef *= Vector_t(gamma * scalefactor, invGamma * scalefactor, gamma * scalefactor);
    
    /// calculate coefficient
    // Relativistic E&M says gamma*v/c^2 = gamma*beta/c = sqrt(gamma*gamma-1)/c
    // but because we already transformed E_trans into the moving frame we have to
    // add 1/gamma so we are using the E_trans from the rest frame -DW
    double betaC = std::sqrt(gamma * gamma - 1.0) * invGamma / Physics::c;
    
    /// calculate B field from E field
    bunch_mp->Bf(0) =  betaC * bunch_mp->Ef(2);
    bunch_mp->Bf(2) = -betaC * bunch_mp->Ef(0);
    
//     std::cout << "self-field computed" << std::endl;
}


void AmrBoxLib::computeSelfFields_cycl(int bin) {
    //TODO
    throw OpalException("AmrBoxLib::computeSelfFields_cycl(int bin) ", "Not yet Implemented.");
}


void AmrBoxLib::updateMesh() {
    //FIXME What about resizing mesh, i.e. geometry?
    const AmrGeomContainer_t& geom = this->Geom();
    
    const Real* tmp = geom[0].CellSize();
    
    Vector_t hr;
    for (int i = 0; i < 3; ++i)
        hr[i] = tmp[i];
    
    bunch_mp->setBaseLevelMeshSpacing(hr);
}

Vektor<int, 3> AmrBoxLib::getBaseLevelGridPoints() {
    const AmrGeomContainer_t& geom = this->Geom();
    const Box& bx = geom[0].Domain();
    
    const IntVect& low = bx.smallEnd();
    const IntVect& high = bx.bigEnd();
    
    return Vektor<int, 3>(high[0] - low[0] + 1,
                          high[1] - low[1] + 1,
                          high[2] - low[2] + 1);
}


int AmrBoxLib::maxLevel() {
    return AmrCore::maxLevel();
}


int AmrBoxLib::finestLevel() {
    return AmrCore::finestLevel();
}


void AmrBoxLib::updateBunch()  {
    const Array<Geometry>& geom = this->Geom();
    const Array<DistributionMapping>& dmap = this->DistributionMap();
    const Array<BoxArray>& ba = this->boxArray();
    const Array<IntVect>& ref_rato = this->refRatio ();
        
        
    Array<int> rr( ref_rato.size() );
    for (unsigned int i = 0; i < rr.size(); ++i) {
        rr[i] = ref_rato[i][0];
    }        
    
    layout_mp->define(geom, ba, dmap, rr);
        
    bunch_mp->update();
        
//     this->setBunch(bunch_mp);
}


void AmrBoxLib::RemakeLevel (int lev, Real time,
                             const BoxArray& new_grids,
                             const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m[lev].reset(new AmrField_t(new_grids, 1, 1, new_dmap));
    nChargePerCell_m[lev]->setVal(0.0);
    
    rho_m[lev].reset(new AmrField_t(new_grids, 1, 1, new_dmap));
    rho_m[lev]->setVal(0.0);
    
    /*
     * particles need to know the BoxArray
     * and DistributionMapping
     */
     layout_mp->SetParticleBoxArray(lev, new_grids);
     layout_mp->SetParticleDistributionMap(lev, new_dmap);
}


void AmrBoxLib::MakeNewLevel (int lev, Real time,
                              const BoxArray& new_grids,
                              const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m[lev].reset(new AmrField_t(new_grids, 1, 1, new_dmap));
    nChargePerCell_m[lev]->setVal(0.0);
    
    rho_m[lev].reset(new AmrField_t(new_grids, 1, 1, new_dmap));
    rho_m[lev]->setVal(0.0);
    
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


void AmrBoxLib::ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow) {
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


void AmrBoxLib::tagForChargeDensity_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    std::cout << "-----------------------" << std::endl
              << "Tagging." << std::endl;
              
    std::cout << "  Finest Level: " << finest_level << std::endl
              << "  Level:        " << lev << std::endl;
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    for (int i = lev; i <= finest_level; ++i) {
        nChargePerCell_m[i]->setVal(0.0);
        // single level scatter
        amrpbase_p->scatter(bunch_mp->Q, *nChargePerCell_m[i], bunch_mp->R, i);
    }
    
    for (int i = finest_level-1; i >= lev; --i) {
        AmrField_t tmp(nChargePerCell_m[i]->boxArray(), 1, 0, nChargePerCell_m[i]->DistributionMap());
        tmp.setVal(0.0);
        BoxLib::average_down(*nChargePerCell_m[i+1], tmp, 0, 1, refRatio(i));
        AmrField_t::Add(*nChargePerCell_m[i], tmp, 0, 0, 1, 0);
    }
    
    for (int i = lev; i <= finest_level; ++i) {
        std::cout << "level " << i << " sum = " << nChargePerCell_m[i]->sum() << " min = "
                  << nChargePerCell_m[i]->min(0) << " max = " << nChargePerCell_m[i]->max(0) << std::endl;
    }
    
    /* BoxLib stores charge per cell volume, we thus need to 
     * divide by the charge per level by the cell volume of this level
     */
    double cellVolume = 1.0;
    for (int i = 0; i < BL_SPACEDIM; ++i) {        
        double dx = geom[lev].CellSize(i);
        cellVolume *= dx;
    }
    double nCharge = nCharge_m / cellVolume;
    
    std::cout << "Cell volume = " << cellVolume << std::endl
              << "tag for charge = " << nCharge_m << std::endl
              << "nCharge / cell volume = " << nCharge << std::endl;
    
    
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        for (MFIter mfi(*nChargePerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
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
                        BL_TO_FORTRAN_3D((*nChargePerCell_m[lev])[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &nCharge);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
    
    std::cout << "Tagging done." << std::endl
              << "----------------------" << std::endl;
}


void AmrBoxLib::tagForPotentialStrength_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the potential is higher than 75 percent of the maximum potential
     * value of this level.
     */
    
    // 1. Assign charge onto grid of level lev
    int baseLevel   = 0;
    int finestLevel = 0;
    
    AmrFieldContainer_t nPartPerCell(1);
    nPartPerCell[9] = std::unique_ptr<AmrField_t>(new AmrField_t(this->boxArray(lev), 1, 1,
                                                  this->DistributionMap(lev)));
    nPartPerCell[baseLevel]->setVal(0.0);
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    // single level scatter
    amrpbase_p->scatter(bunch_mp->Q, *nPartPerCell[baseLevel], bunch_mp->R, lev);
    
    if ( geom[0].isAllPeriodic() ) {
        Real offset = 0.0;
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
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    
    Real minp = 0.0;
    Real maxp = 0.0;
    maxp = scaling_m * phi[baseLevel]->max(0);
    minp = scaling_m * phi[baseLevel]->min(0);
    Real threshold = std::max( std::abs(minp), std::abs(maxp) );
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        for (MFIter mfi(*phi[baseLevel],false); mfi.isValid(); ++mfi) {
            const Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            TagBox&     tagfab  = tags[mfi];
            
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


void AmrBoxLib::tagForEfield_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the efield is higher than 75 percent of the maximum efield
     * value of this level.
     */
    
    // 1. Assign charge onto grid of level lev
    int baseLevel   = 0;
    int finestLevel = 0;
    
    AmrFieldContainer_t nPartPerCell(1);
    nPartPerCell[baseLevel] = std::unique_ptr<AmrField_t>(new AmrField_t(this->boxArray(lev), 1, 1,
                                                          this->DistributionMap(lev)));
    nPartPerCell[baseLevel]->setVal(0.0);
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    // single level scatter
    amrpbase_p->scatter(bunch_mp->Q, *nPartPerCell[baseLevel], bunch_mp->R, lev);
    
    if ( geom[0].isAllPeriodic() ) {
        Real offset = 0.0;
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
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;
    
    Real min[3] = {0.0, 0.0, 0.0};
    Real max[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        max[i] = scaling_m * efield[baseLevel]->max(i);
        min[i] = scaling_m * efield[baseLevel]->min(i);
    }
    Real threshold[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)
        threshold[i] = std::max( std::abs(min[i]), std::abs(max[i]) );
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        // mfi(efield[baseLevel], true)
        for (MFIter mfi(*efield[baseLevel],false); mfi.isValid(); ++mfi) {
            const Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            TagBox&     tagfab  = tags[mfi];
            FArrayBox&  fab     = (*efield[baseLevel])[mfi];
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
                        IntVect iv(i,j,k);
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
    
    IntVect low(0, 0, 0); 
    IntVect high(nGridPts[0] - 1,
                 nGridPts[1] - 1,
                 nGridPts[2] - 1);
    
    const Box bx(low, high);
    BoxArray ba(bx);
    ba.maxSize( this->maxGridSize(0) );
//     this->SetBoxArray(0, ba);
    
    DistributionMapping dmap;
    dmap.define(ba, ParallelDescriptor::NProcs());
//     this->SetDistributionMap(0, dmap);
    
    this->RemakeLevel(0, 0.0, ba, dmap);
    
    for (int i = 0; i < maxLevel() + 1; ++i) {
        std::cout << "level " << i << ": " << this->Geom(i) << std::endl;
    }
    
    layout_mp->define(this->Geom());
    
//     const BoxArray& ba = layout_mp->ParticleBoxArray(0);
//     const DistributionMapping& dm = layout_mp->ParticleDistributionMap(0 /*level*/);
    
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

