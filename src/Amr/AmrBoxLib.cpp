#include "AmrBoxLib.h"

#include "Algorithms/AmrPartBunch.h"
#include "Structure/FieldSolver.h"
#include "Solvers/PoissonSolver.h"

#ifdef AMR_YT_DUMP
    #include "Amr/AmrYtWriter.h"
#endif

#ifdef AMR_PYTHON_DUMP
    #include "Amr/AmrPythonWriter.h"
#endif

#ifdef DBG_SCALARFIELD
    #include "Amr/AmrSliceWriter.h"
#endif
    
#include "AmrBoxLib_F.h"
#include <AMReX_MultiFabUtil.H>

#include <AMReX_ParmParse.H> // used in initialize function

extern Inform* gmsg;

AmrBoxLib::AmrBoxLib() : AmrObject(),
                         amrex::AmrMesh(),
                         bunch_mp(nullptr),
                         layout_mp(nullptr),
                         rho_m(0),
                         phi_m(0),
                         efield_m(0)
{}


AmrBoxLib::AmrBoxLib(TaggingCriteria tagging,
                     double scaling,
                     double nCharge)
    : AmrObject(tagging, scaling, nCharge),
      amrex::AmrMesh(),
      bunch_mp(nullptr),
      layout_mp(nullptr),
      rho_m(0),
      phi_m(0),
      efield_m(0)
{}



AmrBoxLib::AmrBoxLib(const AmrDomain_t& domain,
                     const AmrIntArray_t& nGridPts,
                     short maxLevel)
    : AmrObject(),
      amrex::AmrMesh(&domain, maxLevel, nGridPts, 0 /* cartesian */),
      bunch_mp(nullptr),
      layout_mp(nullptr),
      rho_m(maxLevel + 1),
      phi_m(maxLevel + 1),
      efield_m(maxLevel + 1)
{}


AmrBoxLib::AmrBoxLib(const AmrDomain_t& domain,
                     const AmrIntArray_t& nGridPts,
                     short maxLevel,
                     AmrPartBunch* bunch_p)
    : AmrObject(),
      amrex::AmrMesh(&domain, maxLevel, nGridPts, 0 /* cartesian */),
      bunch_mp(bunch_p),
      layout_mp(static_cast<AmrLayout_t*>(&bunch_p->getLayout())),
      rho_m(maxLevel + 1),
      phi_m(maxLevel + 1),
      efield_m(maxLevel + 1)
{
    /*
     * The layout needs to know how many levels we can make.
     */
    layout_mp->resize(maxLevel);
    
    initBaseLevel_m(nGridPts);
    
    // set mesh spacing of bunch
    updateMesh();
}


std::unique_ptr<AmrBoxLib> AmrBoxLib::create(const AmrInitialInfo& info,
                                             AmrPartBunch* bunch_p)
{
    /* The bunch is initialized first with a Geometry,
     * BoxArray and DistributionMapping on
     * the base level (level = 0). Thus, we take the domain specified there in
     * order to create the Amr object.
     */
    AmrLayout_t* layout_p = static_cast<AmrLayout_t*>(&bunch_p->getLayout());
    AmrDomain_t domain = layout_p->Geom(0).ProbDomain();
    
    AmrIntArray_t nGridPts = {
        info.gridx,
        info.gridy,
        info.gridz
    };
    
    int maxlevel = info.maxlevel;
    
    /*
     * further attributes are given by the BoxLib's ParmParse class.
     */
    initParmParse_m(info, layout_p);
    
    return std::unique_ptr<AmrBoxLib>(new AmrBoxLib(domain,
                                                    nGridPts,
                                                    maxlevel,
                                                    bunch_p
                                                   )
                                     );
}


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
    for (unsigned int lev = 0; lev < efield_m.size(); ++lev) {
        for (std::size_t i = 0; i < bunch_mp->Dimension; ++i) {
            /* calls:
             *  - max(comp, nghost (to search), local)
             *  - min(cmop, nghost, local)
             */
            double max = efield_m[lev]->max(i, 0, false);
            maxE[i] = (maxE[i] < max) ? max : maxE[i];
            
            double min = efield_m[lev]->min(i, 0, false);
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
    
    int baseLevel = 0;
    const int& finestLevel = finest_level; //(&amrpbase_p->getAmrLayout())->finestLevel();
    
    int nLevel = finestLevel + 1;
    rho_m.resize(nLevel);
    phi_m.resize(nLevel);
    efield_m.resize(nLevel);
    
    /// from charge (C) to charge density (C/m^3).
    bunch_mp->Q *= bunch_mp->dt;
    amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, finest_level);
    bunch_mp->Q /= bunch_mp->dt;
    
    //calculating mesh-scale factor
    double gammaz = sum(bunch_mp->P)[2] / bunch_mp->getTotalNum();
    gammaz *= gammaz;
    gammaz = std::sqrt(gammaz + 1.0);
    
    double tmp = 1.0 / bunch_mp->getdT() / gammaz * scalefactor;
    for (int i = 0; i <= finestLevel; ++i) {
        this->rho_m[i]->mult(tmp, 0, 1);
    }
    
    
    for (int i = 0; i <= finestLevel; ++i)
        this->rho_m[i]->mult(-1.0 / Physics::epsilon_0, 0, 1);
    
    // charge density is in rho_m
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    
    IpplTimings::startTimer(bunch_mp->compPotenTimer_m);
    solver->solve(rho_m, phi_m, efield_m, baseLevel, finestLevel);
    IpplTimings::stopTimer(bunch_mp->compPotenTimer_m);
    
    // apply scale of electric-field in order to undo the transformation
    for (int i = 0; i <= finestLevel; ++i)
        this->efield_m[i]->mult(scalefactor, 0, 3);
    
    amrpbase_p->gather(bunch_mp->Ef, this->efield_m, bunch_mp->R, 0, finest_level);
    
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
    
    /*
     * dumping only
     */

#ifdef AMR_YT_DUMP
    INFOMSG("*** START DUMPING FIELDS IN YT FORMAT ***" << endl);
    AmrYtWriter ytWriter(bunch_mp->getLocalTrackStep());
    
    AmrIntArray_t rr(nLevel);
    for (int i = 0; i < nLevel - 1; ++i)
        rr[i] = this->MaxRefRatio(i);
    
    double time = bunch_mp->getT(); // in seconds
    
    // we need to undo coefficient when writing charge density
    for (int i = 0; i <= finest_level; ++i)
        this->rho_m[i]->mult(- Physics::epsilon_0 / scalefactor, 0, 1);
    
    ytWriter.writeFields(rho_m, phi_m, efield_m, rr, this->Geom(), time, scalefactor);
    INFOMSG("*** FINISHED DUMPING FIELDS IN YT FORMAT ***" << endl);
#endif

#ifdef AMR_PYTHON_DUMP
    INFOMSG("*** START DUMPING BUNCH AND GRIDS IN PYTHON FORMAT ***" << endl);
    AmrPythonWriter pyWriter;
    pyWriter.writeBunch(bunch_mp, scalefactor);
    INFOMSG("*** FINISHED DUMPING BUNCH AND GRIDS IN PYTHON FORMAT ***" << endl);
#endif
    
#ifdef DBG_SCALARFIELD
    if ( Ippl::getNodes() > 1 )
        throw OpalException("AmrBoxLib::computeSelfFields() ",
                            "Dumping only in serial execution.");
    
    int step = bunch_mp->getLocalTrackStep();
    AmrSliceWriter sliceWriter;
    
#ifdef AMR_YT_DUMP
    // make sure we undo only once if AMR_YT_DUMP is also enabled
#else
    // we need to undo coefficient when writing charge density
    for (int i = 0; i <= finest_level; ++i)
        this->rho_m[i]->mult(- Physics::epsilon_0 / scalefactor, 0, 1);
#endif
    
    sliceWriter.writeFields(rho_m, phi_m, efield_m,
                            AmrIntArray_t(), this->Geom(), step);
#endif
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
    double scalefactor = layout_mp->domainMapping(*amrpbase_p);
    
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
    
    
    // charge density is in rho_m
    // calculate Possion equation (with coefficient: -1/(eps))
    for (int i = 0; i <= finest_level; ++i) {
        this->rho_m[i]->mult(-1.0 / Physics::epsilon_0, 0, 1);
    }
    
    
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    IpplTimings::startTimer(bunch_mp->compPotenTimer_m);
    solver->solve(rho_m, phi_m, efield_m, baseLevel, finest_level);
    IpplTimings::stopTimer(bunch_mp->compPotenTimer_m);
    
    // apply scale of electric-field in order to undo the transformation
    for (int i = 0; i <= finestLevel(); ++i)
        this->efield_m[i]->mult(scalefactor, 0, 3);
    
    for (int i = 0; i <= finest_level; ++i) {
        if ( this->efield_m[i]->contains_nan(false) )
            throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ",
                                "Ef: NANs at level " + std::to_string(i) + ".");
    }
    
    
    amrpbase_p->gather(bunch_mp->Ef, this->efield_m, bunch_mp->R, 0, finest_level);
    
    // undo domain change
    layout_mp->domainMapping(*amrpbase_p, true);
    
    /// Back Lorentz transformation
    bunch_mp->Ef *= Vector_t(gamma,
                             invGamma,
                             gamma);
    
    /// calculate coefficient
    // Relativistic E&M says gamma*v/c^2 = gamma*beta/c = sqrt(gamma*gamma-1)/c
    // but because we already transformed E_trans into the moving frame we have to
    // add 1/gamma so we are using the E_trans from the rest frame -DW
    double betaC = std::sqrt(gamma * gamma - 1.0) * invGamma / Physics::c;
    
    /// calculate B field from E field
    bunch_mp->Bf(0) =  betaC * bunch_mp->Ef(2);
    bunch_mp->Bf(2) = -betaC * bunch_mp->Ef(0);
    
    /*
     * dumping only
     */
    
#ifdef AMR_YT_DUMP
    INFOMSG("*** START DUMPING FIELDS IN YT FORMAT ***" << endl);
    AmrYtWriter ytWriter(bunch_mp->getLocalTrackStep());
    
    AmrIntArray_t rr(nLevel);
    for (int i = 0; i < nLevel - 1; ++i)
        rr[i] = this->MaxRefRatio(i);
    
    double time = bunch_mp->getT(); // ps
    
    // we need to undo coefficient when writing charge density
    for (int i = 0; i <= finest_level; ++i)
        this->rho_m[i]->mult(- Physics::epsilon_0 / scalefactor, 0, 1);
    
    ytWriter.writeFields(rho_m, phi_m, efield_m, rr, this->Geom(), time, scalefactor);
    INFOMSG("*** FINISHED DUMPING FIELDS IN YT FORMAT ***" << endl);
#endif

#ifdef AMR_PYTHON_DUMP
    INFOMSG("*** START DUMPING BUNCH AND GRIDS IN PYTHON FORMAT ***" << endl);
    AmrPythonWriter pyWriter;
    pyWriter.writeBunch(bunch_mp, scalefactor);
    INFOMSG("*** FINISHED DUMPING BUNCH AND GRIDS IN PYTHON FORMAT ***" << endl);
#endif
    
#ifdef DBG_SCALARFIELD
    if ( Ippl::getNodes() > 1 )
        throw OpalException("AmrBoxLib::computeSelfFields_cycl(double gamma) ",
                            "Dumping only in serial execution.");
    
    int step = bunch_mp->getLocalTrackStep();
    AmrSliceWriter sliceWriter;
    
#ifdef AMR_YT_DUMP
    // make sure we undo only once if AMR_YT_DUMP is also enabled
#else
    // we need to undo coefficient when writing charge density
    for (int i = 0; i <= finest_level; ++i)
        this->rho_m[i]->mult(- Physics::epsilon_0 / scalefactor, 0, 1);
#endif
    
    sliceWriter.writeFields(rho_m, phi_m, efield_m,
                            AmrIntArray_t(), this->Geom(), step);
#endif
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


void AmrBoxLib::RemakeLevel (int lev, AmrReal_t time,
                             const AmrGrid_t& new_grids,
                             const AmrProcMap_t& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    //                                                      #comp  #ghosts cells
    rho_m[lev].reset(new AmrField_t(new_grids, new_dmap,    1,     1));
    phi_m[lev].reset(new AmrField_t(new_grids, new_dmap,    1,     1));
    efield_m[lev].reset(new AmrField_t(new_grids, new_dmap, 3,     1));
    
    // including nghost = 1
    rho_m[lev]->setVal(0.0, 1);
    phi_m[lev]->setVal(0.0, 1);
    efield_m[lev]->setVal(0.0, 1);
    
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
    
    //                                                      #comp  #ghosts cells
    rho_m[lev].reset(new AmrField_t(new_grids, new_dmap,    1,     1));
    phi_m[lev].reset(new AmrField_t(new_grids, new_dmap,    1,     1));
    efield_m[lev].reset(new AmrField_t(new_grids, new_dmap, 3,     1));
    
    // including nghost = 1
    rho_m[lev]->setVal(0.0, 1);
    phi_m[lev]->setVal(0.0, 1);
    efield_m[lev]->setVal(0.0, 1);
    
    
    
    /*
     * particles need to know the BoxArray
     * and DistributionMapping
     */
    layout_mp->SetParticleBoxArray(lev, new_grids);
    layout_mp->SetParticleDistributionMap(lev, new_dmap);
}


void AmrBoxLib::ClearLevel(int lev) {
    rho_m[lev].reset(nullptr);
    phi_m[lev].reset(nullptr);
    efield_m[lev].reset(nullptr);
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
        rho_m[i]->setVal(0.0, 1);
    
    // the new scatter function averages the value also down to the coarsest level
    amrpbase_p->scatter(bunch_mp->Q, rho_m,
                        bunch_mp->R, lev, finest_level);
    
    const double& scalefactor = (&amrpbase_p->getAmrLayout())->getScalingFactor();
    
    for (int i = lev; i <= finest_level; ++i)
        rho_m[i]->mult(scalefactor, 0, 1);
    
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;

    const AmrReal_t* dx      = geom[lev].CellSize();
    const AmrReal_t* prob_lo = geom[lev].ProbLo();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        AmrIntArray_t  itags;
        for (amrex::MFIter mfi(*rho_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
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
                        BL_TO_FORTRAN_3D((*rho_m[lev])[mfi]),
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
    if ( !phi_m[lev]->ok() || (time == 0 && !(phi_m[lev]->norm0(0) > 0)) ) {
        *gmsg << "We need to perform charge tagging in the first time step" << endl;
        this->tagForChargeDensity_m(lev, tags, time, ngrow);
        
        return;
    }
    
    // tag cells for refinement
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;

    const AmrReal_t* dx      = geom[lev].CellSize();
    const AmrReal_t* prob_lo = geom[lev].ProbLo();
    
    AmrReal_t threshold = phi_m[lev]->norm0(0, 0 /*nghost*/, false /*local*/);
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        AmrIntArray_t  itags;
        for (amrex::MFIter mfi(*phi_m[lev],false); mfi.isValid(); ++mfi) {
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
                        BL_TO_FORTRAN_3D((*phi_m[lev])[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &threshold);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
}


void AmrBoxLib::tagForEfield_m(int lev, TagBoxArray_t& tags, AmrReal_t time, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the efield is higher than 75 percent of the maximum efield
     * value of this level.
     */
    
    if ( !efield_m[lev]->ok() || (time == 0 && !(efield_m[lev]->norm0(0) > 0)) ) {
        *gmsg << "We need to perform charge tagging in the first time step" << endl;
        this->tagForChargeDensity_m(lev, tags, time, ngrow);
        
        return;
    }
    
    // tag cells for refinement
    const int clearval = amrex::TagBox::CLEAR;
    const int   tagval = amrex::TagBox::SET;
    
    // obtain maximum absolute value
    AmrIntArray_t comps = {0, 1, 2};
    amrex::Array<AmrReal_t> threshold = efield_m[lev]->norm0(comps,
                                                             0 /*nghosts*/,
                                                             false /*local*/);
    
    for (int i = 0; i < 3; ++i)
        threshold[i] *= scaling_m;
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        AmrIntArray_t  itags;
        // mfi(efield_m[baseLevel], true)
        for (amrex::MFIter mfi(*efield_m[lev],false); mfi.isValid(); ++mfi) {
            const amrex::Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            amrex::TagBox&     tagfab  = tags[mfi];
            amrex::FArrayBox&  fab     = (*efield_m[lev])[mfi];
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
    
    AmrProcMap_t dmap;
    dmap.define(ba, amrex::ParallelDescriptor::NProcs());
    
    this->RemakeLevel(0, 0.0, ba, dmap);
    
    layout_mp->define(this->Geom());
    layout_mp->define(this->ref_ratio);    
}


void AmrBoxLib::initParmParse_m(const AmrInitialInfo& info, AmrLayout_t* layout_p) {
    
    /*
     * All parameters that we set with the OPAL input file
     */
    amrex::ParmParse pAmr("amr");
    pAmr.add("max_grid_size", info.maxgrid);
    
    const int nratios_vect = info.maxlevel*BL_SPACEDIM;
    
    AmrIntArray_t refRatio(nratios_vect);
    
    for (int i = 0; i < info.maxlevel; ++i) {
        refRatio[i * BL_SPACEDIM]     = info.refratx;
        refRatio[i * BL_SPACEDIM + 1] = info.refraty;
        refRatio[i * BL_SPACEDIM + 2] = info.refratz;
    }
    
    pAmr.addarr("ref_ratio_vect", refRatio);
    
    amrex::ParmParse pGeom("geometry");
    amrex::Array<int> isPeriodic = {
        layout_p->Geom(0).isPeriodic(0),
        layout_p->Geom(0).isPeriodic(1),
        layout_p->Geom(0).isPeriodic(2)
    };
    pGeom.addarr("is_periodic", isPeriodic);
    
    
    /*
     * "All" other parameters, we moslty take the
     * defaults of the code
     * 
     * ATTENTION Be careful with default values since they might
     *           be changed by the AMReX developers!
     * 
     * Parmparse parameters not to be set because
     * alreday given due to constructor
     * ----------------------------------
     * 
     *  AMReX_AmrMesh:
     *      - amr.max_level
     * 
     *  AMReX_Geometry:
     *      - geometry.coord_sys
     *      - geometry.prob_lo
     *      - geometry.prob_hi
     *      - geometry.is_periodic
     *      - geometry.spherical_origin_fix [default: 0]
     *        (not set because we're using Cartesian coordinates)
     * 
     * 
     * ParmParse left out
     * ------------------
     *  AMReX_FabArrayBase:
     *      - fabarray.mfiter_tile_size [default: IntVect(1024000,8,8)]
     *      - fabarray.mfghostiter_tile_size [default: IntVect(1024000, 8, 8)]
     *      - fabarray.comm_tile_size [default: IntVect(1024000, 8, 8)]
     *      - fabarray.maxcomp [default: 25]
     *      - fabarray.do_async_sends [default: true]
     * 
     *  AMReX_MemPool:
     *      - fab.init_snan [default: 0, if not BL_TESTING or DEBUG enabled]
     * 
     *  AMReX_MemProfiler:
     *      - amrex.memory.log [default: "memlog"]
     * 
     *  AMReX_VisMF:
     *      - vismf.v [verbose, default: 0]
     *      - vismf.headerversion [default: VisMF::Header::Version_v1 = 1]
     *      - vismf.groupsets [default: false]
     *      - vismf.setbuf [default: true]
     *      - vismf.usesingleread [default: false]
     *      - vismf.usesinglewrite [default: false]
     *      - vismf.checkfilepositions [default: false]
     *      - vismf.usepersistentifstreams [default: true]
     *      - vismf.usesynchronousreads [default: false]
     *      - vismf.usedynamicsetselection [default: true]
     *      - vismf.iobuffersize [default: VisMF::IO_Buffer_Size = 262144 * 8]
     * 
     *  AMReX_ParallelDescriptor:
     *      Needs change of ParallelDescriptor::SetNProcsSidecars in the function
     *      ParallelDescriptor::StartParallel()
     *      - amrex.ncolors [default: m_nCommColors = 1]
     *      - team.size     [default: disabled, enable by BL_USE_UPCXX or BL_USE_MPI3]
     *      - team.reduce   [default: disabled, enable by BL_USE_UPCXX or BL_USE_MPI3]
     *      - mpi.onesided  [default: disabled, enable by BL_USE_UPCXX or BL_USE_MPI3]
     * 
     *  AMReX:
     *      Debugging purposes
     *      - amrex.v [verbose, default: 0]
     *      - amrex.verbose [default: 0]
     *      - amrex.fpe_trap_invalid [invalid, default: 0]
     *      - amrex.fpe_trap_zero [divbyzero, default: 0]
     *      - amrex.fpe_trap_overflow [overflow, 0]
     * 
     *  AMReX_FArrayBox:
     *      For output only
     *      - fab.format [default: FABio::FAB_NATIVE]
     *      For reading in an old FAB
     *      - fab.ordering
     *      - fab.initval [default: quiet_NaN() or max()]
     *      - fab.do_initval [default: false, if not DEBUG or BL_TESTING enabled]
     *      - fab.init_snan [default: false, if not DEBUG or BL_TESTING enabled]
     * 
     *  AMReX_MCMultiGrid:
     *      We do NOT use this solver
     *      - mg.maxiter [default: 40]
     *      - mg.numiter [default: -1]
     *      - mg.nu_0 [default: 1]
     *      - mg.nu_1 [default: 2]
     *      - mg.nu_2 [default: 2]
     *      - mg.nu_f [default: 8]
     *      - mg.v [verbose, default: 0]
     *      - mg.usecg [default: 1]
     *      - mg.rtol_b [default: 0.01]
     *      - mg.bot_atol [default: -1.0]
     *      - mg.nu_b [default: 0]
     *      - mg.numLevelsMAX [default: 1024]
     * 
     *  AMReX_MCLinOp:
     *      We do NOT use this solver
     *      - MCLp.harmavg [default: 0]
     *      - MCLp.v [verbose, default: 0]
     *      - MCLp.maxorder [default: 2]
     * 
     *  AMReX_MCCGSolver:
     *      We do NOT use this solver
     *      - cg.maxiter [default: 40]
     *      - cg.v [verbose, default: 0]
     *      - cg.isExpert [default: 0]
     *      - MCCGSolver::def_unstable_criterion [default: 10]
     * 
     *  AMReX_LinOp:
     *      We do NOT use this class
     *      - Lp.harmavg [default: 0]
     *      - Lp.v [verbose, default: 0]
     *      - Lp.maxorder [default: 2]
     *      - LinOp::LinOp_grow [default: 1]
     * 
     *  AMReX_MultiGrid (single level solver):
     *      We do NOT use this solver, see BoxLibSolvers/FMGPoissonSolver for
     *      explanation of parameters
     *      - mg.v [verbose, default: 0]
     *      - mg.nu_0 [default: 1]
     *      - mg.nu_1 [default: 2]
     *      - mg.nu_2 [default: 2]
     *      - mg.nu_f [default: 8]
     *      - mg.nu_b [default: 0]
     *      - mg.usecg [default: 1]
     *      - mg.rtol_b [default: 0.0001 or 0.01 (if CG_USE_OLD_CONVERGENCE_CRITERIA enabled)]
     *      - mg.verbose [default: 0]
     *      - mg.maxiter [default: 40]
     *      - mg.bot_atol [default: -1.0]
     *      - mg.maxiter_b [default: 120]
     *      - mg.numLevelsMAX [default: 1024]
     *      - mg.smooth_on_cg_unstable [default: 1]
     *      - mg.use_Anorm_for_convergence [default: 1]
     * 
     *  AMReX_CGSolver (single level solver):
     *      We do NOT use this solver
     *      - cg.v [verbose, default: 0]
     *      - cg.SSS [default: SSS_MAX = 4]
     *      - cg.maxiter [default: 80]
     *      - cg.verbose [default: 0]
     *      - cg.variable_SSS [default: true]
     *      - cg.use_jbb_precond [default: 0]
     *      - cg.use_jacobi_precond [default: 0]
     *      - cg.unstable_criterion [default: 10]
     *      - CGSolver::def_cg_solver [default: BiCGStab]
     *      - cg.cg_solver [default: CGSolver::def_cg_solver]
     * 
     *  AMReX_Amr:
     *      We do NOT use this class
     *      - amr.regrid_on_restart [default: 0]
     *      - amr.use_efficient_regrid [default: 0]
     *      - amr.plotfile_on_restart [default: 0]
     *      - amr.checkpoint_on_restart [default: 0]
     *      - amr.compute_new_dt_on_regrid [default: 0]
     *      - amr.mffile_nstreams [default: 1]
     *      - amr.probinit_natonce [default: 32]
     *      - amr.file_name_digits [default: 5]
     *      - amr.initial_grid_file
     *      - amr.regrid_file
     *      - amr.message_int [default: 10]
     *      - amr.run_log [default: not parsed]
     *      - amr.run_log_terse [default: not parsed]
     *      - amr.grid_log [default: not parsed]
     *      - amr.data_log [default: not parsed]
     *      - amr.probin_file [default: not parsed]
     *      - amr.restart
     *      - amr.restart_from_plotfile
     *      - amr.regrid_int [default: 1 at all levels]
     *      - amr.rebalance_grids [default: 0]
     *      - amr.nosub [default: not parsed]
     *      - amr.subcycling_mode [default: "Auto"]
     *      - amr.subcycling_iterations [default: 0]
     *      - amr.checkpoint_files_output
     *      - amr.plot_files_output
     *      - amr.plot_nfiles
     *      - amr.checkpoint_nfiles
     *      - amr.check_file [default: "chk"]
     *      - amr.check_int [default: -1]
     *      - amr.check_per [default: -1.0] 
     *      - amr.plot_file [default: "plt"]
     *      - amr.plot_int [default: -1]
     *      - amr.plot_per [default: -1.0]
     *      - amr.small_plot_file [default: "smallplt"]
     *      - amr.small_plot_int [default: -1]
     *      - amr.small_plot_per [default: -1.0]
     *      - amr.write_plotfile_with_checkpoint [default: 1]
     *      - amr.stream_max_tries [default: 4]
     *      - amr.abort_on_stream_retry_failure [default: false]
     *      - amr.precreateDirectories
     *      - amr.prereadFAHeaders
     *      - amr.plot_headerversion
     *      - amr.checkpoint_headerversion
     * 
     *  AMReX_SlabStat:
     *      We do NOT use this class
     *      - slabstat.boxes
     * 
     *  AMReX_AmrLevel:
     *      We do NOT use this class
     *      - amr.plot_vars [default: not parsed]
     *      - amr.derive_plot_vars [default: not parsed]
     *      - amr.small_plot_vars [default: not parsed]
     * 
     * AMReX_StationData:
     *      We do NOT use this class
     *      - StationData.vars [default: not parsed]
     *      - StationData.coord [default: not parsed]
     */
    
    // verbosity
    pAmr.add("v", 0);
    
    // # cells required for proper nesting.
    pAmr.add("n_proper", 1);
    
    // Grid efficiency. (default: 0.7)
    pAmr.add("grid_eff", 0.95);
    
    int nlev = info.maxlevel + 1;
    
    // Buffer cells around each tagged cell.
    amrex::Array<int> error_buf(nlev, 0);
    pAmr.addarr("n_error_buf", error_buf);  
    
    // Blocking factor in grid generation (by level).
    amrex::Array<int> blockingFactor(nlev, 8);
    
    // chop up grids to have more grids than the number of procs
    pAmr.add("refine_grid_layout", true);
    
    
    /*
     * ParmParse for DistributionMapping
     * 
     * round-robin: FAB i is owned by CPU i%N where N is total number of CPUs
     * knapsack:    FABs are partitioned across CPUs such that the total volume
     *              of the Boxes in the underlying BoxArray are as equal across
     *              CPUs as is possible.
     * SFC:         Is based on a space filling curve (default)
     */
    amrex::ParmParse pDmap("DistributionMapping");
    
    // verbosity
    pDmap.add("v", 0);
    
    // max_efficiency   = 0.9
    pDmap.add("efficiency", 0.9);
    
    // SFC = space filling curve
    pDmap.add("sfc_threshold", 0);
    
    /* if > 0:
     *  nworkers = node_size
     *  nteams   = nprocs / node_size
     *  
     *  if nwokers * nteams != nprocs:
     *      nteams   = nprocs
     *      nworkers = 1
     */
    pDmap.add("node_size", 0);
    
    /* Possibilities:
     *  - ROUNDROBIN
     *  - KNAPSACK
     *  - SFC
     *  - PFC
     *  - RRSFC
     */
    pDmap.add("strategy", "SFC");
}
