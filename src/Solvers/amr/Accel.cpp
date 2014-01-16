#ifdef HAVE_AMR_SOLVER
#include <iomanip>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <Utility.H>
#include <CONSTANTS.H>
#include <Accel.H>
#include <Accel_F.H>
#include <Derive_F.H>
#include <VisMF.H>
#include <TagBox.H>
#include <Particles_F.H>
#include "Electrostatic.H"

static Real fixed_dt    = -1.0;
static Real initial_dt  = -1.0;
static Real dt_cutoff   =  0;

bool Accel::dump_old = false;
int Accel::verbose      = 0;
int Accel::show_timings = 0;

Real Accel::init_shrink = 1.0;
Real Accel::change_max  = 1.1;

Real Accel::opal_c         = -1.0;
Real Accel::opal_coupling  = -1.0;

BCRec Accel::phys_bc;

ErrorList Accel::err_list;

Electrostatic* Accel::e_field_solver  =  0;

int Accel::allow_untagging    = 0;

#ifdef _OPENMP
#include <omp.h>
#endif

int Accel::write_parameters_in_plotfile = true;

// Note: Accel::variableSetUp is in Accel_setup.cpp
void
Accel::variable_cleanup ()
{
    if (e_field_solver != 0)
    {
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "Deleting e_field_solver in variable_cleanup...\n";
        delete e_field_solver;
        e_field_solver = 0;
    }

    desc_lst.clear();
}

void
Accel::read_params ()
{
    static bool done = false;

    if (done) return;  // (caseywstark) when would this happen?

    done = true;  // ?

    ParmParse pp("accel");

    pp.query("v", verbose);
    pp.query("show_timings", show_timings);
    pp.query("fixed_dt", fixed_dt);
    pp.query("initial_dt", initial_dt);
    pp.query("dump_old", dump_old);

    // Note that we *MUST* input these (hence the "get" instead of "query")
    pp.get("opal_c", opal_c);
    pp.get("opal_coupling", opal_coupling);

    // Get boundary conditions
    // Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    // pp.getarr("lo_bc", lo_bc, 0, BL_SPACEDIM);
    // pp.getarr("hi_bc", hi_bc, 0, BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i,0);
        phys_bc.setHi(i,0);
    }

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
#if 0
    if (Geometry::isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (Geometry::isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "Accel::read_params:periodic in direction "
                              << dir
                              << " but low BC is not Interior" << std::endl;
                    BoxLib::Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "Accel::read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior" << std::endl;
                    BoxLib::Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "Accel::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                BoxLib::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "Accel::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                BoxLib::Error();
            }
        }
    }
#endif

    pp.query("allow_untagging", allow_untagging);

    read_particle_params();
    
    pp.query("write_parameter_file",write_parameters_in_plotfile);
}

Accel::Accel ()
{
}

Accel::Accel (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,time)
{
    build_metrics();

    MultiFab& new_grav_mf = get_new_data(Elec_Field_Type);
    new_grav_mf.setVal(0);

    {
        // e_field_solver is a static object, only alloc if not already there
        if (e_field_solver == 0)
        e_field_solver = new Electrostatic(parent, parent->finestLevel(), &phys_bc);

        e_field_solver->install_level(level, this);
   }

    // Set field_n_grow to 3 on init. It'll be reset in advance.
    field_n_grow = 3;
}

Accel::~Accel ()
{
}

void
Accel::restart (Amr&     papa,
              istream& is,
              bool     b_read_special)
{
    AmrLevel::restart(papa, is, b_read_special);

    build_metrics();

    if (level == 0)
    {
        BL_ASSERT(e_field_solver == 0);
        e_field_solver = new Electrostatic(parent, parent->finestLevel(), &phys_bc);
    }
}

void
Accel::build_metrics ()
{
}

void
Accel::set_time_level (Real time,
                     Real dt_old,
                     Real dt_new)
{
    AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

void
Accel::initData ()
{
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Initializing the level data at level " << level << '\n';

    //
    // Initialize this to zero before first solve.
    //
    MultiFab& Phi_new = get_new_data(Elec_Potential_Type);
    Phi_new.setVal(0.);

    if (level == 0)
        init_particles();

//  if (level > 0)
        particle_redistribute();

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Done initializing the level " << level << " data\n";
}

void
Accel::init (AmrLevel& old)
{
    Accel* old_level = (Accel*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new = parent->dtLevel(level);
    Real cur_time = old_level->state[Elec_Field_Type].curTime();
    Real prev_time = old_level->state[Elec_Field_Type].prevTime();

    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    MultiFab& Phi_new = get_new_data(Elec_Potential_Type);
    for (FillPatchIterator fpi(old, Phi_new, 0, cur_time, Elec_Potential_Type, 0, 1);
         fpi.isValid(); ++fpi)
    {
        Phi_new[fpi].copy(fpi());
    }
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
Accel::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = get_level(level-1).state[Elec_Field_Type].curTime();
    Real prev_time = get_level(level-1).state[Elec_Field_Type].prevTime();
    Real dt_old    = (cur_time - prev_time) / (Real)parent->MaxRefRatio(level-1);

    set_time_level(cur_time, dt_old, dt);

    MultiFab& Phi_new = get_new_data(Elec_Potential_Type);
    FillCoarsePatch(Phi_new, 0, cur_time, Elec_Potential_Type, 0, 1);

    // We set dt to be large for this new level to avoid screwing up
    // computeNewDt.
    parent->setDtLevel(1.e100, level);
}

Real
Accel::initial_time_step ()
{
    Real dummy_dt = 0;
    Real init_dt = 0;

    if (initial_dt > 0)
    {
        init_dt = initial_dt;
    }
    else
    {
        init_dt = init_shrink * est_time_step(dummy_dt);
    }

    return init_dt;
}

Real
Accel::est_time_step (Real dt_old)
{
    if (fixed_dt > 0)
        return fixed_dt;

    // This is just a dummy value to start with
    Real est_dt = 1.0e+200;
    particle_est_time_step(est_dt);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Accel::est_time_step at level "
                  << level
                  << ":  estdt = "
                  << est_dt << '\n';

    return est_dt;
}

void
Accel::computeNewDt (int                   finest_level,
                   int                   sub_cycle,
                   Array<int>&           n_cycle,
                   const Array<IntVect>& ref_ratio,
                   Array<Real>&          dt_min,
                   Array<Real>&          dt_level,
                   Real                  stop_time,
                   int                   post_regrid_flag)
{
    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;
    
    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        Accel& adv_level = get_level(i);
        dt_min[i] = adv_level.est_time_step(dt_level[i]);
    }

    if (fixed_dt <= 0.0)
    {
        if (post_regrid_flag == 1)
        {
            //
            // Limit dt's by pre-regrid dt
            //
            for (i = 0; i <= finest_level; i++)
            {
                dt_min[i] = std::min(dt_min[i], dt_level[i]);
            }
            //
            // Find the minimum over all levels
            //
            for (i = 0; i <= finest_level; i++)
            {
                n_factor *= n_cycle[i];
                dt_0 = std::min(dt_0, n_factor * dt_min[i]);
            }
        }
        else
        {
            bool sub_unchanged=true;
            if ((parent->maxLevel() > 0) && (level == 0) &&
                (parent->subcyclingMode() == "Optimal") && 
                (parent->okToRegrid(level) || parent->levelSteps(0) == 0) )
            {
                int new_cycle[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                    new_cycle[i] = n_cycle[i];
                // The max allowable dt
                Real dt_max[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                {
                    dt_max[i] = dt_min[i];
                }
                // find the maximum number of cycles allowed.
                int cycle_max[finest_level+1];
                cycle_max[0] = 1;
                for (i = 1; i <= finest_level; i++)
                {
                    cycle_max[i] = parent->MaxRefRatio(i-1);
                }
                // estimate the amout of work to advance each level.
                Real est_work[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                {
                    est_work[i] = parent->getLevel(i).estimateWork();
                }
                // this value will be used only if the subcycling pattern is changed.
                dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle, dt_max, est_work, cycle_max);
                for (i = 0; i <= finest_level; i++)
                {
                    if (n_cycle[i] != new_cycle[i])
                    {
                        sub_unchanged = false;
                        n_cycle[i] = new_cycle[i];
                    }
                }
                
            }
            
            if (sub_unchanged)
            //
            // Limit dt's by change_max * old dt
            //
            {
                for (i = 0; i <= finest_level; i++)
                {
                    if (verbose && ParallelDescriptor::IOProcessor())
                    {
                        if (dt_min[i] > change_max*dt_level[i])
                        {
                            std::cout << "Accel::compute_new_dt : limiting dt at level "
                                      << i << '\n';
                            std::cout << " ... new dt computed: " << dt_min[i]
                                      << '\n';
                            std::cout << " ... but limiting to: "
                                      << change_max * dt_level[i] << " = " << change_max
                                      << " * " << dt_level[i] << '\n';
                        }
                    }

                    dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
                }
                //
                // Find the minimum over all levels
                //
                for (i = 0; i <= finest_level; i++)
                {
                    n_factor *= n_cycle[i];
                    dt_0 = std::min(dt_0, n_factor * dt_min[i]);
                }
            }
            else
            {
                if (verbose && ParallelDescriptor::IOProcessor())
                {
                   std::cout << "Accel: Changing subcycling pattern. New pattern:\n";
                   for (i = 1; i <= finest_level; i++)
                    std::cout << "   Lev / n_cycle: " << i << " " << n_cycle[i] << '\n';
                }
            }
        }
    }
    else
    {
        dt_0 = fixed_dt;
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
    Real cur_time = state[Elec_Field_Type].curTime();
    if (stop_time >= 0.0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void
Accel::computeInitialDt (int                   finest_level,
                       int                   sub_cycle,
                       Array<int>&           n_cycle,
                       const Array<IntVect>& ref_ratio,
                       Array<Real>&          dt_level,
                       Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    int i;
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    if (parent->subcyclingMode() == "Optimal")
    {
        int new_cycle[finest_level+1];
        for (i = 0; i <= finest_level; i++)
            new_cycle[i] = n_cycle[i];
        Real dt_max[finest_level+1];
        for (i = 0; i <= finest_level; i++)
        {
            dt_max[i] = get_level(i).initial_time_step();
        }
        // Find the maximum number of cycles allowed
        int cycle_max[finest_level+1];
        cycle_max[0] = 1;
        for (i = 1; i <= finest_level; i++)
        {
            cycle_max[i] = parent->MaxRefRatio(i-1);
        }
        // estimate the amout of work to advance each level.
        Real est_work[finest_level+1];
        for (i = 0; i <= finest_level; i++)
        {
            est_work[i] = parent->getLevel(i).estimateWork();
        }
        dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle, dt_max, est_work, cycle_max);
        for (i = 0; i <= finest_level; i++)
        {
            n_cycle[i] = new_cycle[i];
        }
        if (verbose && ParallelDescriptor::IOProcessor() && finest_level > 0)
        {
           std::cout << "Accel: Initial subcycling pattern:\n";
           for (i = 0; i <= finest_level; i++)
               std::cout << "Level " << i << ": " << n_cycle[i] << '\n';
        }
    }
    else
    {
        for (i = 0; i <= finest_level; i++)
        {
            dt_level[i] = get_level(i).initial_time_step();
            n_factor *= n_cycle[i];
            dt_0 = std::min(dt_0, n_factor * dt_level[i]);
        }
    }
    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
    Real cur_time = state[Elec_Field_Type].curTime();
    if (stop_time >= 0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

bool
Accel::writePlotNow ()
{
    if (level > 0)
        BoxLib::Error("Should only call writePlotNow at level 0!");

    return false;
}

void
Accel::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();
    const int ncycle = parent->nCycle(level);
    
    //
    // Remove virtual particles at this level if we have any.
    //
    remove_virtual_particles();

    //
    // Remove Ghost particles on the final iteration 
    //
    if (iteration == ncycle)
        remove_ghost_particles();

    //
    // Redistribute if it is not the last subiteration
    //
    if (iteration < ncycle || level == 0)
    {
         Accel::thePAPC()->Redistribute(false, true, level, field_n_grow);    
    }

    if (level < finest_level)
        average_down();
}

void
Accel::post_restart ()
{
    if (level == 0)
        particle_post_restart(parent->theRestartFile());

    Real cur_time = state[Elec_Field_Type].curTime();

    if (level == 0)
    {

        for (int lev = 0; lev <= parent->finestLevel(); lev++)
        {
            AmrLevel& this_level = get_level(lev);
            e_field_solver->install_level(lev, &this_level);
        }

        // Do multilevel solve here.  We now store phi in the checkpoint file so we can use it
        //  at restart.
        int use_previous_phi_as_guess = 1;
        e_field_solver->multilevel_solve_for_phi(0,parent->finestLevel(),opal_coupling,use_previous_phi_as_guess);
                                                 

        for (int k = 0; k <= parent->finestLevel(); k++)
        {
            const BoxArray& ba = get_level(k).boxArray();
            MultiFab e_field_new(ba, BL_SPACEDIM, 0);
            e_field_solver->get_new_e_field(k, e_field_new, cur_time);
        }
    }
}

void
Accel::postCoarseTimeStep (Real cumtime)
{
}

void
Accel::post_regrid (int lbase,
                  int new_finest)
{
    if (level == lbase)
        particle_redistribute(lbase);


    // Only do solve here if we will be using it in the timestep right after without re-solving
    const Real cur_time = state[Elec_Field_Type].curTime();
    if ((level == lbase) && cur_time > 0)
    {
        int use_previous_phi_as_guess = 1;
        e_field_solver->multilevel_solve_for_phi(level, new_finest, opal_coupling, use_previous_phi_as_guess);
    }
}

void
Accel::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level - 1; k >= 0; k--)
        get_level(k).average_down();

    const Real cur_time = state[Elec_Field_Type].curTime();

    //
    // Solve on full multilevel hierarchy
    //

    e_field_solver->multilevel_solve_for_phi(0, finest_level, opal_coupling);

    // Make this call just to fill the initial state data.
    for (int k = 0; k <= finest_level; k++)
    {
        const BoxArray& ba = get_level(k).boxArray();
        std::cout << ba << std::endl;
        MultiFab e_field_new(ba, BL_SPACEDIM, 0);
        e_field_solver->get_new_e_field(k, e_field_new, cur_time);
    }
}

int
Accel::okToContinue ()
{
    if (level > 0)
        return 1;

    int test = 1;
    if (parent->dtLevel(0) < dt_cutoff)
        test = 0;

    return test;
}

void
Accel::average_down ()
{
    if (level == parent->finestLevel()) return;

    average_down(Elec_Potential_Type);
    average_down(Elec_Field_Type);
}

void
Accel::average_down (int state_index)
{
    if (level == parent->finestLevel()) return;

    Accel& fine_level = get_level(level + 1);

    MultiFab& S_crse = get_new_data(state_index);
    MultiFab& S_fine = fine_level.get_new_data(state_index);
    const int num_comps = S_fine.nComp();

    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_S_fine_BA(S_fine.boxArray().size());

    for (int i = 0; i < S_fine.boxArray().size(); ++i)
    {
        crse_S_fine_BA.set(i, BoxLib::coarsen(S_fine.boxArray()[i],
                                              fine_ratio));
    }

    MultiFab crse_S_fine(crse_S_fine_BA, num_comps, 0);

    {
        for (MFIter mfi(S_fine); mfi.isValid(); ++mfi)
        {
            const int        i        = mfi.index();
            const Box&       overlap  = crse_S_fine_BA[i];
            FArrayBox&       crse_fab = crse_S_fine[i];
            const FArrayBox& fine_fab = S_fine[i];

            BL_FORT_PROC_CALL(FORT_AVGDOWN, fort_avgdown)
                (BL_TO_FORTRAN(crse_fab), num_comps, BL_TO_FORTRAN(fine_fab),
                 overlap.loVect(), overlap.hiVect(), fine_ratio.getVect());
        }
    }

    S_crse.copy(crse_S_fine);
}

void
Accel::alloc_old_data ()
{
    for (int k = 0; k < NUM_STATE_TYPE; k++)
        state[k].allocOldData();
}

void
Accel::remove_old_data ()
{
    AmrLevel::removeOldData();
}

void
Accel::errorEst (TagBoxArray& tags,
               int          clearval,
               int          tagval,
               Real         time,
               int          n_error_buf,
               int          ngrow)
{
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    Array<int> itags;
    Real avg;

    for (int j = 0; j < err_list.size(); j++)
    {
        MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());

        BL_ASSERT(!(mf == 0));

        for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
        {
            int idx = mfi.index();
            const int* lo = mfi.validbox().loVect();
            const int* hi = mfi.validbox().hiVect();

            itags = tags[idx].tags();
            int* tptr = itags.dataPtr();
            const int* tlo = tags[idx].box().loVect();
            const int* thi = tags[idx].box().hiVect();

            Real* dat = (*mf)[mfi].dataPtr();
            const int* dlo = (*mf)[mfi].box().loVect();
            const int* dhi = (*mf)[mfi].box().hiVect();
            const int num_comps = (*mf)[mfi].nComp();

            if (err_list[j].errType() == ErrorRec::Standard)
            {
                RealBox gridloc = RealBox(grids[idx], geom.CellSize(),
                                          geom.ProbLo());
                const Real* xlo = gridloc.lo();
                err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                      &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                      lo, hi, &num_comps, domain_lo, domain_hi,
                                      dx, xlo, prob_lo, &time, &level);
            }
            else if (err_list[j].errType() == ErrorRec::UseAverage)
            {
               err_list[j].errFunc2()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                      &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                      lo, hi, &num_comps, domain_lo, domain_hi,
                                      dx, &level, &avg);
            }

            //
            // Don't forget to set the tags in the TagBox.
            //
            if (allow_untagging == 1)
            {
                tags[idx].tags_and_untags(itags);
            }
            else
            {
                tags[idx].tags(itags);
            }
        }

        delete mf;
    }
}

MultiFab*
Accel::derive (const std::string& name,
             Real               time,
             int                ngrow)
{
    return particle_derive(name, time, ngrow);
}

void
Accel::derive (const std::string& name,
             Real               time,
             MultiFab&          mf,
             int                dcomp)
{
    AmrLevel::derive(name, time, mf, dcomp);
}

void
Accel::GetParticleIDs (Array<int>& part_ids) 
{
    Accel::thePAPC()->GetParticleIDs(part_ids);
}
void
Accel::GetParticleCPU (Array<int>& part_cpu) 
{
    Accel::thePAPC()->GetParticleCPU(part_cpu);
}
void
Accel::GetParticleLocations (Array<Real>& part_locs) 
{
    Accel::thePAPC()->GetParticleLocations(part_locs);
}
void
Accel::GetParticleData (Array<Real>& part_data, int start_comp, int num_comp) 
{
    Accel::thePAPC()->GetParticleData(part_data,start_comp,num_comp);
}
#endif
