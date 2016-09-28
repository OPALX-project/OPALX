#include "Accel.H"
#include "Accel_F.H"
#include "Electrostatic.H"
#include "Electrostatic_F.H"
#include <Particles_F.H>
 
using std::string;
 
Real
Accel::advance (Real time,
                Real dt,
                int  iteration,
                int  ncycle)
{
    // Need ncycle + 1 to catch all relevant ghost particles and 1 more to move them.
    field_n_grow = ncycle + 2;

    // This must match with the component we specify in TrackRun.cpp for Evec.
    // 0: Q, (1,2,3) : (P_x, P_y, P_z), (4,5,6) : (E_x, E_y, E_z)
    int zmom_comp = 3;
    int    v_comp = 1;
    int    e_comp = 4;
    int    b_comp = 7;
 
    const int finest_level = parent->finestLevel();
    int finest_level_to_advance;
    bool nosub = !parent->subCycle();
    
    if (nosub)
    {
        if (level > 0)
            return dt;
            
        finest_level_to_advance = finest_level;
    }
    else
    {
        finest_level_to_advance = level;

        // We must setup virtual and Ghost Particles
        //
        // Setup the virtual particles that represent finer level particles
        // 
        setup_virtual_particles();

        //
        // Setup ghost particles for use in finer levels. Note that Ghost particles
        // that will be used by this level have already been created, the
        // particles being set here are only used by finer levels.
        //
        for(int lev = level; lev <= finest_level_to_advance && lev < finest_level; lev++)
        {
           get_level(lev).setup_ghost_particles();
        }
    }

    Real dt_lev;
    const Real strt = ParallelDescriptor::second();

    //
    // Move current data to previous, clear current.
    // Don't do this if a coarser level has done this already.
    //
    if (level == 0 || iteration > 1)
    {
        for (int lev = level; lev <= finest_level; lev++)
        {
            dt_lev = parent->dtLevel(lev);
            for (int k = 0; k < NUM_STATE_TYPE; k++)
            {
                get_level(lev).state[k].allocOldData();
                get_level(lev).state[k].swapTimeLevels(dt_lev);
            }
        }
    }

    const Real cur_time  = state[Elec_Field_Type].curTime();

    //
    // We now do a multilevel solve for old Electrostatic. This goes to the 
    // finest level regardless of subcycling behavior. 
    // Thus if we are subcycling we skip this step on the 
    // first iteration of finer levels.
    if (level == 0 || iteration > 1)
    {
        // fix fluxes on finer grids
        {
            for (int lev = level; lev < finest_level; lev++)
            {
                e_field_solver->zero_phi_flux_reg(lev + 1);
            }
        }
        
        // swap grav data
        for (int lev = level; lev <= finest_level; lev++)
            get_level(lev).e_field_solver->swap_time_levels(lev);

        if (show_timings)
        {
            const int IOProc = ParallelDescriptor::IOProcessorNumber();
            Real end = ParallelDescriptor::second() - strt;
            ParallelDescriptor::ReduceRealMax(end,IOProc);
            if (ParallelDescriptor::IOProcessor())
               std::cout << "Time before solve for old phi " << end << '\n';
        }

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "\n... old-time level solve at level " << level << '\n';
            
        //
        // Solve for phi
        // If a single-level calculation we can still use the previous phi as a guess.
        int use_previous_phi_as_guess = 1;
        e_field_solver->multilevel_solve_for_old_phi(level, finest_level, opal_coupling, 
                                                     use_previous_phi_as_guess);
    }

    //
    // Advance Particles
    //
    // Advance the particle velocities to the half-time and the positions to the new time
    // We use the cell-centered e_field to correctly interpolate onto particle locations
    ParallelDescriptor::Barrier();
    if (particle_verbose && ParallelDescriptor::IOProcessor())
        std::cout << "MKD ... updating particle positions and velocity\n";
    
    Real gammaz = Accel::thePAPC()->sumParticleMomenta(zmom_comp);
    Real betaC  = gammaz / ( sqrt(gammaz*gammaz+1) * Physics::c );

    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
        // We need field_n_grow grow cells to track boundary particles
        const BoxArray& ba = get_level(lev).get_new_data(Elec_Field_Type).boxArray();
        MultiFab e_field_old(ba, BL_SPACEDIM, field_n_grow);
        get_level(lev).e_field_solver->get_old_e_field(lev, e_field_old, time);
        
        Accel::thePAPC()->FillEvec(e_field_old, lev, betaC, e_comp, b_comp);
        Accel::thePAPC()->MKD     (             lev, dt, field_n_grow, v_comp, e_comp, b_comp);

        // Only need the coarsest virtual particles here.
        if (lev == level && level < finest_level)
        {
            if (Accel::theVirtPC() != 0)
            {
                Accel::theVirtPC()->FillEvec(e_field_old, lev, betaC, e_comp, b_comp);
                Accel::theVirtPC()->MKD     (             lev, dt, field_n_grow, v_comp, e_comp, b_comp);
            }
        }

        // Miiiight need all Ghosts
        if (Accel::theGhostPC() != 0)
        {
            Accel::theGhostPC()->FillEvec(e_field_old, lev, betaC, e_comp, b_comp);
            Accel::theGhostPC()->MKD     (             lev, dt, field_n_grow, v_comp, e_comp, b_comp);
        }
    }

    ParallelDescriptor::Barrier();
    if (particle_verbose && ParallelDescriptor::IOProcessor())
        std::cout << "MKD is complete ...\n";

    //
    // Here we use the "old" phi from the current time step as a guess for this
    // solve
    //
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
        MultiFab::Copy(parent->getLevel(lev).get_new_data(Elec_Potential_Type),
                       parent->getLevel(lev).get_old_data(Elec_Potential_Type),
                       0, 0, 1, 0);
    }

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Time before solve for new phi " << end << '\n';
    }

    // Solve for new Electrostatic
    int use_previous_phi_as_guess = 1;
    if (finest_level_to_advance > level)
    {
        e_field_solver->multilevel_solve_for_new_phi(level, finest_level_to_advance, opal_coupling, 
                                                     use_previous_phi_as_guess);
    }
    else
    {
        int fill_interior = 0;
        e_field_solver->solve_for_new_phi(level,get_new_data(Elec_Potential_Type),
                                          e_field_solver->get_grad_phi_curr(level),
                                          opal_coupling, fill_interior, field_n_grow);
    }

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Time  after solve for new phi " << end << '\n';
    }

    // Advance the particle velocities by dt/2 to the new time. We use the
    // cell-centered E field to correctly interpolate onto particle
    // locations.
    ParallelDescriptor::Barrier();
    if (particle_verbose && ParallelDescriptor::IOProcessor())
        std::cout << "MK ... updating velocity only\n";

    // Compute this again since the momenta changed in the MKD call
    gammaz = Accel::thePAPC()->sumParticleMomenta(zmom_comp);
    betaC  = gammaz / ( sqrt(gammaz*gammaz+1) * Physics::c );

    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
        const BoxArray& ba = get_level(lev).get_new_data(Elec_Field_Type).boxArray();
        MultiFab e_field_new(ba, BL_SPACEDIM, field_n_grow);
        get_level(lev).e_field_solver->get_new_e_field(lev, e_field_new, cur_time);

        Accel::thePAPC()->FillEvec(e_field_new, lev, betaC, e_comp, b_comp);
        Accel::thePAPC()->MK      (             lev, dt   , v_comp, e_comp, b_comp);

        // Virtual particles will be recreated, so we need not kick them.

        // Ghost particles need to be kicked except during the final iteration.
        if (iteration != ncycle)
        {
            if (Accel::theGhostPC() != 0)
            {
                 Accel::theGhostPC()->FillEvec(e_field_new, lev, betaC, e_comp, b_comp);
                 Accel::theGhostPC()->MK      (             lev, dt   , v_comp, e_comp, b_comp);
            }
        }
    } 
    ParallelDescriptor::Barrier();
    if (particle_verbose && ParallelDescriptor::IOProcessor())
        std::cout << "MK is complete ...\n";

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Time  at end of routine " << end << '\n';
    }

    // Redistribution happens in post_timestep
    return dt;
}
