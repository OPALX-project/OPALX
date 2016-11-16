#include <iomanip>
#include <Accel.H>

#include <Electrostatic.H>
#include <Accel_F.H>
#include <Particles_F.H>

namespace
{
    bool virtual_particles_set = false;
    
    std::string ascii_particle_file;
    std::string binary_particle_file;

    const std::string chk_particle_file("PA");

    //
    // Containers for the real "active" Particles
    //
    ElectrostaticParticleContainer* PAPC = 0;
    //
    // Container for temporary, virtual Particles
    //
    ElectrostaticParticleContainer* VirtPC  = 0;
    //
    // Container for temporary, ghost Particles
    //
    ElectrostaticParticleContainer* GhostPC  = 0;

    //
    // We want to call this routine when on exit to clean up particles.
    //
    void RemoveParticlesOnExit ()
    {
        delete PAPC;
        PAPC = 0;

        delete VirtPC;
        VirtPC = 0;

        delete GhostPC;
        GhostPC = 0;
    }
}

int Accel::particle_verbose               = 1;
int Accel::write_particles_in_plotfile    = 1;
int Accel::write_particle_density_at_init = 0;
int Accel::write_coarsened_particles      = 0;
Real Accel::particle_cfl = 0.5;

ElectrostaticParticleContainer*
Accel::thePAPC ()
{
    return PAPC;
}

ElectrostaticParticleContainer*
Accel::theVirtPC ()
{
    return VirtPC;
}
ElectrostaticParticleContainer*
Accel::theGhostPC ()
{
    return GhostPC;
}

void
Accel::read_particle_params ()
{
    ParmParse pp("nyx");

    //
    // Control the verbosity of the Particle class
    //
    ParmParse ppp("particles");
    ppp.query("v", particle_verbose);
    ppp.query("write_in_plotfile", write_particles_in_plotfile);
    //
    // Set the cfl for particle motion (fraction of cell that a particle can
    // move in a timestep).
    //
    ppp.query("cfl", particle_cfl);
}

void
Accel::init_particles ()
{
    if (level > 0)
        return;

    //
    // Need to initialize particles before defining initial electric field.
    //
    BL_ASSERT (PAPC == 0);
    PAPC = new ElectrostaticParticleContainer(parent);

    PAPC->SetAllowParticlesNearBoundary(true);

    if (parent->subCycle())
    {
        VirtPC = new ElectrostaticParticleContainer(parent);
        GhostPC = new ElectrostaticParticleContainer(parent);
    }

    //
    // 2 gives more stuff than 1.
    //
    PAPC->SetVerbose(particle_verbose);
}

void
Accel::addOneParticle (int id_in, int cpu_in, 
                       std::vector<double>& xloc, std::vector<double>& attributes)

{
    if (level > 0)
        return;

    PAPC->addOneParticle(id_in,cpu_in,xloc,attributes);
}

void
Accel::particle_post_restart (const std::string& restart_file)
{
    if (level > 0)
        return;
     
    BL_ASSERT(PAPC == 0);
    PAPC = new ElectrostaticParticleContainer(parent);

    if (parent->subCycle())
    {
        VirtPC = new ElectrostaticParticleContainer(parent);
        GhostPC = new ElectrostaticParticleContainer(parent);
    }

    //
    // Make sure to call RemoveParticlesOnExit() on exit.
    //
    BoxLib::ExecOnFinalize(RemoveParticlesOnExit);
    //
    // 2 gives more stuff than 1.
    //
    PAPC->SetVerbose(particle_verbose);
    PAPC->Restart(restart_file, chk_particle_file);
    //
    // We want the ability to write the particles out to an ascii file.
    //
    ParmParse pp("particles");

    std::string particle_output_file;
    pp.query("particle_output_file", particle_output_file);

    if (!particle_output_file.empty())
    {
        PAPC->WriteAsciiFile(particle_output_file);
    }
}

void
Accel::particle_est_time_step (Real& est_dt)
{
    MultiFab& e_field = get_new_data(Elec_Field_Type);

    long total_num = PAPC->NumberOfParticlesAtLevel(0,true,false);
    long valid_num = PAPC->NumberOfParticlesAtLevel(0,false,false);
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "NUM OF GLOBL VALID PARTS AT LEVEL IN EST_TIME_STEP " << valid_num << std::endl;
        std::cout << "NUM OF GLOBL TOTAL PARTS AT LEVEL IN EST_TIME_STEP " << total_num << std::endl;
    }

    //FIXME BoxLib was updated --> no ParticleContainer<NR,NI>::estTimestep function anymore
    const Real est_dt_particle = 1.0e-5; //PAPC->estTimestep(e_field, level, particle_cfl);

    if (est_dt_particle > 0)
        est_dt = std::min(est_dt, est_dt_particle);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        if (est_dt_particle > 0)
        {
            std::cout << "...estdt from particles at level "
                      << level << ": " << est_dt_particle << '\n';
        }
        else
        {
            std::cout << "...there are no particles at level "
                      << level << '\n';
        }
    }
}

void
Accel::particle_redistribute (int lbase, bool init)
{
    if (PAPC)
    {
        //
        // If we are calling with init = true, then we want to force the redistribute
        //    without checking whether the grids have changed.
        //
        if (init)
        {
            PAPC->Redistribute(false, false, lbase);
            return;
        }
        //
        // These are usually the BoxArray and DistributionMap from the last regridding.
        //
        static Array<BoxArray>            ba;
        static Array<DistributionMapping> dm;

        bool changed = false;

        int flev = parent->finestLevel();
	
        while (parent->getAmrLevels()[flev] == nullptr)
            flev--;
 
        if (int(ba.size()) != flev+1)
        {
            ba.resize(flev+1);
            dm.resize(flev+1);
            changed = true;
        }
        else
        {
            for (int i = 0; i <= flev && !changed; i++)
            {
                if (ba[i] != parent->boxArray(i))
                    //
                    // The BoxArrays have changed in the regridding.
                    //
                    changed = true;

                if (!changed)
                {
                    if (dm[i] != parent->getLevel(i).get_new_data(0).DistributionMap())
                        //
                        // The DistributionMaps have changed in the regridding.
                        //
                        changed = true;
                }
            }
        }

        if (changed)
        {
            //
            // We only need to call Redistribute if the BoxArrays or DistMaps have changed.
	    // We also only call it for particles >= lbase. This is
	    // because of we called redistribute during a subcycle, there may be particles not in
	    // the proper position on coarser levels.
            //
            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "Calling redistribute because changed " << '\n';

            PAPC->Redistribute(false, false, lbase);
            //
            // Use the new BoxArray and DistMap to define ba and dm for next time.
            //
            for (int i = 0; i <= flev; i++)
            {
                ba[i] = parent->boxArray(i);
                dm[i] = parent->getLevel(i).get_new_data(0).DistributionMap();
            }
        }
        else
        {
            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "NOT calling redistribute because NOT changed " << '\n';
        }
    }
}

void
Accel::setup_virtual_particles()
{
    if(Accel::thePAPC() != 0 && !virtual_particles_set)
    {
        ElectrostaticParticleContainer::PBox virts;
        if (level < parent->finestLevel())
        {
    	    get_level(level + 1).setup_virtual_particles();
	    Accel::theVirtPC()->CreateVirtualParticles(level+1, virts);
	    Accel::theVirtPC()->AddParticlesAtLevel(level, virts);
	    Accel::thePAPC()->CreateVirtualParticles(level+1, virts);
	    Accel::theVirtPC()->AddParticlesAtLevel(level, virts);
        }
        virtual_particles_set = true;
    }
}

void
Accel::remove_virtual_particles()
{
    if (Accel::theVirtPC() != 0)
        Accel::theVirtPC()->RemoveParticlesAtLevel(level);
    virtual_particles_set = false;
}

void
Accel::setup_ghost_particles()
{
    BL_ASSERT(level < parent->finestLevel());
    ElectrostaticParticleContainer::PBox ghosts;
    Accel::thePAPC()->CreateGhostParticles(level,Accel::field_n_grow-1, ghosts);
    Accel::theGhostPC()->AddParticlesAtLevel(level+1, ghosts, true);
}

void
Accel::remove_ghost_particles()
{
    if (Accel::theGhostPC() != 0)
        Accel::theGhostPC()->RemoveParticlesAtLevel(level);
}
