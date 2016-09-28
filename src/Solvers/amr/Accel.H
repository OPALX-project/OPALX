
#ifndef _Accel_H_
#define _Accel_H_


#include "Physics/Physics.h"

#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>

#include <ElectroMagneticParticleClass.H>

typedef EMPC<2+3*BL_SPACEDIM> ElectrostaticParticleContainer;

using std::istream;
using std::ostream;

enum StateType {
    Elec_Potential_Type = 0,
    Elec_Field_Type,
    NUM_STATE_TYPE
};

//
// AmrLevel-derived class 
//

class Accel
    :
    public AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    Accel();

    //
    //The basic constructor.
    //
    Accel(Amr& papa, int lev, const Geometry& level_geom, const BoxArray& bl,
        Real time);

    //
    //The destructor.
    //
    virtual ~Accel();

    //
    //Restart from a checkpoint file.
    //
    virtual void restart(Amr& papa, istream& is, bool b_read_special=false);

    //
    //Call AmrLevel::checkPoint and then add radiation info
    //
    virtual void checkPoint(const std::string& dir, std::ostream& os,
                            VisMF::How how, bool dump_old);

    // A string written as the first item in `write_plot_file()` at level zero.
    // It is so we can distinguish between different types of plot files. For
    // Accel it has the form: Accel-Vnnn.
    virtual std::string thePlotFileType() const;

    virtual void setPlotVariables();

    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile(const std::string& dir, ostream& os, VisMF::How how);

    //
    //Write MultiFab as plot file
    //
    void writeMultiFabAsPlotFile(const std::string& pltfile, const MultiFab& mf, std::string componentName);

    //
    //Write all parameters into specified directory.
    //
    static int write_parameters_in_plotfile;
    virtual void write_parameter_file(const std::string& dir);

    //
    // If true then print the warnings from the Fortran routines
    //
    static int print_fortran_warnings;

    //
    //Define data descriptors.
    //
    static void variable_setup();
    static void setup();

    //
    //Define tagging functions.
    //
    static void error_setup();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variable_cleanup();

    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData();

    //
    // Read particle-related inputs
    //
    static void read_particle_params();

    //
    // Write particles in checkpoint directories
    //
    void particle_check_point(const std::string& dir);

    //
    // Write particles in plotfile directories
    //
    void particle_plot_file(const std::string& dir);

    //
    // How to initialize at restart
    //
    void particle_post_restart(const std::string& restart_file);

    //
    // Redistribute
    //
    void particle_redistribute(int lbase = 0, bool init = false);

    //
    // Initialize particle locations and velocities (and strengths if relevant)
    //
    virtual void init_particles();

    //
    // Setup virtual particles if necessary
    //
    void setup_virtual_particles();
    
    //
    // Remove virtual particles if necessary
    //
    void remove_virtual_particles();
    //
    // Setup ghost particles (for finer levels) if necessary
    //
    void setup_ghost_particles();
    
    //
    // Remove ghost particles (for this level) if necessary
    //
    void remove_ghost_particles();

    //
    // Time step control based on particles
    //
    void particle_est_time_step(Real& est_dt);

    //
    // Derived quantities associated with particles
    //
    MultiFab* particle_derive (const std::string& name, Real time, int ngrow);

    //
    // Default verbosity of Particle class
    //
    static int particle_verbose;

    //
    // Default cfl of particles in Particle class
    //
    static Real particle_cfl;
    
    //
    // Should we write particles into plotfiles?
    //
    static int write_particles_in_plotfile;

    //
    // Shall we write the initial single-level particle density into a multifab
    //   called "ParticleDensity"?
    //
    static int write_particle_density_at_init;

    //
    // Shall we write an ascii file with only the particles in every other cell --
    //   this is a cheap way of creating a "coarser" particle file
    //
    static int write_coarsened_particles;

    //
    //Set time levels of state data.
    //
    virtual void set_time_level (Real time, Real dt_old, Real dt_new);

    //
    //Initialize data on this level from another Accel (during regrid).
    //
    virtual void init(AmrLevel& old);

    //
    // Initialize data on this level after regridding if old level did not
    // previously exist
    //
    virtual void init();

    //
    // Proceed with next timestep?
    //
    virtual int okToContinue();

    //
    // Tell Amr to write a plotfile now
    //
    bool writePlotNow();

    //
    // Advance grids at this level in time.
    //
    virtual Real advance(Real time, Real dt, int iteration, int ncycle);

    //
    // Estimate time step.
    //
    Real est_time_step(Real dt_old);

    //
    // Compute initial time step.
    //
    Real initial_time_step();

    //
    // Compute initial `dt'.
    //
    virtual void computeInitialDt(int finest_level, int sub_cycle,
                                  Array<int>& n_cycle,
                                  const Array<IntVect>& ref_ratio,
                                  Array<Real>& dt_level, Real stop_time);
    //
    // Compute new `dt'.
    //
    virtual void computeNewDt(int finest_level, int sub_cycle,
                              Array<int>& n_cycle,
                              const Array<IntVect>& ref_ratio,
                              Array<Real>& dt_min, Array<Real>& dt_level,
                              Real stop_time, int post_regrid_flag);

    //
    // Allocate data at old time.
    //
    virtual void alloc_old_data();

    //
    // Remove data at old time.
    //
    virtual void remove_old_data();

    //
    // Do work after timestep().
    //
    virtual void post_timestep(int iteration);

    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep(Real cumtime);

    //
    // Do work after `regrid()`.
    //
    virtual void post_regrid(int lbase, int new_finest);

    //
    // Do work after a `restart()`.
    //
    virtual void post_restart();

    //
    // Do work after `init()`.
    //
    virtual void post_init(Real stop_time);

    //
    // Error estimation for regridding.
    //
    virtual void errorEst(TagBoxArray& tb, int clearval, int tagval, Real time,
                          int n_error_buf=0, int ngrow=0);

    //
    // Called in grid_places after other tagging routines to modify
    //   the list of tagged points
    //
    virtual void manual_tags_placement (TagBoxArray&    tags,
                                        Array<IntVect>& bf_lev);

    // Returns a MultiFab containing the derived data for this level. The user
    // is responsible for deleting this pointer when done with it. If
    // `ngrow` > 0 the MultiFab is built on the appropriately grown BoxArray.
    MultiFab* derive(const std::string& name, Real time, int ngrow);

    // This version of `derive()` fills the dcomp'th component of mf with the
    // derived quantity.
    void derive(const std::string& name, Real time, MultiFab& mf, int dcomp);

    static int num_grow();

    static ElectrostaticParticleContainer* thePAPC();
    static ElectrostaticParticleContainer* theVirtPC();
    static ElectrostaticParticleContainer* theGhostPC();

    void addOneParticle (int id_in, int cpu_in, 
                         std::vector<double>& xloc, std::vector<double>& attributes);

    void GetParticleIDs       (Array<int>&  part_ids);
    void GetParticleCPU       (Array<int>&  part_cpu);
    void GetParticleLocations (Array<Real>& part_locs);
    void GetParticleData      (Array<Real>& part_data, int start_comp, int num_comp);

protected:

    // 
    // Number of grow cells required for particle motion.
    //
    int field_n_grow;

    static void read_params();

    Accel& get_level(int lev);

    void average_down();
    void average_down(int state_indx);

    void build_metrics();

    //
    // Static data members.
    //
    static bool dump_old;
    static int verbose;
    static int show_timings;
    static Real cfl;
    static Real init_shrink;
    static Real change_max;
    static ErrorList err_list;
    static BCRec phys_bc;
    static int NUM_GROW;

    static int allow_untagging;
    static int do_special_tagging;

    bool FillPatchedOldState_ok;

    // There can be only one Electrostatic object, it covers all levels:
    static class Electrostatic *e_field_solver;
};

//
// Inlines.
//

inline
int
Accel::num_grow()
{
    return NUM_GROW;
}

inline
Accel&
Accel::get_level(int level)
{
    return *(Accel *) &parent->getLevel(level);
}

#endif /*_Accel_H_*/
