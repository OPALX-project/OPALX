#ifndef _Electrostatic_H_
#define _Electrostatic_H_

#include <AmrLevel.H>
#include <FluxRegister.H>
#include <Particles.H>

#include "TrilinosSolver.h"

class Electrostatic {

public:

    Electrostatic(Amr* Parent, int _finest_level, BCRec* _phys_bc);
    virtual ~Electrostatic();

    void read_params();

    void install_level(int level, AmrLevel* level_data_to_install);

    int get_no_sync();
    int get_no_composite();

    PArray<MultiFab>& get_grad_phi_prev(int level);
    PArray<MultiFab>& get_grad_phi_curr(int level);

    void plus_phi_curr(int level, MultiFab& addend);
    void plus_grad_phi_curr(int level, PArray<MultiFab>& addend);

    // 
    // Moves curr phi data to prev and resets curr to 1.e50
    //
    void swap_time_levels(int level);

    void solve_for_old_phi(int level, MultiFab& phi, PArray<MultiFab>& grad_phi,
                           Real opal_coupling, int fill_interior);
    void solve_for_new_phi(int level, MultiFab& phi, PArray<MultiFab>& grad_phi,
                           Real opal_coupling, int fill_interior, int field_n_grow = 1);
    void solve_for_phi(int level, MultiFab& Rhs, MultiFab& phi,
                       PArray<MultiFab>& grad_phi, 
                       Real opal_coupling, Real time,
                       int fill_interior);
    void solve_for_delta_phi(int crse_level, int fine_level, MultiFab& CrseRhs,
                             PArray<MultiFab>& delta_phi,
                             PArray<PArray<MultiFab> >& grad_delta_phi);

    void e_field_sync(int crse_level, int fine_level,
                      const MultiFab& drho_and_drhoU, const MultiFab& dphi,
                      PArray<MultiFab>& grad_delta_phi_cc);

    void multilevel_solve_for_old_phi(int level, int finest_level, Real opal_coupling,
                                      int use_previous_phi_as_guess=0);
    void multilevel_solve_for_new_phi(int level, int finest_level,Real opal_coupling,
                                      int use_previous_phi_as_guess=0);
    void multilevel_solve_for_phi(int level, int finest_level, Real opal_coupling,
                                  int use_previous_phi_as_guess=0);
    void actual_multilevel_solve(int level, int finest_level, Real opal_coupling,
                                 Array<PArray<MultiFab> >& grad_phi, int is_new,
                                 int use_previous_phi_as_guess=0);

    Real compute_level_average(int level, MultiFab* mf);
    Real compute_multilevel_average(int level, MultiFab* mf, int flev = -1);

    void get_crse_grad_phi(int level, PArray<MultiFab>& grad_phi_crse,
                           Real time);
    void get_crse_phi(int level, MultiFab& phi_crse, Real time);

    //
    // Sets phi_flux registers to 0
    //
    void zero_phi_flux_reg(int level);

    void get_old_e_field(int level, MultiFab& e_field, Real time);
    void get_new_e_field(int level, MultiFab& e_field, Real time);

    static void average_down(MultiFab& crse, const MultiFab& fine,
                             const IntVect& ratio);

    void average_fine_ec_onto_crse_ec(int level, int is_new);

    void compute_face_fluxes(BoxArray& these_grids, PArray<MultiFab>& grad_phi, MultiFab& phi, const Real* dx);

    void add_to_fluxes(int level, int iteration, int ncycle);

    void reflux_phi(int level, MultiFab& dphi);

//  void make_mg_bc();

    void set_dirichlet_bcs(int level, MultiFab* phi);

protected:
    //
    // Pointers to amr,amrlevel.
    //
    Amr*             parent;
    PArray<AmrLevel> level_data;
    //
    // Pointers to grad_phi at previous and current time
    //
    Array< PArray<MultiFab> > grad_phi_curr;
    Array< PArray<MultiFab> > grad_phi_prev;

    PArray<FluxRegister> phi_flux_reg;
    //
    // BoxArray at each level
    //
    Array<BoxArray> grids;
    //
    // Resnorm at each level
    //
    Array<Real> level_solver_resnorm;

    int density;
    int finest_level;
    int finest_level_allocated;

    int mg_bc[2*BL_SPACEDIM];

    BCRec* phys_bc;

    static int verbose;
    static int show_timings;
    static int no_sync;
    static int no_composite;

    void fill_ec_grow(int level, PArray<MultiFab>& ecF,
                      const PArray<MultiFab>& ecC) const;

    void AddParticlesToRhs(int level, MultiFab& Rhs, int ngrow);
    void AddGhostParticlesToRhs(int level, MultiFab& Rhs, int ngrow);
    void AddVirtualParticlesToRhs(int level, MultiFab& Rhs, int ngrow);

    void AddParticlesToRhs(int base_level, int finest_level, PArray<MultiFab>& Rhs_particles);
    void AddGhostParticlesToRhs(int level, PArray<MultiFab>& Rhs_particles);
    void AddVirtualParticlesToRhs(int finest_level, PArray<MultiFab>& Rhs_particles);

    Solver* init_trilinos(PArray<MultiFab>& rhs, PArray<MultiFab>& soln, const Real* dx);

};
#endif

