#ifndef _Electrostatic_F_H_
#define _Electrostatic_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(FORT_AVGDOWN_PHI, fort_avgdown_phi)
    (BL_FORT_FAB_ARG(crse_fab),
     const BL_FORT_FAB_ARG(fine_fab),
     const int ovlo[], const int ovhi[], const int rat[]);

BL_FORT_PROC_DECL(FORT_EDGE_INTERP, fort_edge_interp)
    (const int* flo, const int* fhi,
     const int* nc, const int* refRatio, const int* dir,
     BL_FORT_FAB_ARG(fine));

BL_FORT_PROC_DECL(FORT_PC_EDGE_INTERP, fort_pc_edge_interp)
    (const int* lo, const int* hi,
     const int* nc, const int* refRatio, const int* dir,
     const BL_FORT_FAB_ARG(crse),
     BL_FORT_FAB_ARG(fine));

BL_FORT_PROC_DECL(CA_AVG_EC_TO_CC, fort_avg_ec_to_cc)
    (const int* lo, const int* hi,
     const int* bc_lo, const int* symmetry_type,
     BL_FORT_FAB_ARG(phi_cc),
     D_DECL(const BL_FORT_FAB_ARG(xgrad),
            const BL_FORT_FAB_ARG(ygrad),
            const BL_FORT_FAB_ARG(zgrad)),
     const Real* problo);

BL_FORT_PROC_DECL(CA_COMPUTE_EC, fort_compute_ec)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(phi_cc),
     D_DECL(const BL_FORT_FAB_ARG(xgrad),
            const BL_FORT_FAB_ARG(ygrad),
            const BL_FORT_FAB_ARG(zgrad)),
     const Real* dx);

BL_FORT_PROC_DECL(CA_AVERAGE_EC, fort_average_ec)
    (D_DECL(const BL_FORT_FAB_ARG(xgrad_fine),
            const BL_FORT_FAB_ARG(ygrad_fine),
            const BL_FORT_FAB_ARG(zgrad_fine)),
     D_DECL(const BL_FORT_FAB_ARG(xgrad_crse),
            const BL_FORT_FAB_ARG(ygrad_crse),
            const BL_FORT_FAB_ARG(zgrad_crse)),
     const int ovlo[], const int ovhi[], const int rat[]);

BL_FORT_PROC_DECL(CA_TEST_RESIDUAL, ca_test_residual)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(rhs),
     D_DECL(const BL_FORT_FAB_ARG(xgrad),
            const BL_FORT_FAB_ARG(ygrad),
            const BL_FORT_FAB_ARG(zgrad)),
     const Real* dx, const Real* problo,
     const int* coord_type);

BL_FORT_PROC_DECL(FORT_SET_HOMOG_BCS, fort_set_homog_bcs)
    (const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     BL_FORT_FAB_ARG(phi),
     const Real* dx);

BL_FORT_PROC_DECL(FORT_ADD_MONOPOLE_BCS, fort_add_monopole_bcs)
    (const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const int* nparticles,
     const Real* part_data,
     BL_FORT_FAB_ARG(phi), const Real* dx);
#endif

