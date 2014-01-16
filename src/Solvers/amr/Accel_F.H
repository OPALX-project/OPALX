#ifndef _Accel_F_H_
#define _Accel_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(PARALLEL_INIT_FORTRAN, parallel_init_fortran)();

BL_FORT_PROC_DECL(FILCC, filcc)
    (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
     const int * domlo, const int * domhi,
     const Real * dx_crse, const Real * xlo,
     const int * bc);

BL_FORT_PROC_DECL(GENERIC_FILL, generic_fill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(SUM_OVER_LEVEL, sum_over_level)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real dx[], Real* sum);

BL_FORT_PROC_DECL(FORT_AVGDOWN, fort_avgdown)
    (BL_FORT_FAB_ARG(crse_fab), const int& nc,
     const BL_FORT_FAB_ARG(fine_fab),
     const int ovlo[], const int ovhi[], const int rat[]);
#endif

