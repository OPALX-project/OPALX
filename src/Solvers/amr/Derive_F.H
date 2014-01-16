#ifndef _Derive_F_H_
#define _Derive_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(CA_DERLAPVAR,ca_derlapvar)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERMAGGRAV,ca_dermaggrav)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERNULL,ca_dernull)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERSTATE,ca_derstate)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec,
     const int* level, const int* grid_no);
#endif
