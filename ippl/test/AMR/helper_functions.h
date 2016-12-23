#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>

#ifdef UNIQUE_PTR
typedef Array<std::unique_ptr<MultiFab> > container_t;
#else
#include <PArray.H>
typedef PArray<MultiFab> container_t;
#endif


void writePotential(const container_t& phi, double dx, double dlow, std::string filename) {
    // potential and efield a long x-direction (y = 0, z = 0)
    int lidx = 0;
#ifdef UNIQUE_PTR
    for (MFIter mfi(*phi[0 /*level*/]); mfi.isValid(); ++mfi) {
#else
    for (MFIter mfi(phi[0 /*level*/]); mfi.isValid(); ++mfi) {
#endif
        const Box& bx = mfi.validbox();
#ifdef UNIQUE_PTR
        const FArrayBox& lhs = (*phi[0])[mfi];
#else
        const FArrayBox& lhs = (phi[0])[mfi];
#endif
        /* Write the potential of all levels to file in the format: x, y, z, \phi
         */
        for (int proc = 0; proc < ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == ParallelDescriptor::MyProc() ) {
                std::string outfile = filename + std::to_string(0);
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                    out.open(outfile, std::ios_base::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
//                     for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
//                         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                            IntVect ivec(i, j, k);
                            // add one in order to have same convention as PartBunch::computeSelfField()
                            out << i + 1 << " " << j + 1 << " " << k + 1 << " " << dlow + dx * i << " "
                                << lhs(ivec, 0)  << std::endl;
//                         }
//                     }
                }
                out.close();
            }
            ParallelDescriptor::Barrier();
        }
    }
}

void writeElectricField(const container_t& efield, double dx, double dlow) {
    // potential and efield a long x-direction (y = 0, z = 0)
    int lidx = 0;
#ifdef UNIQUE_PTR
    for (MFIter mfi(*efield[0 /*level*/]); mfi.isValid(); ++mfi) {
#else
    for (MFIter mfi(efield[0 /*level*/]); mfi.isValid(); ++mfi) {
#endif
        const Box& bx = mfi.validbox();
#ifdef UNIQUE_PTR
        const FArrayBox& lhs = (*efield[0])[mfi];
#else
        const FArrayBox& lhs = (efield[0])[mfi];
#endif
        /* Write the potential of all levels to file in the format: x, y, z, ex, ey, ez
         */
        for (int proc = 0; proc < ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == ParallelDescriptor::MyProc() ) {
                std::string outfile = "amr-efield_scalar-level-" + std::to_string(0);
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                    out.open(outfile, std::ios_base::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
//                     for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
//                         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                            IntVect ivec(i, j, k);
                            // add one in order to have same convention as PartBunch::computeSelfField()
                            out << i + 1 << " " << j + 1 << " " << k + 1 << " " << dlow + dx * i << " "
                                << lhs(ivec, 0) << " " << lhs(ivec, 1) << " " << lhs(ivec, 2) << std::endl;
//                         }
//                     }
                }
                out.close();
            }
            ParallelDescriptor::Barrier();
        }
    }
}


double totalFieldEnergy(container_t& efield, const Array<int>& rr) {
    
    for (int lev = efield.size() - 2; lev >= 0; lev--)
        BoxLib::average_down(efield[lev+1], efield[lev], 0, 3, rr[lev]);
    
    double fieldEnergy = 0.0;
//     long volume = 0;
//     for (int l = 0; l < efield.size(); ++l) {
        long nPoints = 0;

        int l = 0;
#ifdef UNIQUE_PTR
        fieldEnergy += MultiFab::Dot(*efield[l], 0, *efield[l], 0, 3, 0);
        for (unsigned int b = 0; b < efield[l]->boxArray().size(); ++b) {
//             volume += efield[l]->box(b).volume(); // cell-centered --> volume() == numPts
            nPoints += efield[l]->box(b).numPts();
#else
        fieldEnergy += MultiFab::Dot(efield[l], 0, efield[l], 0, 3, 0);
        for (unsigned int b = 0; b < efield[l].boxArray().size(); ++b) {
            nPoints += efield[l].box(b).numPts();
#endif
        }
        fieldEnergy /= double(nPoints);
//     }
    
    return 0.5 * /*Physics::epsilon_0 * */fieldEnergy;
}

#endif