/* Cosmic initial conditions, roughly following MC2 initialization routines.

   Final output (particle positions and velocities) are in grid units. 

                        Zarija Lukic, May 2009
                           zarija@lanl.gov
*/

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TypesAndDefs.h"
#include "Parallelization.h"
#include "DataBase.h"
#include "MT_Random.h"
#include "Cosmology.h"
#include "PerfMon.h"
#ifndef IPPLCREATE
#include "distribution.h"
#endif
#include "Basedata.h"
#include "Initializer.h"
#ifndef IPPLCREATE
#include "distribution.h"
#endif


#ifdef FFTW2
#include "fftw_mpi.h"
#endif
#ifdef FFTW3
#include <fftw3-mpi.h>
#endif

#define pi M_PI
#define Ndim 3  /* Dimension of the problem */

#ifdef INITIALIZER_NAMESPACE
namespace initializer {
#endif

ParallelTools Parallel;
GlobalStuff DataBase;
CosmoClass Cosmology;

real *Pk;
my_fftw_complex *rho;
long My_Ng;  // lenght of the above arrays = number of my grid points

integer nx1, nx2, ny1, ny2, nz1, nz2; // end points of my domain
integer ngx, ngy, ngz;                // # of zones in my domain

integer MyPE, MasterPE, NumPEs;   // MPI stuff

#ifdef FFTW2
fftwnd_mpi_plan plan;
int lnx, lxs, lnyt, lyst, lsize;
#endif
#ifdef FFTW3
fftw_plan plan;
#endif


void CreateMPI_FFTW_COMPLEX(MPI_Datatype *MPI_FFTW_COMPLEX){
   MPI_Datatype typelist[3] = {MPI_FLOAT, MPI_FLOAT, MPI_UB};
   MPI_Aint displ[3];
   int blocklen[3] = {1, 1, 1};
   int i, base;
      
   /* compute displacements of structure components */	
   MPI_Address(&rho[0].re, displ);
   MPI_Address(&rho[0].im, displ+1);
   MPI_Address(rho+1, displ+2);
   base = displ[0];
   for (i=0; i<3; ++i) displ[i] -= base;
      
   /* build datatype describing structure */
   MPI_Type_struct(3, blocklen, displ, typelist, MPI_FFTW_COMPLEX);
      //MPI_Type_contiguous(2, MPI_FLOAT, MPI_FFTW_COMPLEX);
   MPI_Type_commit(MPI_FFTW_COMPLEX);
   return;
}


void initspec(){
   const real tpi=2.0*pi;
   const real n_s = DataBase.n_s;
   const real sigma8=DataBase.Sigma_8;
   const integer ngrid=DataBase.ngrid;
   const integer nq=ngrid/2;
   const real tpiL=tpi/DataBase.box_size; // k0, physical units
   long i, j, k, index, k_i, k_j, k_k;
   real kk, trans_f, s8, norm;
      
   /* Set P(k)=T^2*k^n array.
   Linear array, taking care of the mirror symmetry
   (reality condition for the density field). */
   for (i=0; i<ngx; ++i){
      k_i = i+nx1;
      if (k_i >= nq) {k_i = -MOD(ngrid-k_i,ngrid);}
      for (j=0; j<ngy; ++j){
         k_j = j+ny1;
         if (k_j >= nq) {k_j = -MOD(ngrid-k_j,ngrid);}
         for (k=0; k<ngz; ++k){
            k_k = k+nz1;
            if (k_k >= nq) {k_k = -MOD(ngrid-k_k,ngrid);}
            index = (i*ngy+j)*ngz+k;
            kk = tpiL*sqrt(k_i*k_i+k_j*k_j+k_k*k_k);
            trans_f = Cosmology.TransferFunction(kk);
            Pk[index] = trans_f*trans_f*pow(kk, n_s);
         }
      }
   }
      
   // Sigma_8 normalization:
   s8 = Cosmology.Sigma_r(8.0, 1.0);
   norm = sigma8*sigma8/(s8*s8);
   s8 = Cosmology.Sigma_r(8.0, norm);
   if (MyPE == MasterPE)
      printf("\nsigma_8 = %f, target was %f\n", s8, sigma8);
      
   // Finally, P(k)=A*T^2*k^n:
   for (i=0; i<My_Ng; ++i) Pk[i] *= norm;

   return;
}


void indens1D(){
   const real tpi=2.0*pi;
   const integer ngrid=DataBase.ngrid;
   const integer nq=ngrid/2;
   long i, j, k, index, index_cc;
   int proc;
   integer is_conj, rn_size;
   real amp, rn;
   real *am_plane, *ph_plane;

   is_conj=Parallel.GetMirrorFlag();
      
   if (is_conj == 0){
      proc = Parallel.GetMyPE();
      init_random(DataBase.seed, proc);
      for (i=0; i<My_Ng; ++i){
         rn  = genrand_real();  // Amplitude
         amp = sqrt(-1.0*log(rn)*Pk[i]);
         rn  = genrand_real();  // Phase
         rho[i].re = amp*cos(tpi*rn);
         rho[i].im = amp*sin(tpi*rn);
      }
   }
   else if (is_conj == 1) {
      rn_size = ngy*ngz;
      am_plane = (real *)malloc(rn_size*sizeof(real));
      ph_plane = (real *)malloc(rn_size*sizeof(real));

      // First plane (i=0):
      proc = Parallel.GetMirrorPE()+1;
      init_random(DataBase.seed, proc);
      for (i=0; i<rn_size; ++i){
         am_plane[i] = genrand_real();
         ph_plane[i] = genrand_real();
      }
      for (j=0; j<ngy; ++j){
         for (k=0; k<ngz; ++k){
            index = j*ngz+k;
            index_cc = MOD(ngrid-j,ngrid)*ngz + MOD(ngrid-k,ngrid);
            amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
            rn = ph_plane[index_cc];
            rho[index].re = amp*cos(tpi*rn);
            rho[index].im = -1.0*amp*sin(tpi*rn);
         }
      }
               
      // Other planes:
      proc = Parallel.GetMirrorPE();
      init_random(DataBase.seed, proc);
      for (i=0; i<rn_size; ++i){ // Skip the first plane
         am_plane[i] = genrand_real();
         ph_plane[i] = genrand_real();
      }
      for (i=ngx-1; i>0; --i){ // Go through other planes
         for (k=0; k<rn_size; ++k){
            am_plane[k] = genrand_real();
            ph_plane[k] = genrand_real();
         }
         for (j=0; j<ngy; ++j){
            for (k=0; k<ngz; ++k){
               index = (i*ngy+j)*ngz+k;
               index_cc = MOD(ngrid-j,ngrid)*ngz + MOD(ngrid-k,ngrid);
               amp = sqrt(-1.0*log(am_plane[index_cc])*Pk[index]);
               rn = ph_plane[index_cc];
               rho[index].re = amp*cos(tpi*rn);
               rho[index].im = -1.0*amp*sin(tpi*rn);
            }
         }
      }
      free(am_plane);
      free(ph_plane);
   }
   else {
      Parallel.ParallelError("indens: wrong mirror flag!", "stop");
   }
      
   // 0 and nq planes:
   if (nx1 == 0 || nx1 == nq){
      i = 0;
      for (j=0; j<ngy; ++j){
         for (k=0; k<ngz; ++k){
            index = (i*ngy + j)*ngz + k;
            if (j < nq){
               rn  = genrand_real();  // Amplitude
               amp = sqrt(-1.0*log(rn)*Pk[index]);
               rn  = genrand_real();  // Phase
               rho[index].re = amp*cos(tpi*rn);
               rho[index].im = amp*sin(tpi*rn);
            }
            else {
               index_cc = (i*ngy + MOD(ngrid-j,ngrid))*ngz + MOD(ngrid-k,ngrid);
               rho[index].re = rho[index_cc].re;
               rho[index].im = -1.0*rho[index_cc].im;
            }
         }
      }
               
      for (j=0; j<2; ++j){
         for (k=0; k<ngz; ++k){
            index = (i*ngy + j*nq)*ngz + k;
            if (k < nq){
               rn  = genrand_real();  // Amplitude
               amp = sqrt(-1.0*log(rn)*Pk[index]);
               rn  = genrand_real();  // Phase
               rho[index].re = genrand_real();
               rho[index].im = genrand_real();
            }
            else{
               index_cc = (i*ngy + MOD(ngrid-j*nq,ngrid))*ngz + MOD(ngrid-k,ngrid);
               rho[index].re = rho[index_cc].re;
               rho[index].im = -1.0*rho[index_cc].im;
            }
         }
      }
               
      // Real values:
      index = (i*ngy + nq)*ngz + 0;
      rho[index].im = 0.0;
      index = (i*ngy + 0)*ngz + 0;
      rho[index].im = 0.0;
      index = (i*ngy + 0)*ngz + nq;
      rho[index].im = 0.0;
      index = (i*ngy + nq)*ngz + nq;
      rho[index].im = 0.0;
      if (nx1 == 0) rho[0].re = 0.0;
   }
      
   return;
}


void indens2D(){
   const real tpi=2.0*pi;
   const integer ng=DataBase.ngrid;
   const integer nq=ng/2;
   long i, j, k, index, index_cc;
   int proc;
      
   Parallel.ParallelError("indens: 2D not implemented yet!", "stop");
      
   return;
}


void indens3D(){
   const real tpi=2.0*pi;
   const integer ng=DataBase.ngrid;
   const integer nq=ng/2;
   long i, j, k, index, index_cc;
   int proc;
      
   Parallel.ParallelError("indens: 3D not implemented yet!", "stop");

   return;
}


void indens_single(){
   const real tpi=2.0*pi;
   const integer ng=DataBase.ngrid;
   const integer nq=ng/2;
   long i, j, k, index, index_cc;
   long tmp_i, tmp_j;
   int proc;
   real amp, rn;
      
   proc = Parallel.GetMyPE();
   init_random(DataBase.seed, proc);
      
   // Set two halves of the domain:
   for (i=0; i<ng; ++i){
      for (j=0; j<ng; ++j){
         for (k=0; k<ng; ++k){
            index = (i*ng + j)*ng + k;
            if (i < nq){
               rn  = genrand_real();  // Amplitude
               amp = sqrt(-1.0*log(rn)*Pk[index]);
               rn  = genrand_real();  // Phase
               rho[index].re = amp*cos(tpi*rn);
               rho[index].im = amp*sin(tpi*rn);
            }
            else {
               index_cc = (MOD(ng-i,ng)*ng + MOD(ng-j,ng))*ng + MOD(ng-k,ng);
               rho[index].re = rho[index_cc].re;
               rho[index].im = -1.0*rho[index_cc].im;
            }
         }
      }
   }

   // 0 & nq planes:
   for (i=0; i<2; ++i) {
      tmp_i = i*nq;
      for (j=0; j<ng; ++j){
         for (k=0; k<ng; ++k){
            index = (tmp_i*ng + j)*ng + k;
            if (j < nq){
               rn  = genrand_real();  // Amplitude
               amp = sqrt(-1.0*log(rn)*Pk[index]);
               rn  = genrand_real();  // Phase
               rho[index].re = amp*cos(tpi*rn);
               rho[index].im = amp*sin(tpi*rn);
            }
            else {
               index_cc = (tmp_i*ng + MOD(ng-j,ng))*ng + MOD(ng-k,ng);
               rho[index].re = rho[index_cc].re;
               rho[index].im = -1.0*rho[index_cc].im;
            }
         }
      }
      for (j=0; j<2; ++j){
         tmp_j = j*nq;
         for (k=0; k<ng; ++k){
            index = (tmp_i*ng + tmp_j)*ng + k;
            if (k < nq){
               rn  = genrand_real();  // Amplitude
               amp = sqrt(-1.0*log(rn)*Pk[index]);
               rn  = genrand_real();  // Phase
               rho[index].re = amp*cos(tpi*rn);
               rho[index].im = amp*sin(tpi*rn);
            }
            else{
               index_cc = (tmp_i*ng + tmp_j)*ng + MOD(ng-k,ng);
               rho[index].re = rho[index_cc].re;
               rho[index].im = -1.0*rho[index_cc].im;
            }
         }
      }
   }
      
   // Real values:
   rho[0].re = 0.0;
   rho[0].im = 0.0;
   index = (nq*ng + 0)*ng + 0;
   rho[index].im = 0.0;
   index = (0*ng + nq)*ng + 0;
   rho[index].im = 0.0;
   index = (0*ng + 0)*ng + nq;
   rho[index].im = 0.0;
   index = (nq*ng + nq)*ng + 0;
   rho[index].im = 0.0;
   index = (0*ng + nq)*ng + nq;
   rho[index].im = 0.0;
   index = (nq*ng + 0)*ng + nq;
   rho[index].im = 0.0;
   index = (nq*ng + nq)*ng + nq;
   rho[index].im = 0.0;
      
   return;
}


void test_reality(MPI_Comm comm){
   long i;
   const integer ngrid=DataBase.ngrid;
      
   float test_zero, max_rho, min_rho, ave_rho, scal;
      
   min_rho = 1.0e30;
   max_rho = -1.0e30;
   ave_rho = 0.0;
   for (i=0; i<My_Ng; ++i) {
      ave_rho += (real)rho[i].re;
      if ((real)rho[i].re > max_rho) max_rho = (real)rho[i].re;
      if ((real)rho[i].re < min_rho) min_rho = (real)rho[i].re;
   }
   ave_rho   = ave_rho/My_Ng;
      
   MPI_Allreduce(MPI_IN_PLACE, &max_rho, 1, MPI_FLOAT, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, &min_rho, 1, MPI_FLOAT, MPI_MIN, comm);
   MPI_Allreduce(MPI_IN_PLACE, &ave_rho, 1, MPI_FLOAT, MPI_SUM, comm);
   if (MyPE == MasterPE) {
      std::cout << std::endl << "Min and max value of density in k space: "
            << min_rho << " " << max_rho << std::endl;
      ave_rho = ave_rho/NumPEs;
      std::cout << "Average value of density in k space: "
            << ave_rho << std::endl << std::flush;
   }

#if defined (FFTW2) || defined (FFTW3)
   MPI_Barrier(comm);
#ifdef FFTW2
   fftwnd_mpi(plan, 1, (fftw_complex *)rho, NULL, FFTW_TRANSPOSED_ORDER);
#endif
#ifdef FFTW3
   fftw_execute(plan);
#endif
   if (MyPE == MasterPE) std::cout << std::endl << "FFT DONE!" << std::endl;

#ifdef IPPLCREATE
   scal = std::pow(DataBase.box_size, 1.5f); // inverse FFT scaling
#else
   scal = pow(DataBase.box_size, 1.5f); // inverse FFT scaling
#endif
   for (i=0; i< My_Ng; ++i) rho[i].re = rho[i].re/scal;

   min_rho = 1.0e30;
   max_rho = -1.0e30;
   test_zero = 0.0;
   ave_rho = 0.0;
   for (i=0; i<My_Ng; ++i) {
      if (test_zero < fabs((real)rho[i].im))
         test_zero = fabs((real)rho[i].im);
      ave_rho += (real)rho[i].re;
      if ((real)rho[i].re > max_rho) max_rho = (real)rho[i].re;
      if ((real)rho[i].re < min_rho) min_rho = (real)rho[i].re;
   }
   ave_rho = ave_rho/My_Ng;


   MPI_Allreduce(MPI_IN_PLACE, &test_zero, 1, MPI_FLOAT, MPI_MAX, comm);	
   MPI_Allreduce(MPI_IN_PLACE, &max_rho, 1, MPI_FLOAT, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, &min_rho, 1, MPI_FLOAT, MPI_MIN, comm);
   MPI_Allreduce(MPI_IN_PLACE, &ave_rho, 1, MPI_FLOAT, MPI_SUM, comm);
   if (MyPE == MasterPE) {
      std::cout << std::endl << "Max value of the imaginary part of density is "
            << test_zero << std::endl;
      std::cout << "Min and max value of density in code units: "
            << min_rho << " " << max_rho << std::endl;
      ave_rho = ave_rho/NumPEs;
      std::cout << "Average value of density in code units: "
            << ave_rho << std::endl << std::flush;
   }

#else
   Parallel.ParallelError("Cannot test rho(r) without FFT", "no stop");
#endif
      
   return;
}


void solve_gravity(integer axis, MPI_Comm comm){ // 0=x, 1=y, 2=z
   long i, j, k, k_i, k_j, k_k, index;
   const real tpi=2.0*pi;
   const integer ngrid=DataBase.ngrid;
   const integer nq=ngrid/2;
   real kk, Grad, Green, scal;
   real xx, yy, zz;
   real min_F, max_F, ave_F, test_zero;
   my_fftw_complex temp;
      
   /* Multiply by the Green's function (-1/k^2), 
   and take the gradient in k-space (-ik_i) */
   switch (axis) {
      case 0:
         xx = 1.0;
         yy = 0.0;
         zz = 0.0;
         break;
      case 1:
         xx = 0.0;
         yy = 1.0;
         zz = 0.0;
         break;
      case 2:
         xx = 0.0;
         yy = 0.0;
         zz = 1.0;
         break;
   }
   for (i=0; i<ngx; ++i){
      k_i = i+nx1;
      if (k_i >= nq) {k_i = -MOD(ngrid-k_i,ngrid);}
      for (j=0; j<ngy; ++j){
         k_j = j+ny1;
         if (k_j >= nq) {k_j = -MOD(ngrid-k_j,ngrid);}
         for (k=0; k<ngz; ++k){
            k_k = k+nz1;
            if (k_k >= nq) {k_k = -MOD(ngrid-k_k,ngrid);}
            Grad = -tpi/ngrid*(k_i*xx+k_j*yy+k_k*zz);
            index = (i*ngy+j)*ngz+k;
            kk = tpi/ngrid*sqrt(k_i*k_i+k_j*k_j+k_k*k_k);
            Green = -1.0/(kk*kk);
            if (kk == 0.0) Green = 0.0;
            temp = rho[index];
            rho[index].re = -Grad * Green * temp.im;
            rho[index].im = Grad * Green * temp.re;
         }
      }
   }
   // Fix k=0 point:
   if (nx1 == 0 && ny1 == 0 && nz1 == 0) rho[0].re = 0.0;
      
#if defined (FFTW2)	|| defined (FFTW3)
   MPI_Barrier(comm);
#ifdef FFTW2	
   fftwnd_mpi(plan, 1, (fftw_complex *)rho, NULL, FFTW_NORMAL_ORDER);
#endif
#ifdef FFTW3
   fftw_execute(plan);
#endif
   if (MyPE == MasterPE) std::cout << std::endl << "FFT DONE!" << std::endl;

   /* compensation for inverse FFT scaling */
#ifdef IPPLCREATE
   scal = std::pow(DataBase.box_size, 1.5f);
#else
   scal = pow(DataBase.box_size, 1.5f);
#endif
   for (i=0; i<My_Ng; ++i) rho[i].re = rho[i].re/scal;
#else
   Parallel.ParallelError("FFT is not done", "no stop");
#endif
      
#ifdef TESTING
   OnClock("Tests");
   min_F = 1.0e30;
   max_F = -1.0e30;
   test_zero = 0.0;
   ave_F = 0.0;
   for (i=0; i<My_Ng; ++i) { 
      if (test_zero < fabs((real)rho[i].im))
         test_zero = fabs((real)rho[i].im);  
      ave_F += (real)rho[i].re;
      if ((real)rho[i].re > max_F) max_F = (real)rho[i].re;
      if ((real)rho[i].re < min_F) min_F = (real)rho[i].re;
   }
   ave_F = ave_F/My_Ng;

   MPI_Allreduce(MPI_IN_PLACE, &test_zero, 1, MPI_FLOAT, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, &max_F, 1, MPI_FLOAT, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, &min_F, 1, MPI_FLOAT, MPI_MIN, comm);
   MPI_Allreduce(MPI_IN_PLACE, &ave_F, 1, MPI_FLOAT, MPI_SUM, comm);
   if (MyPE == MasterPE) {
      std::cout << std::endl << "Max value of the imaginary part of the force is "
            << test_zero << std::endl;
      std::cout << "Min and max value of force in code units: "
            << min_F << " " << max_F << std::endl;
      ave_F = ave_F/NumPEs;
      std::cout << "Average value of force in code units: "
            << ave_F << std::endl << std::flush;
   }
   OffClock("Tests");
#endif
      
   return;
}


void set_particles(real z_in, real d_z, real ddot, integer axis, MPI_Comm comm){ // 0=x, 1=y, 2=z
   const integer ngrid=DataBase.ngrid;
   long i, j, k, index;
   real pos_0, xx, yy, zz, move, max_move, ave_move;
   my_fftw_complex *F_i;
      
   F_i = rho; /* After gravity solve, this complex array 
   holds F_i component of the force, not density */
      
   switch (axis) {
      case 0:
         xx = 1.0;
         yy = 0.0;
         zz = 0.0;
         break;
      case 1:
         xx = 0.0;
         yy = 1.0;
         zz = 0.0;
         break;
      case 2:
         xx = 0.0;
         yy = 0.0;
         zz = 1.0;
         break;
   }

   /* Move particles according to Zeldovich approximation 
   particles will be put in rho array; real array will 
   hold positions, imaginary will hold velocities */
#ifdef TESTING
   OnClock("Tests");
   max_move = 0.0;
   ave_move = 0.0;
   for (i=0; i<My_Ng; ++i) {
      move = fabs(d_z*F_i[i].re);
      ave_move += move;
      if (move > max_move) max_move = move;
   }
   ave_move = ave_move/My_Ng;

   MPI_Allreduce(MPI_IN_PLACE, &max_move, 1, MPI_FLOAT, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, &ave_move, 1, MPI_FLOAT, MPI_SUM, comm);
   if (MyPE == MasterPE) {
      if (axis == 0) {std::cout << "Max move in X: " << max_move << std::endl;}
      else if (axis == 1) {std::cout << "Max move in Y: " << max_move << std::endl;}
      else {std::cout << "Max move in Z: " << max_move << std::endl;}
      ave_move = ave_move/NumPEs;
      if (axis == 0) {std::cout << "Average move in X: " << ave_move << std::endl;}
      else if (axis == 1) {std::cout << "Average move in Y: " << ave_move << std::endl;}
      else {std::cout << "Average move in Z: " << ave_move << std::endl;}
   }
   OffClock("Tests");
#endif

   for (i=0; i<ngx; ++i){
      for (j=0; j<ngy; ++j){
         for (k=0; k<ngz; ++k){
            index = (i*ngy+j)*ngz+k;
            pos_0 = (i+nx1)*xx + (j+ny1)*yy + (k+nz1)*zz;
            F_i[index].im = -ddot*F_i[index].re;
            F_i[index].re = pos_0 - d_z*F_i[index].re;
            if ((float)F_i[index].re < 0.0) {F_i[index].re += (float)ngrid;}
            else if ((float)F_i[index].re > (float)ngrid) {F_i[index].re -= (float)ngrid;}
         }
      }
   }
      
   return;
}


#ifndef IPPLCREATE
void output(integer axis, float* pos, float* vel, MPI_Comm comm){ // 0=x, 1=y, 2=z
   long i, j;
   const char *fname;
   FILE *pfile;
   std::ofstream OutFile;
   MPI_Status status;
   MPI_Datatype MPI_FFTW_COMPLEX;
   my_fftw_complex *data;
   distribution_t ddd[1];
   int n_dist[3];
   int padding[3] = {0, 0, 0};
      
   // David's stuff:
   if (NumPEs > 1) {
      n_dist[0] = (int)DataBase.ngrid;
      n_dist[1] = (int)DataBase.ngrid;
      n_dist[2] = (int)DataBase.ngrid;
      distribution_init(comm, n_dist, padding, ddd, false);
      distribution_assert_commensurate(ddd);
      data  = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
      distribution_1_to_3(rho, data, ddd);
      for (i=0; i<My_Ng; ++i){rho[i] = data[i];}
      distribution_fini(ddd);
      free(data);
   }
   // End of David's stuff.
      
   if (DataBase.PrintFlag == 0){
      // No printing, particles are redirected somewhere
      for (i=0; i<My_Ng; ++i){
         pos[i] = (float)rho[i].re;
         vel[i] = (float)rho[i].im;
      }
   }
   else if (DataBase.PrintFlag == 1){
      // Serial print through the master processor into an ASCII file
      CreateMPI_FFTW_COMPLEX(&MPI_FFTW_COMPLEX);
      if (MyPE == MasterPE){
         fname = "part.dat";
         OutFile.open(fname, std::ios::out | std::ios::app);
         for (i=0; i<My_Ng; ++i) 
            {OutFile << (float)rho[i].re << " " << (float)rho[i].im << std::endl;}
         for (j=0; j<NumPEs; ++j){  // Get data from processor j:
            if (j == MasterPE) continue;
            MPI_Recv(rho, My_Ng, MPI_FFTW_COMPLEX, j, 101, comm, &status);
            for (i=0; i<My_Ng; ++i) 
               OutFile << (float)rho[i].re << " " << (float)rho[i].im << std::endl;
         }
         OutFile.close();
      }
      else{
         MPI_Send(rho, My_Ng, MPI_FFTW_COMPLEX, MasterPE, 101, comm);
      }
   }
   else if (DataBase.PrintFlag == 2){
      // Serial print through the master processor into a binary file
      CreateMPI_FFTW_COMPLEX(&MPI_FFTW_COMPLEX);
      if (MyPE == MasterPE) {
         fname = "part.bin";
         pfile = fopen(fname, "aw");
         fwrite(rho, sizeof(rho[0]), My_Ng, pfile);
         for (j=0; j<NumPEs; ++j){  // Get data from processor j:
            if (j == MasterPE) continue;
            MPI_Recv(rho, My_Ng, MPI_FFTW_COMPLEX, j, 101, comm, &status);
            for (i=0; i<My_Ng; ++i){fwrite(rho, sizeof(rho[0]), My_Ng, pfile);}
         }
         fclose(pfile);
      }
      else{
         MPI_Send(rho, My_Ng, MPI_FFTW_COMPLEX, MasterPE, 101, comm);
      }
   }
   else if (DataBase.PrintFlag == 3){
      // Parallel print -- each processor dumps a (binary) file
      
   }
      
   return;
}

void init_particles(float* pos_x, float* pos_y, float* pos_z, 
		    float* vel_x, float* vel_y, float* vel_z, 
		    Basedata& bdata, const char *tfName, MPI_Comm comm) {
   integer i;
   real d_z, ddot;
   double t1, t2;
   void (*indens)();  /* pointer to indens routine */
   
   // Preliminaries:
   StartMonitor();
   OnClock("Initialization");
   Parallel.InitMPI(&t1, comm);
   MyPE = Parallel.GetMyPE();
   MasterPE = Parallel.GetMasterPE();
   NumPEs = Parallel.GetNumPEs();
   DataBase.GetParameters(bdata, Parallel);
   Cosmology.SetParameters(DataBase, tfName);
   Parallel.DomainDecompose(DataBase.dim, DataBase.ngrid);
   if (NumPEs == 1){ indens = indens_single; }
   else {
      if (DataBase.dim == 1) { indens = indens1D; }
      else if (DataBase.dim == 2) { indens = indens2D; }
      else if (DataBase.dim == 3) { indens = indens3D; }
      else {Parallel.ParallelError("Main: Decomposition has to be 1-3", "stop");}
   }
      
   // Find what part of the global grid i got:
   nx1 = Parallel.GetNx1();
   nx2 = Parallel.GetNx2();
   ny1 = Parallel.GetNy1();
   ny2 = Parallel.GetNy2();
   nz1 = Parallel.GetNz1();
   nz2 = Parallel.GetNz2();
   ngx = nx2-nx1+1;
   ngy = ny2-ny1+1;
   ngz = nz2-nz1+1;
   if (DataBase.dim==3 && (ngz != ngy || ngx != ngz)) 
      Parallel.AbortMPI("Main: Domain incorrectly decomposed");
   if (DataBase.dim==2 && (ngx != ngy || ngz != DataBase.ngrid))
      Parallel.AbortMPI("Main: Domain incorrectly decomposed");
   if (DataBase.dim==1 && (ngz != ngy || ngz != DataBase.ngrid))
      Parallel.AbortMPI("Main: Domain incorrectly decomposed");
   My_Ng = ngx*ngy*ngz;

   // Allocate arays:
   Pk   = (real *)malloc(My_Ng*sizeof(real));
   rho  = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
      
   // Initialize FFTW:
#ifdef FFTW2
   const integer ngrid=DataBase.ngrid;
   plan = fftw3d_mpi_create_plan(comm, ngrid, ngrid, ngrid, 
         FFTW_BACKWARD, FFTW_ESTIMATE);
   fftwnd_mpi_local_sizes(plan, &lnx, &lxs, &lnyt, &lyst, &lsize);
   if (lnx != ngx){
      std::cout << "FFTW has different decomposition!";
      std::cout << "  Mine: " << ngx << "  FFTW's: " << lnx << std::endl;
      Parallel.AbortMPI(" ");
   }
#endif
#ifdef FFTW3
   const integer ngrid=DataBase.ngrid;
   fftw_mpi_init();
   plan = fftw_mpi_plan_dft_3d(ngrid, ngrid, ngrid, (fftw_complex *)rho, (fftw_complex *)rho, comm, 
         FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
   OffClock("Initialization");
      
   // Set P(k) array:
   OnClock("Creating Pk(x,y,z)");
   initspec();
   OffClock("Creating Pk(x,y,z)");
      
   // Growth factor of the initial redshift:
   d_z = Cosmology.GrowthFactor(DataBase.z_in);
   ddot = Cosmology.GrowthDeriv(DataBase.z_in);
   if (MyPE == MasterPE){
      printf("redshift: %f; growth factor = %f; ", DataBase.z_in, d_z);
      printf("derivative = %f\n\n", ddot);
   }
      
#ifdef TESTING
   // Check if density field is real:
   OnClock("Tests");
   indens();
   test_reality(comm);
   OffClock("Tests");
#endif
      
   // Set particles:
   for (i=0; i<Ndim; ++i){
      OnClock("Creating rho(x,y,z)");
      indens();          // Set density in k-space:
      OffClock("Creating rho(x,y,z)");
      OnClock("Poisson solve");
      solve_gravity(i, comm);  // Calculate i-th component of the force
      OffClock("Poisson solve");
      OnClock("Particle move");
      set_particles(DataBase.z_in, d_z, ddot, i, comm); // Zeldovich move
      OffClock("Particle move");
      
      OnClock("Output");
      if (i == 0) {output(i, pos_x, vel_x, comm);}   // Put particle data
      else if (i == 1) {output(i, pos_y, vel_y, comm);}
      else {output(i, pos_z, vel_z, comm);} 
      OffClock("Output");
   }
      
   // Clean all the allocations:
   free(Pk);
   free(rho);
   
   // ada if calling this Cosmology.~CosmoClass() a double free is reported AFTER mpi FINALIZE
   //Cosmology.~CosmoClass();

   DataBase.~GlobalStuff();
   Parallel.~ParallelTools();
#ifdef FFTW2
   fftwnd_mpi_destroy_plan(plan);
#endif
#ifdef FFTW3
   fftw_destroy_plan(plan);
#endif
      
   MPI_Barrier(comm);
   if (MyPE == MasterPE) PrintClockSummary(stdout);
   Parallel.FinalMPI(&t2);	

}
#endif

#ifdef IPPLCREATE
//template <class T, unsigned int Dim>
void init_particles_ippl(ParticleAttrib< Vektor<T, 3> > &pos, 
        ParticleAttrib< Vektor<T, 3> > &vel, 
        Basedata& bdata, const char *tfName, MPI_Comm comm) {

   integer i;
   real d_z, ddot;
   double t1, t2;
   void (*indens)();  // pointer to indens routine
   
   // Preliminaries:
   StartMonitor();
   OnClock("Initialization");
   Parallel.InitMPI(&t1, comm);
   MyPE = Parallel.GetMyPE();
   MasterPE = Parallel.GetMasterPE();
   NumPEs = Parallel.GetNumPEs();
   DataBase.GetParameters(bdata, Parallel);
   Cosmology.SetParameters(DataBase, tfName);
   Parallel.DomainDecompose(DataBase.dim, DataBase.ngrid);
   if (NumPEs == 1){ indens = indens_single; }
   else {
      if (DataBase.dim == 1) { indens = indens1D; }
      else if (DataBase.dim == 2) { indens = indens2D; }
      else if (DataBase.dim == 3) { indens = indens3D; }
      else {Parallel.ParallelError("Main: Decomposition has to be 1-3", "stop");}
   }
      
   // Find what part of the global grid i got:
   nx1 = Parallel.GetNx1();
   nx2 = Parallel.GetNx2();
   ny1 = Parallel.GetNy1();
   ny2 = Parallel.GetNy2();
   nz1 = Parallel.GetNz1();
   nz2 = Parallel.GetNz2();
   ngx = nx2-nx1+1;
   ngy = ny2-ny1+1;
   ngz = nz2-nz1+1;
   if (DataBase.dim==3 && (ngz != ngy || ngx != ngz)) 
      Parallel.AbortMPI("Main: Domain incorrectly decomposed");
   if (DataBase.dim==2 && (ngx != ngy || ngz != DataBase.ngrid))
      Parallel.AbortMPI("Main: Domain incorrectly decomposed");
   if (DataBase.dim==1 && (ngz != ngy || ngz != DataBase.ngrid))
      Parallel.AbortMPI("Main: Domain incorrectly decomposed");
   My_Ng = ngx*ngy*ngz;

   // Allocate arays:
   Pk   = (real *)malloc(My_Ng*sizeof(real));
   rho  = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
      
   // Initialize FFTW:
#ifdef FFTW2
   const integer ngrid=DataBase.ngrid;
   plan = fftw3d_mpi_create_plan(comm, ngrid, ngrid, ngrid, 
         FFTW_BACKWARD, FFTW_ESTIMATE);
   fftwnd_mpi_local_sizes(plan, &lnx, &lxs, &lnyt, &lyst, &lsize);
   if (lnx != ngx){
      std::cout << "FFTW has different decomposition!";
      std::cout << "  Mine: " << ngx << "  FFTW's: " << lnx << std::endl;
      Parallel.AbortMPI(" ");
   }
#endif
#ifdef FFTW3
   const integer ngrid=DataBase.ngrid;
   fftw_mpi_init();
   plan = fftw_mpi_plan_dft_3d(ngrid, ngrid, ngrid, (fftw_complex *)rho, (fftw_complex *)rho, comm, 
         FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
   OffClock("Initialization");
      
   // Set P(k) array:
   OnClock("Creating Pk(x,y,z)");
   initspec();
   OffClock("Creating Pk(x,y,z)");
      
   // Growth factor of the initial redshift:
   d_z = Cosmology.GrowthFactor(DataBase.z_in);
   ddot = Cosmology.GrowthDeriv(DataBase.z_in);
   if (MyPE == MasterPE){
      printf("redshift: %f; growth factor = %f; ", DataBase.z_in, d_z);
      printf("derivative = %f\n\n", ddot);
   }
      
#ifdef TESTING
   // Check if density field is real:
   OnClock("Tests");
   indens();
   test_reality(comm);
   OffClock("Tests");
#endif
      
   // Set particles:
   for (i=0; i<Ndim; ++i){
      OnClock("Creating rho(x,y,z)");
      indens();          // Set density in k-space:
      OffClock("Creating rho(x,y,z)");
      OnClock("Poisson solve");
      solve_gravity(i, comm);  // Calculate i-th component of the force
      OffClock("Poisson solve");
      OnClock("Particle move");
      set_particles(DataBase.z_in, d_z, ddot, i, comm); // Zeldovich move
      OffClock("Particle move");
      
      OnClock("Output");
      //MPI_Datatype MPI_FFTW_COMPLEX;
      //my_fftw_complex *data;
      //distribution_t ddd[1];
      //int n_dist[3];
      //int padding[3] = {0, 0, 0};

      //FIXME: WHAT THE HECK IS DAVIDS STUFF?
      // David's stuff: we dont do it here -> crashing

      for(size_t j=0; j<My_Ng; ++j){
          pos[j](i) = (double)rho[j].re;
          vel[j](i) = (double)rho[j].im;
      }
      OffClock("Output");
   }
      
   // Clean all the allocations:
   free(Pk);
   free(rho);
   
   // ada if calling this Cosmology.~CosmoClass() a double free is reported AFTER mpi FINALIZE
   //Cosmology.~CosmoClass();

   DataBase.~GlobalStuff();
   Parallel.~ParallelTools();
#ifdef FFTW2
   fftwnd_mpi_destroy_plan(plan);
#endif
#ifdef FFTW3
   fftw_destroy_plan(plan);
#endif
      
   MPI_Barrier(comm);
   if (MyPE == MasterPE) PrintClockSummary(stdout);
   Parallel.FinalMPI(&t2);	


}
#endif //IPPLCREATE

#ifdef INITIALIZER_NAMESPACE
}
#endif

#undef Ndim
#undef pi
