/* Cosmic initial conditions, currently for flat cosmologies only.

   Things (particle positions and velocities) are in grid units here. 

                        Zarija Lukic, May 2009
                           zarija@lanl.gov
*/

#include "Initializer.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>


#ifdef USENAMESPACE
namespace initializer {
#endif

void Initializer::CreateMPI_FFTW_COMPLEX(MPI_Datatype *MPI_FFTW_COMPLEX){
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


void Initializer::initspec(){
   const real tpi=2.0*pi;
   const real n_s = DataBase.n_s;
   const real sigma8=DataBase.Sigma_8;
   const real f_NL=DataBase.f_NL;
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
      std::cout << "\nsigma_8=" << s8 << ", target was " << sigma8 << std::endl;

   // For non-Gaussian initial conditions: P(k)=A*k^n, trasfer function will come later
   if (f_NL != 0.0){
      if (MyPE == MasterPE)
         std::cout << "Non-Gaussian initial conditions, f_NL=" << f_NL << std::endl;
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
               Pk[index] = norm*pow(kk, n_s);
            }
         }
      }
      return;
   }
   
   // For Gaussian initial conditions: P(k)=A*T^2*k^n
   if (MyPE == MasterPE)
         std::cout << "Gaussian initial conditions, f_NL=" << f_NL << std::endl;
   for (i=0; i<My_Ng; ++i) Pk[i] *= norm;

   return;
}

void Initializer::indens(){
	const integer ng=DataBase.ngrid;
	long i;
	int proc;
	double rn, rn1, rn2, scal;
	
	// Generate Gaussian white noise field:
	proc = Parallel.GetMyPE();
	init_random(DataBase.seed, proc);
	for (i=0; i<My_Ng; ++i){
		do {
			rn1 = -1 + 2*genrand_real();
			rn2 = -1 + 2*genrand_real();
			rn = rn1*rn1 + rn2*rn2;
		} while (rn > 1.0 || rn == 0);
		rho[i].re = rn2 * sqrt(-2.0*log(rn)/rn);
		rho[i].im = 0.0;
	}
	
	// Go to the k-space:
#if defined (FFTW2)	|| defined (FFTW3)
   MPI_Barrier(Parallel.GetMpiComm());
#ifdef FFTW2
   fftwnd_mpi(kplan, 1, (fftw_complex *)rho, NULL, FFTW_NORMAL_ORDER);
#endif
#ifdef FFTW3
   fftw_execute(kplan);
#endif
   if (MyPE == MasterPE) std::cout << "FFT DONE!" << std::endl;
#else
   Parallel.ParallelError("WARNING: FFT is not performed", "no stop");
#endif
	
	// Multiply by the power spectrum:
   scal = pow(ng, 1.5);
	for (i=0; i<My_Ng; ++i){
		rho[i].re *= sqrt(Pk[i])/scal;
		rho[i].im *= sqrt(Pk[i])/scal;
	}
   // Fix k=0 point:
   if (nx1 == 0 && ny1 == 0 && nz1 == 0) rho[0].re = 0.0;
	
	return;
}


void Initializer::test_reality(){
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
      
   MPI_Allreduce(MPI_IN_PLACE, &max_rho, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &min_rho, 1, MPI_FLOAT, MPI_MIN, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &ave_rho, 1, MPI_FLOAT, MPI_SUM, Parallel.GetMpiComm());
   if (MyPE == MasterPE) {
      std::cout << std::endl << "Min and max value of density in k space: "
            << min_rho << " " << max_rho << std::endl;
      ave_rho = ave_rho/NumPEs;
      std::cout << "Average value of density in k space: "
            << ave_rho << std::endl << std::flush;
   }

#if defined (FFTW2) || defined (FFTW3)
   MPI_Barrier(Parallel.GetMpiComm());
#ifdef FFTW2
   fftwnd_mpi(plan, 1, (fftw_complex *)rho, NULL, FFTW_TRANSPOSED_ORDER);
#endif
#ifdef FFTW3
   fftw_execute(plan);
#endif
   if (MyPE == MasterPE) std::cout << "FFT DONE!" << std::endl;

   scal = pow(DataBase.box_size, 1.5); // inverse FFT scaling
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

   MPI_Allreduce(MPI_IN_PLACE, &test_zero, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());	
   MPI_Allreduce(MPI_IN_PLACE, &max_rho, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &min_rho, 1, MPI_FLOAT, MPI_MIN, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &ave_rho, 1, MPI_FLOAT, MPI_SUM, Parallel.GetMpiComm());
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


void Initializer::gravity_force(integer axis){ // 0=x, 1=y, 2=z
   long i, j, k, k_i, k_j, k_k, index;
   const real tpi=2.0*pi;
   const integer ngrid=DataBase.ngrid;
   const integer nq=ngrid/2;
   real kk, Grad, Green, scal;
   real xx, yy, zz;
   float min_F, max_F, ave_F, test_zero;
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
            rho[index].im =  Grad * Green * temp.re;
         }
      }
   }
   // Fix k=0 point:
   if (nx1 == 0 && ny1 == 0 && nz1 == 0) rho[0].re = 0.0;
      
#if defined (FFTW2)	|| defined (FFTW3)
   MPI_Barrier(Parallel.GetMpiComm());
#ifdef FFTW2	
   fftwnd_mpi(plan, 1, (fftw_complex *)rho, NULL, FFTW_NORMAL_ORDER);
#endif
#ifdef FFTW3
   fftw_execute(plan);
#endif
   if (MyPE == MasterPE) std::cout << "FFT DONE!" << std::endl;

   /* compensation for inverse FFT scaling */
   scal = pow(DataBase.box_size, 1.5);
   for (i=0; i<My_Ng; ++i) rho[i].re = rho[i].re/scal;
#else
   Parallel.ParallelError("WARNING: FFT is not performed", "no stop");
#endif
      
#ifdef TESTING
   OnClock("Tests");
   min_F = 1.0e30;
   max_F = -1.0e30;
   test_zero = 0.0;
   ave_F = 0.0;
   for (i=0; i<My_Ng; ++i) { 
      if (test_zero < fabs((float)rho[i].im))
         test_zero = fabs((float)rho[i].im);  
      ave_F += (float)rho[i].re;
      if ((float)rho[i].re > max_F) max_F = (float)rho[i].re;
      if ((float)rho[i].re < min_F) min_F = (float)rho[i].re;
   }
   ave_F = ave_F/My_Ng;

   MPI_Allreduce(MPI_IN_PLACE, &test_zero, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &max_F, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &min_F, 1, MPI_FLOAT, MPI_MIN, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &ave_F, 1, MPI_FLOAT, MPI_SUM, Parallel.GetMpiComm());
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

   
void Initializer::gravity_potential(){ 
   long i, j, k, k_i, k_j, k_k, index;
   const real tpi=2.0*pi;
   const integer ngrid=DataBase.ngrid;
   const integer nq=ngrid/2;
   real kk, Green, scal;
   float min_phi, max_phi, ave_phi, test_zero;
   
   /* Multiply by the Green's function (-1/k^2) */
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
            kk = tpi/ngrid*sqrt(k_i*k_i+k_j*k_j+k_k*k_k);
            Green = -1.0/(kk*kk);
            if (kk == 0.0) Green = 0.0;
            rho[index].re *= Green;
            rho[index].im *= Green;
         }
      }
   }
   // Fix k=0 point:
   if (nx1 == 0 && ny1 == 0 && nz1 == 0) rho[0].re = 0.0;
   
#if defined (FFTW2)	|| defined (FFTW3)
   MPI_Barrier(Parallel.GetMpiComm());
#ifdef FFTW2	
   fftwnd_mpi(plan, 1, (fftw_complex *)rho, NULL, FFTW_NORMAL_ORDER);
#endif
#ifdef FFTW3
   fftw_execute(plan);
#endif
   if (MyPE == MasterPE) std::cout << "FFT DONE!" << std::endl;
   
   /* compensation for inverse FFT scaling */
   scal = pow(DataBase.box_size, 1.5);
   for (i=0; i<My_Ng; ++i) rho[i].re = rho[i].re/scal;
#else
   Parallel.ParallelError("WARNING: FFT is not performed", "no stop");
#endif
   
#ifdef TESTING
   OnClock("Tests");
   min_phi = 1.0e30;
   max_phi = -1.0e30;
   test_zero = 0.0;
   ave_phi = 0.0;
   for (i=0; i<My_Ng; ++i) { 
      if (test_zero < fabs((float)rho[i].im))
         test_zero = fabs((float)rho[i].im);  
      ave_phi += (float)rho[i].re;
      if ((float)rho[i].re > max_phi) max_phi = (float)rho[i].re;
      if ((float)rho[i].re < min_phi) min_phi = (float)rho[i].re;
   }
   ave_phi = ave_phi/My_Ng;
   
   MPI_Allreduce(MPI_IN_PLACE, &test_zero, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &max_phi, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &min_phi, 1, MPI_FLOAT, MPI_MIN, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &ave_phi, 1, MPI_FLOAT, MPI_SUM, Parallel.GetMpiComm());
   if (MyPE == MasterPE) {
      std::cout << std::endl << "Max value of the imaginary part of the potential is "
      << test_zero << std::endl;
      std::cout << "Min and max value of the potential in code units: "
      << min_phi << " " << max_phi << std::endl;
      ave_phi = ave_phi/NumPEs;
      std::cout << "Average value of potential in code units: "
      << ave_phi << std::endl << std::flush;
   }
   OffClock("Tests");
#endif
   
   return;
}

   
void Initializer::non_gaussian(integer axis){ // 0=x, 1=y, 2=z 
   const real f_NL=DataBase.f_NL;
   const integer ngrid=DataBase.ngrid;
   const integer nq=ngrid/2;
   const real tpi=2.0*pi;
   const real tpiL=tpi/DataBase.box_size; // k0, physical units
   long i, j, k, k_i, k_j, k_k, index;
   float ave_phi2;
   float xx, yy, zz;
   real Grad, kk, trans_f, fff;
   my_fftw_complex temp, *phi;

   phi = rho; /* After gravity solve this array holds potential */
   // Find average value for Phi^2:
   ave_phi2 = 0.0;
   for (i=0; i<My_Ng; ++i) {ave_phi2 += (float)(phi[i].re*phi[i].re);}
   ave_phi2 = ave_phi2/My_Ng;
   MPI_Allreduce(MPI_IN_PLACE, &ave_phi2, 1, MPI_FLOAT, MPI_SUM, Parallel.GetMpiComm());
   ave_phi2 = ave_phi2/NumPEs;

   // Add f_NL:
   fff = 1.0;
   for (i=0; i<My_Ng; ++i) {phi[i].re += f_NL*(phi[i].re*phi[i].re - ave_phi2)*fff;}
   
   // Go back to k-space (no forward scaling):
#if defined (FFTW2)	|| defined (FFTW3)
   MPI_Barrier(Parallel.GetMpiComm());
#ifdef FFTW2
   fftwnd_mpi(kplan, 1, (fftw_complex *)phi, NULL, FFTW_NORMAL_ORDER);
#endif
#ifdef FFTW3
   fftw_execute(kplan);
#endif
   if (MyPE == MasterPE) std::cout << "FFT DONE!" << std::endl;
#else
   Parallel.ParallelError("WARNING: FFT is not performed", "no stop");
#endif
   
   // Add transfer function and gradient in k-space:
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
            index = (i*ngy+j)*ngz+k;
            kk = tpiL*sqrt(k_i*k_i+k_j*k_j+k_k*k_k);
            trans_f = Cosmology.TransferFunction(kk);
            Grad = -tpi/ngrid*(k_i*xx+k_j*yy+k_k*zz);
            temp = phi[index];
            phi[index].re = -Grad * trans_f * temp.im;
            phi[index].im =  Grad * trans_f * temp.re;
         }
      }
   }
   if (nx1 == 0 && ny1 == 0 && nz1 == 0) phi[0].re = 0.0; // fix k=0
   
   // FFT to real space to get the force:
#if defined (FFTW2)	|| defined (FFTW3)
   MPI_Barrier(Parallel.GetMpiComm());
#ifdef FFTW2	
   fftwnd_mpi(plan, 1, (fftw_complex *)phi, NULL, FFTW_NORMAL_ORDER);
#endif
#ifdef FFTW3
   fftw_execute(plan);
#endif
   if (MyPE == MasterPE) std::cout << "FFT DONE!" << std::endl;
   
   real scal = pow(ngrid, 3.0);  /* inverse FFT normalization */
   for (i=0; i<My_Ng; ++i) phi[i].re = phi[i].re/scal;
#else
   Parallel.ParallelError("WARNING: FFT is not performed", "no stop");
#endif
   
   return;
}

   
void Initializer::set_particles(real z_in, real d_z, real ddot, integer axis){ // 0=x, 1=y, 2=z
   const integer ngrid=DataBase.ngrid;
   long i, j, k, index;
   real pos_0, xx, yy, zz;
   float  move, max_move, ave_move;
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
      move = fabs(d_z*(float)F_i[i].re);
      ave_move += move;
      if (move > max_move) max_move = move;
   }
   ave_move = ave_move/My_Ng;

   MPI_Allreduce(MPI_IN_PLACE, &max_move, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());
   MPI_Allreduce(MPI_IN_PLACE, &ave_move, 1, MPI_FLOAT, MPI_SUM, Parallel.GetMpiComm());
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
         }
      }
   }
      
   return;
}


void Initializer::output(integer axis, real* pos, real* vel){ // 0=x, 1=y, 2=z
   const integer ngrid=DataBase.ngrid;
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
      
   /*// David's stuff:
   if (NumPEs > 1) {
      n_dist[0] = (int)DataBase.ngrid;
      n_dist[1] = (int)DataBase.ngrid;
      n_dist[2] = (int)DataBase.ngrid;
      distribution_init(Parallel.GetMpiComm(), n_dist, padding, ddd, false);
      distribution_assert_commensurate(ddd);
      data  = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
      distribution_1_to_3(rho, data, ddd);
      for (i=0; i<My_Ng; ++i){rho[i] = data[i];}
      distribution_fini(ddd);
      free(data);
   }
   // End of David's stuff.*/
   
   for (i=0; i<My_Ng; ++i){
      pos[i] = (real)rho[i].re;
      vel[i] = (real)rho[i].im;
   }

   /*if (DataBase.PrintFormat == 0){
      // No printing, particles are redirected somewhere
   }
   else if (DataBase.PrintFormat == 1){
      // Serial print through the master processor into an ASCII file
      CreateMPI_FFTW_COMPLEX(&MPI_FFTW_COMPLEX);
      if (MyPE == MasterPE){
         fname = "part.dat";
         OutFile.open(fname, std::ios::out | std::ios::app);
         for (i=0; i<My_Ng; ++i) 
            {OutFile << (float)rho[i].re << " " << (float)rho[i].im << std::endl;}
         for (j=0; j<NumPEs; ++j){  // Get data from processor j:
            if (j == MasterPE) continue;
            MPI_Recv(rho, My_Ng, MPI_FFTW_COMPLEX, j, 101, Parallel.GetMpiComm(), &status);
            for (i=0; i<My_Ng; ++i) 
               OutFile << (float)rho[i].re << " " << (float)rho[i].im << std::endl;
         }
         OutFile.close();
      }
      else{
         MPI_Send(rho, My_Ng, MPI_FFTW_COMPLEX, MasterPE, 101, Parallel.GetMpiComm());
      }
   }
   else if (DataBase.PrintFormat == 2){
      // Serial print through the master processor into a binary file
      CreateMPI_FFTW_COMPLEX(&MPI_FFTW_COMPLEX);
      if (MyPE == MasterPE) {
         fname = "part.bin";
         pfile = fopen(fname, "aw");
         fwrite(rho, sizeof(rho[0]), My_Ng, pfile);
         for (j=0; j<NumPEs; ++j){  // Get data from processor j:
            if (j == MasterPE) continue;
            MPI_Recv(rho, My_Ng, MPI_FFTW_COMPLEX, j, 101, Parallel.GetMpiComm(), &status);
            for (i=0; i<My_Ng; ++i){fwrite(rho, sizeof(rho[0]), My_Ng, pfile);}
         }
         fclose(pfile);
      }
      else{
         MPI_Send(rho, My_Ng, MPI_FFTW_COMPLEX, MasterPE, 101, Parallel.GetMpiComm());
      }
   }
   else if (DataBase.PrintFormat == 3){
      // Parallel print -- each processor dumps a (binary) file
      
   }*/
      
   return;
}


void Initializer::exchange_buffers(real *pos, real *vel){
   long j, k, index;
   int recv_proc, send_proc;
   MPI_Status status;
   
   // Data to be sent -- leftmost plane in x:
   for (j=0; j<ngy; ++j){
      for (k=0; k<ngz; ++k){
         index = j*ngz+k;
         buf_pos[index] = pos[index];
         buf_vel[index] = vel[index];
      }
   }
   
   if (NumPEs > 1) { // Send data to process on the left, get data from
                     // process on the right side. Boundaries are periodic. 
      send_proc = MyPE-1;
      recv_proc = MyPE+1;
      if(MyPE == 0) send_proc = NumPEs-1;
      if(MyPE == NumPEs-1) recv_proc = 0;
      
      MPI_Sendrecv_replace(buf_pos, ngy*ngz, MY_MPI_REAL, send_proc, 0, 
                           recv_proc, 0, Parallel.GetMpiComm(), &status);
      MPI_Sendrecv_replace(buf_vel, ngy*ngz, MY_MPI_REAL, send_proc, 0, 
                           recv_proc, 0, Parallel.GetMpiComm(), &status);
   }
   return;
}


void Initializer::inject_neutrinos(real *pos, real *vel, integer axis){
   const integer ngrid=DataBase.ngrid;
   long i, j, k, index, iplus, jplus, kplus;
   int nu_pairs, nn;
   int jwrap, kwrap;
   real xx, yy, zz, pos0;
   real weighted_pos, weighted_vel;
   
   nu_pairs=DataBase.nu_pairs;
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
   
   for (i=0; i<ngx-1; ++i){
      for (j=0; j<ngy; ++j){
         if (j == ngy-1) {jplus=0; jwrap=1;}
         else {jplus=j+1; jwrap=0;}
         for (k=0; k<ngz; ++k){
            if (k = ngz-1) {kplus=0; kwrap=1;}
            else {kplus=k+1; kwrap=0;}
            
            index = ((i+1)*ngy+j)*ngz+k;
            pos0 = (i+1+nx1)*xx + (j+ny1)*yy + (k+nz1)*zz;
            weighted_pos = (pos[index]-pos0);
            weighted_vel = vel[index];
            index = (i*ngy+jplus)*ngz+k;
            pos0 = (i+nx1)*xx + (jplus+ny1)*yy + (k+nz1)*zz;
            weighted_pos += (pos[index]-pos0);
            weighted_vel += vel[index];
            index = (i*ngy+j)*ngz+kplus;
            pos0 = (i+nx1)*xx + (j+ny1)*yy + (kplus+nz1)*zz;
            weighted_pos += (pos[index]-pos0);
            weighted_vel += vel[index];
            index = ((i+1)*ngy+jplus)*ngz+k;
            pos0 = (i+1+nx1)*xx + (jplus+ny1)*yy + (k+nz1)*zz;
            weighted_pos += (pos[index]-pos0);
            weighted_vel += vel[index];
            index = ((i+1)*ngy+j)*ngz+kplus;
            pos0 = (i+1+nx1)*xx + (j+ny1)*yy + (kplus+nz1)*zz;
            weighted_pos += (pos[index]-pos0);
            weighted_vel += vel[index];
            index = (i*ngy+jplus)*ngz+kplus;
            pos0 = (i+nx1)*xx + (jplus+ny1)*yy + (kplus+nz1)*zz;
            weighted_pos += (pos[index]-pos0);
            weighted_vel += vel[index];
            index = ((i+1)*ngy+jplus)*ngz+kplus;
            pos0 = (i+1+nx1)*xx + (jplus+ny1)*yy + (kplus+nz1)*zz;
            weighted_pos += (pos[index]-pos0);
            weighted_vel += vel[index];
            index = (i*ngy+j)*ngz+k;
            pos0 = (i+nx1)*xx + (j+ny1)*yy + (k+nz1)*zz;
            weighted_pos += (pos[index]-pos0);
            weighted_vel += vel[index];
            weighted_pos = pos0+0.5+(jwrap*yy+kwrap*zz)*(real)ngrid + weighted_pos/8.0;
            weighted_vel /= 8.0;
            
            for (nn=0; nn<2*nu_pairs; ++nn){
               pos[My_Ng+index*2*nu_pairs+nn] = weighted_pos;
               vel[My_Ng+index*2*nu_pairs+nn] = weighted_vel;
            }
         }
      }
   }
   
   i = ngx-1;
   iplus = 0;
   for (j=0; j<ngy; ++j){
      if (j == ngy-1) {jplus = 0;}
      else {jplus = j+1;}
      for (k=0; k<ngz; ++k){
         if (k = ngz-1) {kplus = 0;}
         else {kplus = k+1;}
         
         index = (iplus*ngy+j)*ngz+k;
         pos0 = (iplus+nx1)*xx + (j+ny1)*yy + (k+nz1)*zz;
         weighted_pos = (buf_pos[index]-pos0);
         weighted_vel = buf_vel[index];
         index = (i*ngy+jplus)*ngz+k;
         pos0 = (i+nx1)*xx + (jplus+ny1)*yy + (k+nz1)*zz;
         weighted_pos += (pos[index]-pos0);
         weighted_vel += vel[index];
         index = (i*ngy+j)*ngz+kplus;
         pos0 = (i+nx1)*xx + (j+ny1)*yy + (kplus+nz1)*zz;
         weighted_pos += (pos[index]-pos0);
         weighted_vel += vel[index];
         index = (iplus*ngy+jplus)*ngz+k;
         pos0 = (iplus+nx1)*xx + (jplus+ny1)*yy + (k+nz1)*zz;
         weighted_pos += (buf_pos[index]-pos0);
         weighted_vel += buf_vel[index];
         index = (iplus*ngy+j)*ngz+kplus;
         pos0 = (iplus+nx1)*xx + (j+ny1)*yy + (kplus+nz1)*zz;
         weighted_pos += (buf_pos[index]-pos0);
         weighted_vel += buf_vel[index];
         index = (i*ngy+jplus)*ngz+kplus;
         pos0 = (i+nx1)*xx + (jplus+ny1)*yy + (kplus+nz1)*zz;
         weighted_pos += (pos[index]-pos0);
         weighted_vel += vel[index];
         index = (iplus*ngy+jplus)*ngz+kplus;
         pos0 = (iplus+nx1)*xx + (jplus+ny1)*yy + (kplus+nz1)*zz;
         weighted_pos += (buf_pos[index]-pos0);
         weighted_vel += buf_vel[index];
         index = (i*ngy+j)*ngz+k;
         pos0 = (i+nx1)*xx + (j+ny1)*yy + (k+nz1)*zz;
         weighted_pos += (pos[index]-pos0);
         weighted_vel += vel[index];
         weighted_pos = pos0+0.5+(xx*ngrid) + weighted_pos/8.0;
         weighted_vel /= 8.0;
         
         for (nn=0; nn<2*nu_pairs; ++nn){
            pos[My_Ng+index*2*nu_pairs+nn] = weighted_pos;
            vel[My_Ng+index*2*nu_pairs+nn] = weighted_vel;
         }
      }
   }
   
   return;
}


void Initializer::thermal_vel(real *vx, real *vy, real *vz){
   long i;
   int nu_pairs;
   real tx, ty, tz;
   
   nu_pairs=DataBase.nu_pairs;
   for (i=My_Ng; i<(2*nu_pairs+1)*My_Ng; i+=2) {
      Cosmology.FDVelocity(tx, ty, tz);
      vx[i]   += tx;
      vx[i+1] = -vx[i];
      vy[i]   += ty;
      vy[i+1] = -vy[i];
      vz[i]   += tz;
      vz[i+1] = -vz[i];
   }
   
   return;
}


void Initializer::apply_periodic(real *pos){
   const integer ngrid=DataBase.ngrid;
   long i, np;
   
   np = My_Ng;
   if (DataBase.Omega_nu > 0.0) {
      np = (2*DataBase.nu_pairs+1)*My_Ng;
   }
   for (i=0; i<np; ++i){
      if (pos[i] < 0.0) {pos[i] += (real)ngrid;}
      else if (pos[i] >= ngrid) {pos[i] -= (real)ngrid;}
   }
   return;
}


void Initializer::init_particles(real* pos_x, real* pos_y, real* pos_z, 
		    real* vel_x, real* vel_y, real* vel_z, 
		    InputParser& par, const char *tfName) {
   integer i;
   real d_z, ddot, f_NL;
   double t1, t2;
   
   // Preliminaries:
   OnClock("Initialization");
   Parallel.InitMPI(&t1);
   MyPE = Parallel.GetMyPE();
   MasterPE = Parallel.GetMasterPE();
   NumPEs = Parallel.GetNumPEs();

   DataBase.GetParameters(par, Parallel);
   Cosmology.SetParameters(DataBase, tfName);
   Parallel.DomainDecompose(DataBase.dim, DataBase.ngrid);
      
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
   if (DataBase.dim==2 && ngz != DataBase.ngrid)
      Parallel.AbortMPI("Main: Domain incorrectly decomposed");
   if (DataBase.dim==1 && (ngz != ngy || ngz != DataBase.ngrid))
      Parallel.AbortMPI("Main: Domain incorrectly decomposed");
   My_Ng = ngx*ngy*ngz;

   // Allocate arays:
   Pk   = (real *)malloc(My_Ng*sizeof(real));
   rho  = (my_fftw_complex *)malloc(My_Ng*sizeof(my_fftw_complex));
      
   // Initialize FFTW:
#ifdef FFTW2
   if (DataBase.dim != 1) 
      Parallel.AbortMPI("FFTW works only with 1D decomposition!");
   const integer ngrid=DataBase.ngrid;
   plan = fftw3d_mpi_create_plan(Parallel.GetMpiComm(), ngrid, ngrid, ngrid, 
                                 FFTW_BACKWARD, FFTW_ESTIMATE);
   fftwnd_mpi_local_sizes(plan, &lnx, &lxs, &lnyt, &lyst, &lsize);
   if (lnx != ngx){
       std::cout << "FFTW has different decomposition!";
       std::cout << "  Mine: " << ngx << "  FFTW's: " << lnx << std::endl;
       Parallel.AbortMPI(" ");
   }
   kplan = fftw3d_mpi_create_plan(Parallel.GetMpiComm(), ngrid, ngrid, ngrid, 
                                  FFTW_FORWARD, FFTW_ESTIMATE);
#endif
#ifdef FFTW3
   if (DataBase.dim != 1) 
      Parallel.AbortMPI("FFTW works only with 1D decomposition!");
   const integer ngrid=DataBase.ngrid;
   fftw_mpi_init();
   plan = fftw_mpi_plan_dft_3d(ngrid, ngrid, ngrid, (fftw_complex *)rho, 
                               (fftw_complex *)rho, Parallel.GetMpiComm(), 
                               FFTW_BACKWARD, FFTW_ESTIMATE);
   kplan = fftw_mpi_plan_dft_3d(ngrid, ngrid, ngrid, (fftw_complex *)rho, 
                               (fftw_complex *)rho, Parallel.GetMpiComm(), 
                               FFTW_FORWARD, FFTW_ESTIMATE);
#endif
   OffClock("Initialization");

   // Set P(k) array:
   OnClock("Creating Pk(x,y,z)");
   initspec();
   OffClock("Creating Pk(x,y,z)");
      
   // Growth factor for the initial redshift:
   Cosmology.GrowthFactor(DataBase.z_in, &d_z, &ddot);
   if (MyPE == MasterPE){
       std::cout << "redshift: " << DataBase.z_in << "; growth factor=" << d_z << "; ";
       std::cout << "derivative=" << ddot << "\n" << std::endl;
   }
      
#ifdef TESTING
   // Check if density field is real:
   OnClock("Tests");
   indens();
   test_reality();
   OffClock("Tests");
#endif
      
   // Set particles:
   f_NL=DataBase.f_NL;
   for (i=0; i<Ndim; ++i){
      OnClock("Creating rho(x,y,z)");
      indens();          // Set density in k-space:
      OffClock("Creating rho(x,y,z)");
      
      if (f_NL == 0.0){ // Solve for the force immediately 
         OnClock("Poisson solve");
         gravity_force(i);  // Calculate i-th component of the force
         OffClock("Poisson solve");
      }
      else { // Solve for Phi, apply f_NL, apply T(k), solve for the force
         OnClock("Poisson solve");
         gravity_potential();  // Calculate gravitational potential
         OffClock("Poisson solve");
         OnClock("Non-Gaussianity");
         non_gaussian(i);  // Create non-gaussian force
         OffClock("Non-Gaussianity");
      }

      OnClock("Particle move");
      set_particles(DataBase.z_in, d_z, ddot, i); // Zeldovich move
      OffClock("Particle move");
      
      if (i == 0) {
         OnClock("Output");
         output(i, pos_x, vel_x);
         OffClock("Output");
         if (DataBase.Omega_nu > 0.0){
            OnClock("Neutrinos");
            buf_pos = (real *)malloc(ngy*ngz*sizeof(real));
            buf_vel = (real *)malloc(ngy*ngz*sizeof(real));
            exchange_buffers(pos_x, vel_x);
            inject_neutrinos(pos_x, vel_x, i);
            OffClock("Neutrinos");
         }
         apply_periodic(pos_x);
      }   // Put particle data
      else if (i == 1) {
         OnClock("Output");
         output(i, pos_y, vel_y);
         OffClock("Output");
         if (DataBase.Omega_nu > 0.0){
            OnClock("Neutrinos");
            exchange_buffers(pos_y, vel_y);
            inject_neutrinos(pos_y, vel_y, i);
            OffClock("Neutrinos");
         }
         apply_periodic(pos_y);
      }
      else {
         OnClock("Output");
         output(i, pos_z, vel_z);
         OffClock("Output");
         if (DataBase.Omega_nu > 0.0){
            OnClock("Neutrinos");
            exchange_buffers(pos_z, vel_z);
            inject_neutrinos(pos_z, vel_z, i);
            free(buf_pos);
            free(buf_vel);
            OffClock("Neutrinos");
         }
         apply_periodic(pos_z);
      }
   }
   // Add thermal velocity to neutrino particles:
   if (DataBase.Omega_nu > 0.0) thermal_vel(vel_x, vel_y, vel_z);
   
   // Clean all allocations:
   free(Pk);
   free(rho);
#ifdef FFTW2
   fftwnd_mpi_destroy_plan(plan);
   fftwnd_mpi_destroy_plan(kplan);
#endif
#ifdef FFTW3
   fftw_destroy_plan(plan);
   fftw_destroy_plan(kplan);
#endif
      
   MPI_Barrier(Parallel.GetMpiComm());
   //if (MyPE == MasterPE) PrintClockSummary(stdout);
   Parallel.FinalMPI(&t2);	

}

#ifdef USENAMESPACE
}
#endif
