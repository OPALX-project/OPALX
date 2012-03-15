/*
   Initializer:
   Cosmology.cpp

      This file has implementation for routines:

         GetParameters(DataBase)
      which gets cosmology parameters from the
      DataBase (defined in DataBase.h), and stores them 
      with this class as well. 

         Sigma_r(R, AA)
      which returns mass variance on a scale R, of the 
      linear density field with normalization AA.

         TransferFunction(k)
      which returns transfer function for the mode k in 
      Fourier space.

         GrowthFactor(z)
      which returns the linear growth factor at redshift z.

         GrowthDeriv(z)
      which returns the derivative of the linear growth factor 
      at redshift z.

                        Zarija Lukic, February 2009
                           zarija@lanl.gov
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "TypesAndDefs.h"
#include "DataBase.h"
#include "Cosmology.h"

#ifdef INITIALIZER_NAMESPACE
namespace initializer {
#endif

void CosmoClass::SetParameters(GlobalStuff& DataBase, const char *tfName){
   std::ifstream inputFile;
   integer i, ln;
   real tmp, tfcdm, tfbar, norm;
   real akh1, akh2, alpha;
      
   Omega_m = DataBase.Omega_m;
   Omega_bar = DataBase.Omega_bar;
   h = DataBase.Hubble;
   n_s = DataBase.n_s;
   TFFlag = DataBase.TFFlag;
   w_de = DataBase.w_de;
      
   cobe_temp=2.728;  // COBE/FIRAS CMB temperature in K
   tt=cobe_temp/2.7 * cobe_temp/2.7;
      
   if (TFFlag == 0) { // CMBFast transfer function
               //inputFile.open("cmb.tf", std::ios::in);
      inputFile.open(tfName, std::ios::in);
      if (!inputFile){
         printf("Cannot open CMBFast file!");
         abort();
      }
      ln = 0;
      while (! inputFile.eof()) {
         inputFile >> tmp >> tfcdm >> tfbar >> tmp >> tmp;
                        //std::cout << tfcdm << " " << tfbar << std::endl;
         ++ln;
      }
      table_size = ln-1; // It also reads newline...
      table_kk = (real *)malloc(table_size*sizeof(real));
      table_tf = (real *)malloc(table_size*sizeof(real));
      inputFile.close();
      inputFile.open("cmb.tf", std::ios::in);
               //inputFile.seekg(0, std::ios::beg);
      inputFile >> table_kk[0] >> tfcdm >> tfbar >> tmp >> tmp;
      table_tf[0] = tfbar*Omega_bar/Omega_m + 
            tfcdm*(Omega_m-Omega_bar)/Omega_m;
      norm = table_tf[0];
      for (i=1; i<table_size; ++i){
         inputFile >> table_kk[i] >> tfcdm >> tfbar >> tmp >> tmp;
         table_tf[i] = tfbar*Omega_bar/Omega_m + 
               tfcdm*(Omega_m-Omega_bar)/Omega_m;
         table_tf[i] = table_tf[i]/norm; // Normalization
                        //std::cout << table_kk[0] << " " << table_tf[0] << std::endl;
      }
      last_k = table_kk[table_size-1];
      inputFile.close();
      tan_f=&CosmoClass::cmbfast;
   }
   else if (TFFlag == 1) { // Klypin-Holtzmann transfer function
      akh1=pow(46.9*Omega_m*h*h, 0.670)*(1.0+pow(32.1*Omega_m*h*h, -0.532));
      akh2=pow(12.0*Omega_m*h*h, 0.424)*(1.0+pow(45.0*Omega_m*h*h, -0.582));
      alpha=pow(akh1, -1.0*Omega_bar/Omega_m)*pow(akh2, pow(-1.0*Omega_bar/Omega_m, 3.0));
      kh_tmp = Omega_m*h*sqrt(alpha)*pow(1.0-Omega_bar/Omega_m, 0.6);
      tan_f=&CosmoClass::klypin_holtzmann;
      last_k = 10.0;
   }
   else if (TFFlag == 2) { // Hu-Sugiyama transfer function
      akh1=pow(46.9*Omega_m*h*h, 0.670)*(1.0+pow(32.1*Omega_m*h*h, -0.532));
      akh2=pow(12.0*Omega_m*h*h, 0.424)*(1.0+pow(45.0*Omega_m*h*h, -0.582));
      alpha=pow(akh1, -1.0*Omega_bar/Omega_m)*pow(akh2, pow(-1.0*Omega_bar/Omega_m, 3.0));
      kh_tmp = Omega_m*h*sqrt(alpha);
      tan_f = &CosmoClass::hu_sugiyama;
      last_k = 10.0;
   }
   else if (TFFlag == 3) { // Peacock-Dodds transfer function
      kh_tmp = Omega_m*h*exp(-2.0*Omega_bar);
      tan_f = &CosmoClass::peacock_dodds;
      last_k = 10.0;
   }
   else if (TFFlag == 4) { // BBKS transfer function
      kh_tmp = Omega_m*h;
      tan_f = &CosmoClass::bbks;
      last_k = 10.0;
   }
   else{
      printf("Cosmology::SetParameters: TFFlag has to be 0, 1, 2, 3 or 4!\n");
      abort();
   }
      
   return;
}

CosmoClass::~CosmoClass(){
   free(table_kk);
   free(table_tf);
}

#define pi M_PI

/* Transfer functions */

real CosmoClass::cmbfast(real k){
   real t_f;
   t_f = interpolate(table_kk, table_tf, table_size, k);
   return(t_f);
}

real CosmoClass::klypin_holtzmann(real k){
   real qkh, t_f;
      
   if (k == 0.0) return(0.0);
   qkh = k*tt/kh_tmp;
      /* NOTE: the following line has 0/0 for k=0.
   This was taken care of at the beginning of the routine. */
   t_f = log(1.0+2.34*qkh)/(2.34*qkh) * pow(1.0+13.0*qkh+
         pow(10.5*qkh, 2.0)+pow(10.4*qkh, 3.0)+pow(6.51*qkh, 4.0), -0.25);
   return(t_f);
}

real CosmoClass::hu_sugiyama(real k){
   real qkh, t_f;
      
   if (k == 0.0) return(0.0);
   qkh = k*tt/kh_tmp;
      /* NOTE: the following line has 0/0 for k=0.
   This was taken care of at the beginning of the routine. */
   t_f = log(1.0+2.34*qkh)/(2.34*qkh) * pow(1.0+3.89*qkh+
         pow(16.1*qkh, 2.0)+pow(5.46*qkh, 3.0)+pow(6.71*qkh, 4.0), -0.25);
   return(t_f);
}

real CosmoClass::peacock_dodds(real k){
   real qkh, t_f;
      
   if (k == 0.0) return(0.0);
   qkh = k/kh_tmp;
      /* NOTE: the following line has 0/0 for k=0.
   This was taken care of at the beginning of the routine. */
   t_f = log(1.0+2.34*qkh)/(2.34*qkh) * pow(1.0+3.89*qkh+
         pow(16.1*qkh, 2.0)+pow(5.46*qkh, 3.0)+pow(6.71*qkh, 4.0), -0.25);
   return(t_f);
}

real CosmoClass::bbks(real k){
   real qkh, t_f;
      
   if (k == 0.0) return(0.0);
   qkh = k/kh_tmp;
      /* NOTE: the following line has 0/0 for k=0.
   This was taken care of at the beginning of the routine. */
   t_f = log(1.0+2.34*qkh)/(2.34*qkh) * pow(1.0+3.89*qkh+
         pow(16.1*qkh, 2.0)+pow(5.46*qkh, 3.0)+pow(6.71*qkh, 4.0), -0.25);
   return(t_f);
}

real CosmoClass::TransferFunction(real k){
   return((this->*tan_f)(k));
}

/* Sigma(R) routines */

real CosmoClass::sigma2(real k){
   real prim_ps, w_f, t_f, s2;
      
   prim_ps = Pk_norm*pow(k, n_s);
   w_f = 3.0*(sin(k*R_M) - k*R_M*cos(k*R_M))/pow(k*R_M, 3.0);
   t_f = TransferFunction(k);
   s2 = k*k * t_f*t_f * w_f*w_f * prim_ps;
   return(s2);
}

real CosmoClass::Sigma_r(real R, real AA){
   real sigma;
   const real k_min=0.0, k_max=last_k;
      
   R_M = R;
   Pk_norm = AA;
      
   sigma = 1.0/(2.0*pi*pi) * integrate(&CosmoClass::sigma2, k_min, k_max);
   sigma = sqrt(sigma);
   return(sigma);
}

#undef pi

/* Routines which calculate growth factor */

real CosmoClass::gf_func(real xx){
   real tmp1;
   tmp1 = pow(xx, 6.0/5.0);
   return(pow(1.0+gf_tmp_x*tmp1, -3.0/2.0));
}

real CosmoClass::LCDM_growth(real z){
   const real a=0.0, b=1.0;
   const real om_in=1.0/Omega_m;
   real g_int, gf;
      
   gf_tmp_x = (om_in-1.0)/pow(1.0+z,3.0);
   g_int = integrate(&CosmoClass::gf_func, a, b);
   gf = 1.0/(1.0+z) * pow((1.0+gf_tmp_x), 0.5) * g_int;
   return(gf);
}

real CosmoClass::GrowthFactor(real z){
   real growth, norm;
      
   if (Omega_m == 1.0){  // CDM
      growth = 1.0/(1.0+z);
   }		
   else if (Omega_m < 1.0 && w_de == -1.0){   // LCDM
      norm   = LCDM_growth(0.0);
      growth = LCDM_growth(z)/norm;
   }
   else if (Omega_m < 1.0 && w_de != -1.0){ // wCDM
      printf("GrowthFactor: wCDM not implemented yet!");
      abort();
   }
   else{
      printf("GrowthFactor: unknown cosmology!");
      abort();
   }
   return(growth);
}

/* Routines for calculating the derivative of the 
   growth factor */

real CosmoClass::gf_deriv(real xx){
   real tmp1;
   tmp1 = pow(xx, 6.0/5.0);
   return(tmp1*pow(1.0+gf_tmp_x*tmp1, -5.0/2.0));
}

real CosmoClass::GrowthDeriv(real z){
   const real a=0.0, b=1.0;
   const real om_in=1.0/Omega_m;
   real growth, ddot, term1, term2, term3, ddint;
   real norm, hzin, zdot;
      
   if (Omega_m == 1.0){  // CDM
      growth = 1.0/(1.0+z);
   }		
   else if (Omega_m < 1.0 && w_de == -1.0){   // LCDM
      gf_tmp_x = (om_in-1.0)/pow(1.0+z,3.0);
      norm = LCDM_growth(0.0);
      growth = LCDM_growth(z);
      term1 = growth/(1.0+z);
      term2 = -1.5/(gf_tmp_x+1.0)*gf_tmp_x/(1.0+z)*growth;
      ddint = integrate(&CosmoClass::gf_deriv, a, b);
      term3 = 4.5*gf_tmp_x/pow(1.0+z,2.0)*pow((1.0+gf_tmp_x), 0.5)*ddint;
      hzin  = pow((1.0+Omega_m*(pow(1.0+z,3.0)-1.0)), 0.5);
      zdot  = (1.0+z)*hzin;
      ddot  = (term1 + term2 + term3)*zdot/norm;
   }
   else if (Omega_m < 1.0 && w_de != -1.0){ // wCDM
      printf("GrowthDeriv: wCDM not implemented yet!");
      abort();
   }
   else{
      printf("GrowthDeriv: unknown cosmology!");
      abort();
   }
   return(ddot);
}

/* Numerical integration routines from Numerical Recipes 
   calling:
      x = integrate(&CosmoClass::func, a, b)

   where
      func -- is the (real) function whose integral is calculated, 
                                 and has to be member of the CosmoClass
      a    -- is the (real) lower boundary for integration
      b    -- is the (real) upper boundary for integration
*/

#define FUNC(x) ((this->*func)(x))
real CosmoClass::midpoint(real (CosmoClass::*func)(real), 
                          real a, real b, integer n){
                             real x,tnm,sum,del,ddel;
                             static real s;
                             integer it,j;
      
                             if (n == 1) {
                                return (s=(b-a)*FUNC(0.5*(a+b)));
                             } 
                             else {
                                it = 1;
                                for(j=1; j<n-1; ++j) it *= 3;
                                tnm=it;
                                del=(b-a)/(3.0*tnm);
                                ddel=del+del;
                                x=a+0.5*del;
                                sum=0.0;
                                for (j=1; j<=it; ++j){
                                   sum += FUNC(x);
                                   x += ddel;
                                   sum += FUNC(x);
                                   x += del;
                                }
                                s=(s+(b-a)*sum/tnm)/3.0;
                                return (s);
                             }
                          }
#undef FUNC

#define EPS 1.0e-4
#define JMAX 20
#define JMIN 5
                          real CosmoClass::integrate(real (CosmoClass::*func)(real), real a, real b){
                             real sol, old;
                             integer j;
      
                             sol = 0.0;
                             old = -1.0e26;
                             for (j=1; j<=JMAX; ++j) {
                                sol=midpoint(func, a, b, j);
                                if (j > JMIN){
                                   if (fabs(sol-old) < EPS*fabs(old) || 
                                       (sol==0.0 && old==0.0)) return(sol);
                                }
                                old = sol;
                             }
                             printf("integrate: no convergence!");
                             abort();
                          }
#undef EPS
#undef JMAX
#undef JMIN

/* Linear interpolation routine
   calling:
                          y = interpolate(xx, yy, n, x)

                          where
                          xx -- is the (real) x array of ordered data
                          yy -- is the (real) y array
                          n  -- size of above arrays
                          x  -- is the (real) point whose y value should be interpolated
*/

                          real CosmoClass::interpolate(real xx[], real yy[], unsigned long n, real x){
                             real y, dx, dy;
                             unsigned long jlo;
      
                             locate(xx, n, x, &jlo);
      // Linear interpolation:
                             dx = xx[jlo] - xx[jlo+1];
                             dy = yy[jlo] - yy[jlo+1];
                             y = dy/dx*(x-xx[jlo]) + yy[jlo];
      
                             return(y);
                          }

/* Routines for locating a value in an ordered table 
                          from Numerical Recipes */

                          void CosmoClass::locate(real xx[], unsigned long n, real x, unsigned long *j){
                             unsigned long ju,jm,jl;
                             int ascnd;
      
                             jl=0;
                             ju=n-1;
                             ascnd=(xx[n-1] >= xx[0]);
                             while (ju-jl > 1) {
                                jm=(ju+jl)/2;
                                if ((x >= xx[jm]) == ascnd)
                                   jl=jm;
                                else
                                   ju=jm;
                             }
                             if (x == xx[0]) *j=0;
                             else if(x == xx[n-1]) *j=n-2;
                             else *j=jl;
                             return;
                          }

                          void CosmoClass::hunt(real xx[], unsigned long n, real x, unsigned long *jlo){
                             unsigned long jm,jhi,inc;
                             int ascnd;
      
                             ascnd=(xx[n-1] >= xx[0]);
                             if (*jlo <= 0 || *jlo > n) {
                                *jlo=0;
                                jhi=n+1;
                             } else {
                                inc=1;
                                if (x >= xx[*jlo] == ascnd) {
                                   if (*jlo == n) return;
                                   jhi=(*jlo)+1;
                                   while (x >= xx[jhi] == ascnd) {
                                      *jlo=jhi;
                                      inc += inc;
                                      jhi=(*jlo)+inc;
                                      if (jhi > n) {
                                         jhi=n+1;
                                         break;
                                      }
                                   }
                                } else {
                                   if (*jlo == 1) {
                                      *jlo=0;
                                      return;
                                   }
                                   jhi=(*jlo)--;
                                   while (x < xx[*jlo] == ascnd) {
                                      jhi=(*jlo);
                                      inc <<= 1;
                                      if (inc >= jhi) {
                                         *jlo=0;
                                         break;
                                      }
                                      else *jlo=jhi-inc;
                                   }
                                }
                             }
                             while (jhi-(*jlo) != 1) {
                                jm=(jhi+(*jlo)) >> 1;
                                if (x >= xx[jm] == ascnd)
                                   *jlo=jm;
                                else
                                   jhi=jm;
                             }
                             if (x == xx[n]) *jlo=n-1;
                             if (x == xx[1]) *jlo=1;
                          }

#ifdef INITIALIZER_NAMESPACE
}
#endif

