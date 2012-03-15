/*
   Initializer:
   Cosmology.h

      Defines CosmoStuff class. This will be used 
      for calculating many quantities -- transfer function, 
      mass variance on certain scale, growth factor. 
      See Cosmology.cpp for more explanations and 
      actual implementation. 

                        Zarija Lukic, February 2009
                              zarija@lanl.gov
*/

#ifndef Cosmology_Header_Included
#define Cosmology_Header_Included

#include "TypesAndDefs.h"

#ifdef INITIALIZER_NAMESPACE
namespace initializer {
#endif

class GlobalStuff;

class CosmoClass{
      
   public :
      void SetParameters(GlobalStuff& DataBase, const char *tfName);
      real Sigma_r(real R, real norm);
      real GrowthFactor(real z);
      real GrowthDeriv(real z);	
      real TransferFunction(real k);
      ~CosmoClass();
      
   private :
      real cobe_temp, tt;
      real Omega_m, Omega_bar, h, n_s, w_de;
      integer TFFlag;
      real *table_kk, *table_tf;
      unsigned long table_size;

      real kh_tmp, R_M, Pk_norm, gf_tmp_x, last_k;
      
      real (CosmoClass::*tan_f)(real k);
      real integrate(real (CosmoClass::*func)(real), real a, real b);
      real midpoint(real (CosmoClass::*func)(real), real a, real b, integer n);
      real interpolate(real xx[], real yy[], unsigned long n, real x);
      void locate(real xx[], unsigned long n, real x, unsigned long *j);
      void hunt(real xx[], unsigned long n, real x, unsigned long *jlo);
      real klypin_holtzmann(real k);
      real hu_sugiyama(real k);
      real peacock_dodds(real k);
      real bbks(real k);
      real cmbfast(real k);
      real LCDM_growth(real z);
      real gf_func(real xx);
      real gf_deriv(real xx);
      real sigma2(real k);

};

#ifdef INITIALIZER_NAMESPACE
}
#endif

#endif
