/*
   Initializer:
   DataBase.cpp

      Has only 
         GetParameters(Basedata& bdata, ParallelTools& Parallel)
      routine which reads in simulation and cosmology 
      parameters from the main code's Basedata and stors them into 
      DataBase class (defined in DataBase.h).

                     Zarija Lukic, February 2009
                           zarija@lanl.gov
*/

#include <iostream>
#include <fstream>
#include "TypesAndDefs.h"
#include "DataBase.h"
#include "Parallelization.h"
//#include "Basedata.h"
#include "InputParser.h"

#ifdef USENAMESPACE
namespace initializer {
#endif

void GlobalStuff::GetParameters(InputParser& par, ParallelTools& Parallel){
   int MyPE, MasterPE;

   MyPE = Parallel.GetMyPE();
   MasterPE = Parallel.GetMasterPE();
   
   int np;
   if(!par.getByName("np", np) && MyPE==MasterPE)  // Grid size is equal to np^3
       std::cout << "Error: np not found!" << std::endl;
   ngrid = np;
   
   if(!par.getByName("box_size", box_size) && MyPE==MasterPE)
       std::cout << "Error: box_size not found!" << std::endl;
   
   int decomposition_type=1;
   /*
   if(!par.getByName("dim", decomposition_type)){
      if (MyPE==MasterPE){
         std::cout << "Warning: dim not found, assuming slab (1D) decomposition";
         std::cout << std::endl;
      }
   }
   */
   dim = decomposition_type;

   int sseed;
   if(!par.getByName("seed", sseed) && MyPE==MasterPE)
       std::cout << "Error: seed not found!" << std::endl;
   seed = sseed;

   if(!par.getByName("z_in", z_in) && MyPE==MasterPE)
       std::cout << "Error: z_in not found!" << std::endl;
      
   // Get cosmology parameters:
   if(!par.getByName("hubble", Hubble) && MyPE==MasterPE)
       std::cout << "Error: hubble not found!" << std::endl;

   if(!par.getByName("Omega_m", Omega_m) && MyPE==MasterPE)
       std::cout << "Error: Omega_m not found!" << std::endl;

   if(!par.getByName("Omega_bar", Omega_bar) && MyPE==MasterPE)
       std::cout << "Error: Omega_bar not found!" << std::endl;

   if(!par.getByName("Omega_nu", Omega_nu) && MyPE==MasterPE)
       std::cout << "Error: Omega_nu not found!" << std::endl;

   if(!par.getByName("Sigma_8", Sigma_8) && MyPE==MasterPE)
       std::cout << "Error: Sigma_8 not found!" << std::endl;
   
   if(!par.getByName("n_s", n_s) && MyPE==MasterPE)
       std::cout << "Error: n_s not found!" << std::endl;
   
   if(!par.getByName("w_de", w_de) && MyPE==MasterPE)
       std::cout << "Error: w_de not found!" << std::endl;
   
   if(!par.getByName("TFFlag", TFFlag) && MyPE==MasterPE)
       std::cout << "Error: TFFlag not found!" << std::endl;

   // Get neutrino parameters:
   if(!par.getByName("N_nu", N_nu) && MyPE==MasterPE)
       std::cout << "Error: N_nu not found!" << std::endl;
   
   if(!par.getByName("nu_pairs", nu_pairs) && MyPE==MasterPE)
       std::cout << "Error: nu_pairs not found!" << std::endl;
   
   // Get non-gaussianity parameters:
   if(!par.getByName("f_NL", f_NL)){
      if (MyPE==MasterPE)
         std::cout << "Warning: f_NL not found, assuming it is zero!" << std::endl;
      f_NL = 0.0;
   }
   
   // Get code parameters:
   if(!par.getByName("PrintFormat", PrintFormat) && MyPE==MasterPE)
      std::cout << "Error: PrintFormat not found!" << std::endl;

   
   // Make some sanity checks:
   // Errors:
   if (Omega_m > 1.0)
      Parallel.ParallelError("Cannot initialize closed universe!", "stop");
   if (Omega_bar > Omega_m)
      Parallel.ParallelError("More baryons than total matter content!", "stop");
   if (dim < 1 || dim > 3)
      Parallel.ParallelError("1, 2 or 3D decomposition is possible!", "stop");
   if (Omega_nu > 0 && nu_pairs < 1)
      Parallel.ParallelError("Omega_nu > 0, yet nu_pairs < 1!", "stop");
      
   // Warnings:
   if (z_in > 1000.0)
      Parallel.ParallelError
            ("Starting before CMB epoch. Make sure you want that.", " ");
   if (Sigma_8 > 1.0)
      Parallel.ParallelError("Sigma(8) > 1. Make sure you want that.", " ");
   if (Hubble > 1.0)
      Parallel.ParallelError("H > 100 km/s/Mpc. Make sure you want that.", " ");


   if (MyPE == MasterPE) 
      std::cout << "Done loading parameters." << std::endl;

   return;
}

#ifdef USENAMESPACE
}
#endif
