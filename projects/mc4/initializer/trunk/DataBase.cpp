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
#include "Basedata.h"

#ifdef INITIALIZER_NAMESPACE
namespace initializer {
#endif

void GlobalStuff::GetParameters(Basedata& bdata, ParallelTools& Parallel){
   int MyPE, MasterPE;

   MyPE = Parallel.GetMyPE();
   MasterPE = Parallel.GetMasterPE();
      
      // Get simulation parameters:
   ngrid = bdata.np();
   box_size = bdata.rL();
   dim = 1;
   seed = bdata.iseed();
   z_in = bdata.zin();
      
      // Get cosmology parameters:
   Omega_m = bdata.omegatot();
   Omega_bar = bdata.deut()/bdata.hubble()/bdata.hubble();
   Hubble = bdata.hubble();
   Sigma_8 = bdata.ss8();
   n_s = bdata.ns();
   w_de = -1.0;
   TFFlag = bdata.trans();
      
      // Get code parameters:
   PrintFlag = 0;
      
      // Make some sanity checks:
      // Errors:
   if (Omega_m > 1.0)
      Parallel.ParallelError("Cannot initialize closed universe!", "stop");
   if (Omega_bar > Omega_m)
      Parallel.ParallelError("More baryons than total matter content!", "stop");
   if (dim < 1 || dim > 3)
      Parallel.ParallelError("1, 2 or 3D decomposition is possible!", "stop");
      
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

#ifdef INITIALIZER_NAMESPACE
}
#endif
