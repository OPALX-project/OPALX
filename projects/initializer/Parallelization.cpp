/*
   Initializer:
   Parallelization.cpp

      Contains routines which handle mpi parallelization (see
         Paralelization.h for the interface). Important ones
         are:
               InitMPI(double *w_time)
                  which should be called at the start of the program; it
         starts mpi, sets process number, master process number,
         and number of processes. It also returns wall clock time
         via w_time variable. Similar to that,
               FinalMPI(double *wtime)
         closes mpi session, and returns end time.

         The largest and most important routine is
               DomainDecompose(integer dim, integer ngrid)
         where dim defines dimension of decomposition:
               dim = 1 is slab decomposition (x axis is sliced)
               dim = 2 is rod decomposition (rods go along z axis)
               dim = 3 is brick decomposition
         and ngrid is number of grid zones per dimension.
         Routine sets indicies for the arrays, needed for the loops.

         For example, full array:
               rho(0:ngrid-1, 0:ngrid-1, 0:ngrid-1)
         is divided such that each processor holds:
               rho(nx1:nx2, ny1:ny2, nz1:nz2)

         Serial job (running on 1 processor) is also handled, and
                  dim parameters is obviously not important there.

                           Zarija Lukic, February 2009
                                 zarija@lanl.gov
*/

#include "mpi.h"
#include <iostream>
#include <math.h>

#include <string.h>

#include "TypesAndDefs.h"
#include "Parallelization.h"

#ifdef USENAMESPACE
namespace initializer {
#endif

/* Routine should be used for errors which *ALL* processes encounter */
void ParallelTools::ParallelError(const char *message, const char *stop){
   int MyProc, MasterProc, rc;
   MyProc = GetMyPE();
   MasterProc = GetMasterPE();
   if (MyProc == MasterProc)
      std::cout << message << std::endl << std::flush;
   MPI_Barrier(comm);
   if (strcmp(stop, "stop") == 0 || strcmp(stop, "Stop") == 0 ||
       strcmp(stop, "STOP") == 0){
      std::cout << std::endl;
      rc = 911;
      MPI_Abort(comm, rc);
   }
   return;
}


/* Handling errors which *SOME* processes encounter */
void ParallelTools::AbortMPI(const char *message){
   int rc;
   std::cout << message << std::endl << std::endl << std::flush;
   rc = 912;
   MPI_Abort(comm, rc);
   return;
}


void ParallelTools::InitMPI(double *w_time, MPI_Comm comm){
   //char **argv;
   //int argc, rc;
   int MasterProc, MyProc, NumProcs;
//	rc = MPI_Init(&argc, &argv);
//	if (rc != MPI_SUCCESS)
//		ParallelError("Could not start MPI.", "stop");
   *w_time = MPI_Wtime();
   MPI_Comm_size(comm, &NumProcs);
   MPI_Comm_rank(comm, &MyProc);
   MasterProc = 0;
   SetMyPE(MyProc);
   SetMasterPE(MasterProc);
   SetNumPEs(NumProcs);
   this->comm = comm;
   return;
}


void ParallelTools::FinalMPI(double *w_time){
   *w_time = MPI_Wtime();
//	MPI_Finalize();
   return;
}


void ParallelTools::DomainDecompose(integer dim, integer ngrid){
   int MyProc, NumProcs, MasterProc;
   integer nx1, nx2, ny1, ny2, nz1, nz2, ngx, ngy, ngz;
   //integer x_pos, y_pos, z_pos;
   int reorder=1;
   int dim_sizes[3], is_periodic[3], coords[3];
   MPI_Comm MPI_CART_COMM;

   MyProc = GetMyPE();
   NumProcs = GetNumPEs();
   MasterProc = GetMasterPE();
      if (NumProcs == 1) {         // Serial job:
         std::cout << std::endl;
         std::cout << "Running as a serial job." << std::endl;
         std::cout << ngrid << "^3 grid" << std::endl;
         std::cout << "No decomposition" << std::flush;
         SetNx1(0);
         SetNx2(ngrid-1);
         SetNy1(0);
         SetNy2(ngrid-1);
         SetNz1(0);
         SetNz2(ngrid-1);
         std::cout << "...............done" << std::endl << std::flush;
         return;
      }

      is_periodic[0] = 1;
      is_periodic[1] = 1;
      is_periodic[2] = 1;
      if (MyProc == MasterProc){
          std::cout << std::endl;
          std::cout << "Code runs using " << NumProcs << " processors,";
          std::cout << " on " << ngrid << " X " << ngrid << " X ";
          std::cout << ngrid << " domain." << std::endl << std::flush;
      }
      switch (dim) {
         case 1:
             if (MyProc == MasterProc)
                 std::cout << "Decomposing into slabs: " <<std::endl << std::flush;
             dim_sizes[0] = 0;
             dim_sizes[1] = dim_sizes[2] = 1;
            break;
         case 2:
            if (MyProc == MasterProc)
                std::cout << "Decomposing into pencils: " <<std::endl << std::flush;
            dim_sizes[0] = dim_sizes[1] = 0;
            dim_sizes[2] = 1;
            break;
         case 3:
            if (MyProc == MasterProc)
                std::cout << "Decomposing into (3D) bricks: " <<std::endl << std::flush;
            dim_sizes[0] = dim_sizes[1] = dim_sizes[2] = 0;
            break;
         default:
            ParallelError("DomainDecompose: dim has to be 1, 2 or 3!", "stop");
            break;
      }

      MPI_Dims_create(NumProcs, 3, dim_sizes);
      MPI_Cart_create(MPI_COMM_WORLD, 3, dim_sizes, is_periodic, reorder, &MPI_CART_COMM);
      ngx = ngrid/dim_sizes[0];
      if (ngrid % dim_sizes[0] != 0)
          ParallelError("DomainDecompose: Cannot evenly decompose in X direction", "stop");
      ngy = ngrid/dim_sizes[1];
      if (ngrid % dim_sizes[1] != 0)
          ParallelError("DomainDecompose: Cannot evenly decompose in Y direction", "stop");
      ngz = ngrid/dim_sizes[2];
      if (ngrid % dim_sizes[2] != 0)
          ParallelError("DomainDecompose: Cannot evenly decompose in Z direction", "stop");

      if (MyProc == MasterProc){
          std::cout << "   Number of points in x: " << ngx <<std::endl;
          std::cout << "   Number of points in y: " << ngy <<std::endl;
          std::cout << "   Number of points in z: " << ngz <<std::endl << std::flush;
      }
      MPI_Cart_get(MPI_CART_COMM, 3, dim_sizes, is_periodic, coords);
      nx1 = coords[0]*ngrid/dim_sizes[0];
      nx2 = (coords[0]+1)*ngrid/dim_sizes[0] - 1;
      ny1 = coords[1]*ngrid/dim_sizes[1];
      ny2 = (coords[1]+1)*ngrid/dim_sizes[1] - 1;
      nz1 = coords[2]*ngrid/dim_sizes[2];
      nz2 = (coords[2]+1)*ngrid/dim_sizes[2] - 1;

      /* Set domain */
      SetNx1(nx1);
      SetNx2(nx2);
      SetNy1(ny1);
      SetNy2(ny2);
      SetNz1(nz1);
      SetNz2(nz2);

      return;
}

#ifdef USENAMESPACE
}
#endif

