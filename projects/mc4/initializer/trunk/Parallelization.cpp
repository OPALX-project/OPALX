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
               dim = 3 is cube decomposition
         and ngrid is number of grid zones per dimension. 
         Routine sets indicies for the arrays, but also determines 
         which processor hold mirror-symmetric data. This is needed 
         for enforcing reality condition of the density field.

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
#include "TypesAndDefs.h"
#include "Parallelization.h"

#ifdef INITIALIZER_NAMESPACE
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
   if (stop == "stop" or stop == "Stop" or stop == "STOP"){
      std::cout << std::endl;
      MPI_Abort(comm, rc);
   }
   return;
}


/* Handling errors which *SOME* processes encounter */
void ParallelTools::AbortMPI(const char *message){
   int rc;
   std::cout << message << std::endl << std::endl << std::flush;
   MPI_Abort(comm, rc);
   return;
}


void ParallelTools::InitMPI(double *w_time, MPI_Comm comm){
   char **argv;
   int argc, rc;
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
   int MyProc, NumProcs, MasterProc, MirrorProc, MirrorFlag;
   integer ndiv;
   integer nx1, nx2, ny1, ny2, nz1, nz2;
   integer x_pos, y_pos, z_pos;
      
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
         SetMirrorPE(0);
         SetMirrorFlag(0);
         std::cout << "...............done" << std::endl << std::flush;
         return;
      }
      if (MyProc == MasterProc){
	std::cout << std::endl;
	std::cout << "Initializer will use " << NumProcs << " processors." << std::endl;
	std::cout << ngrid << "^3 grid" << std::endl;
      }
      switch (dim) {
         case 1:
            ndiv = NumProcs;
            if (ngrid%ndiv != 0)
               {ParallelError("Ngrid not divisible by the number of lunched processes!", "stop");}

            if ((float)(ndiv/2) != ndiv/2.0) {
               ParallelError("Code has to run on even number processors for 1D decomposition!", "no stop");
               ndiv = ndiv - 1;
            }
            NumProcs = ndiv;
            if (MyProc == MasterProc){
               std::cout << std::endl;
               std::cout << "Initializer will use " << NumProcs << " processors." << std::endl;
               std::cout << ngrid << "^3 grid" << std::endl;
               std::cout << "Decomposing into slabs" << std::flush;
            }
            SetNumPEs(NumProcs);
            nx1 = MyProc*ngrid/ndiv;
            nx2 = (MyProc + 1)*ngrid/ndiv - 1;
            ny1 = 0;
            ny2 = ngrid - 1;
            nz1 = 0;
            nz2 = ngrid - 1;
            break;
         case 2:
            ndiv = (integer)sqrt((double)NumProcs);
            if ((float)(ndiv/2) != ndiv/2.0) {
               ParallelError("Code has to run on (even number)^2 processors for 2D decomposition!", "no stop");
               ndiv = ndiv - 1;
            }
            NumProcs = (integer)pow(ndiv, 2.0);
            if (NumProcs < 4)
               ParallelError("DomainDecompose: need at least 4 processors for 2D decomposition.", "stop");
            if (MyProc == MasterProc){
               std::cout << std::endl;
               std::cout << "Initializer will use " << NumProcs << " processors." << std::endl;
               std::cout << ngrid << "^3 grid" << std::endl;
               std::cout << "Decomposing into pencils" << std::flush;
            }
            SetNumPEs(NumProcs);
            y_pos = MyProc/ndiv;
            x_pos = MyProc - y_pos*ndiv;
            nx1 = x_pos*ngrid/ndiv;
            nx2 = (x_pos+1)*ngrid/ndiv - 1;
            ny1 = y_pos*ngrid/ndiv;
            ny2 = (y_pos+1)*ngrid/ndiv - 1;
            nz1 = 0;
            nz2 = ngrid - 1;
            break;
         case 3:
            ndiv = (integer)pow(NumProcs, 0.333334);
            if ((float)(ndiv/2) != ndiv/2.0) {
               ParallelError("Code has to run on (even number)^3 processors for 3D decompostion!", "no stop");
               ndiv = ndiv - 1;
            }
            NumProcs = (integer)pow(ndiv, 3.0);
            if (NumProcs < 8)
               ParallelError("DomainDecompose: need at least 8 processors for 3D decomposition.", "stop");
            if (MyProc == MasterProc){
               std::cout << std::endl;
               std::cout << "Initializer will use " << NumProcs << " processors." << std::endl;
               std::cout << ngrid << "^3 grid" << std::endl;
               std::cout << "Decomposing into cubes" << std::flush;
            }
            SetNumPEs(NumProcs);
            z_pos = MyProc/(ndiv*ndiv);
            y_pos = (MyProc - z_pos*ndiv*ndiv)/ndiv;
            x_pos = MyProc - z_pos*ndiv*ndiv - y_pos*ndiv;
            nx1 = x_pos*ngrid/ndiv;
            nx2 = (x_pos+1)*ngrid/ndiv - 1;
            ny1 = y_pos*ngrid/ndiv;
            ny2 = (y_pos+1)*ngrid/ndiv - 1;
            nz1 = z_pos*ngrid/ndiv;
            nz2 = (z_pos+1)*ngrid/ndiv - 1;
            break;
         default:
            ParallelError("DomainDecompose: dim has to be 1, 2 or 3!", "stop");
            break;
      }
      /* Set domain */
      SetNx1(nx1);
      SetNx2(nx2);
      SetNy1(ny1);
      SetNy2(ny2);
      SetNz1(nz1);
      SetNz2(nz2);
      
      /* Set mirror symmetric processor data */
      if (MyProc > NumProcs-1) {MirrorProc = 0;}
      else {MirrorProc = (integer)pow(ndiv, dim) - 1 - MyProc;}
      if (MyProc < NumProcs/2) {
         MirrorFlag = 0;
      }
      else {
         MirrorFlag = 1;
      }
      SetMirrorPE(MirrorProc);
      SetMirrorFlag(MirrorFlag);
      
      if (MyProc == MasterProc)
         std::cout << "...............done" << std::endl << std::flush;
      
      return;
}

#ifdef INITIALIZER_NAMESPACE
}
#endif
