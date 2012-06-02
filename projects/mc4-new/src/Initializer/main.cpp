/* Example driver for initializer */

#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "TypesAndDefs.h"
#include "InputParser.h"
#include "Initializer.h"
#include "Output.h"
#include "PerfMon.h"
#ifdef H5PART
  #include <hdf5.h>
  #include "H5Part.h"
  extern "C" {
    #include "H5Block.h"
    #include "H5BlockTypes.h"
  }
#endif

using namespace initializer;

Initializer init;

int main(int argc, char * const argv[]) {
  
   if(argc < 4) {
      fprintf(stderr,"USAGE: init <indatName> <tfName> <outBase>\n");
      exit(-1);
   }
   string indatName = argv[1];
   string tfName = argv[2];
   string outBase = argv[3];
   
   StartMonitor();
   InputParser par(indatName);
   
   real *pos_x, *pos_y, *pos_z, *vel_x, *vel_y, *vel_z;
   long Npart, i;
   int NumProcs, MyPE;
   integer *ID;
   Output out;
   //Initializer init;
   
   MPI::Init();
   
   MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &MyPE);
   
   int np = 0;
   par.getByName("np", np);
   Npart = (np/NumProcs)*np*np;
   
   real Omega_nu;
   par.getByName("Omega_nu", Omega_nu);
   if (Omega_nu > 0.0) {
      int neut_pairs = 0;
      par.getByName("nu_pairs", neut_pairs);
      Npart = (2*neut_pairs+1)*Npart;
   }
   
   pos_x = (real *)malloc(Npart*sizeof(real));
   pos_y = (real *)malloc(Npart*sizeof(real));
   pos_z = (real *)malloc(Npart*sizeof(real));
   vel_x = (real *)malloc(Npart*sizeof(real));
   vel_y = (real *)malloc(Npart*sizeof(real));
   vel_z = (real *)malloc(Npart*sizeof(real));
   
   init.init_particles(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, 
                       par, tfName.c_str());
   
   ID = (integer *)malloc(Npart*sizeof(integer)); 	 
   for(i=0; i<Npart; i++) { 	 
      ID[i] = Npart*MyPE + i; 	 
   }
   
   real rL = 0.0;
   par.getByName("box_size", rL);
   out.grid2phys(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, 	 
                 Npart, np, rL );
   int PrintFormat;
   par.getByName("PrintFormat", PrintFormat);
   
   OnClock("Output");
   if (MyPE == 0) {std::cout << std::endl << "Writing output................" 
                             << std::flush;}
   if (PrintFormat == 0) {  // ASCII, the same format as MC2
      out.write_hcosmo_ascii(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
                            Npart, outBase);
   }
   else if (PrintFormat == 1) { // Binary serial
      out.write_hcosmo_serial(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
                              ID, Npart, outBase);
   }
   else if (PrintFormat == 2) { // Binary parallel 
      out.write_hcosmo_parallel(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
                                ID, Npart, outBase);
   }
   else if (PrintFormat == 3) { // HDF5
#ifdef H5PART
      h5part_int64_t ng = bdata.ng();
      h5part_int64_t ng2d = bdata.ng2d();
      h5part_int64_t np = bdata.np();
      h5part_int64_t rL = bdata.rL();
      np *= np * np;
      out.write_hcosmo_h5(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
                          ID, Npart, outBase, ng, ng2d, np, rL, indatName);
#else
      if (MyPE == 0) {std::cout << "You must enable HDF5 in the Makefile!" 
                      << std::endl;}
#endif
   }
   else if (PrintFormat == 4) {  // Gadget format
      out.write_gadget(par, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
                       ID, Npart, outBase);
   }
   else {
		if (MyPE == 0)
   	   std::cout << "no output saved on disk...";
	}
   if (MyPE == 0) {std::cout << "done" << std::endl << std::flush;}
   OffClock("Output");
   
   free(pos_x);
   free(pos_y);
   free(pos_z);
   free(vel_x);
   free(vel_y);
   free(vel_z);
   free(ID);
   
   MPI_Barrier(MPI_COMM_WORLD);
   if (MyPE == 0) PrintClockSummary(stdout);
   MPI_Finalize();
   return(0);
}
