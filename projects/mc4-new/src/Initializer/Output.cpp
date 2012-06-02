/*
 *  Output.cpp
 *  Initializer
 *
 *  Created by Zarija on 1/27/10.
 *
 */

#include <mpi.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "TypesAndDefs.h"
#include "Output.h"
#include "InputParser.h"

#ifdef H5PART
#include <hdf5.h>
#include "H5Part.h"
extern "C" {
#include "H5Block.h"
#include "H5BlockTypes.h"
}
#endif


#ifdef USENAMESPACE
namespace initializer {
#endif
using namespace std;

void Output::grid2phys(real *pos_x, real *pos_y, real *pos_z, 	 
               real *vel_x, real *vel_y, real *vel_z, 	 
               int Npart, int np, float rL) { 	 
   long i; 	 
   
   float grid2phys_pos = 1.0*rL/np; 	 
   float grid2phys_vel = 100.0*rL/np; 	 
   
   for(i=0; i<Npart; i++) { 	 
      pos_x[i] *= grid2phys_pos; 	 
      pos_y[i] *= grid2phys_pos; 	 
      pos_z[i] *= grid2phys_pos; 	 
      vel_x[i] *= grid2phys_vel; 	 
      vel_y[i] *= grid2phys_vel; 	 
      vel_z[i] *= grid2phys_vel; 	 
   } 	 
   
   return; 	 
} 	 

   
void Output::write_hcosmo_serial(real *pos_x, real *pos_y, real *pos_z, 	 
                  real *vel_x, real *vel_y, real *vel_z, 	 
                  integer *id, int Npart, string outBase) { 	 
   FILE *outFile;
   ostringstream outName;
   long i;
   int MyPE, NumPEs, MasterPE, proc;
   MPI_Status status;
   
   MPI_Comm_size(MPI_COMM_WORLD, &NumPEs);
   MPI_Comm_rank(MPI_COMM_WORLD, &MyPE);
   MasterPE = 0;
   
   outName << outBase << ".bin";
   if (MyPE == MasterPE){
      outFile = fopen(outName.str().c_str(), "wb");
      for(i=0; i<Npart; i++) { 
         fwrite(&pos_x[i], sizeof(real), 1, outFile); 	 
         fwrite(&vel_x[i], sizeof(real), 1, outFile); 	 
         fwrite(&pos_y[i], sizeof(real), 1, outFile); 	 
         fwrite(&vel_y[i], sizeof(real), 1, outFile); 	 
         fwrite(&pos_z[i], sizeof(real), 1, outFile); 	 
         fwrite(&vel_z[i], sizeof(real), 1, outFile); 	 
         fwrite(&id[i], sizeof(integer), 1, outFile); 	 
      }
      for (proc=0; proc<NumPEs; ++proc){  // Get data from other processors:
         if (proc == MasterPE) continue;
         MPI_Recv(pos_x, Npart, MY_MPI_REAL, proc, 101, MPI_COMM_WORLD, &status);
         MPI_Recv(pos_y, Npart, MY_MPI_REAL, proc, 102, MPI_COMM_WORLD, &status);
         MPI_Recv(pos_z, Npart, MY_MPI_REAL, proc, 103, MPI_COMM_WORLD, &status);
         MPI_Recv(vel_x, Npart, MY_MPI_REAL, proc, 104, MPI_COMM_WORLD, &status);
         MPI_Recv(vel_y, Npart, MY_MPI_REAL, proc, 105, MPI_COMM_WORLD, &status);
         MPI_Recv(vel_z, Npart, MY_MPI_REAL, proc, 106, MPI_COMM_WORLD, &status);
         MPI_Recv(id, Npart, MY_MPI_INTEGER, proc, 107, MPI_COMM_WORLD, &status);
         for (i=0; i<Npart; ++i) {
            fwrite(&pos_x[i], sizeof(real), 1, outFile); 	 
            fwrite(&vel_x[i], sizeof(real), 1, outFile); 	 
            fwrite(&pos_y[i], sizeof(real), 1, outFile); 	 
            fwrite(&vel_y[i], sizeof(real), 1, outFile); 	 
            fwrite(&pos_z[i], sizeof(real), 1, outFile); 	 
            fwrite(&vel_z[i], sizeof(real), 1, outFile); 	 
            fwrite(&id[i], sizeof(integer), 1, outFile);
         }
      }
      fclose(outFile);
   }
   else {  // Send data to master processor:
      MPI_Send(pos_x, Npart, MY_MPI_REAL, MasterPE, 101, MPI_COMM_WORLD);
      MPI_Send(pos_y, Npart, MY_MPI_REAL, MasterPE, 102, MPI_COMM_WORLD);
      MPI_Send(pos_z, Npart, MY_MPI_REAL, MasterPE, 103, MPI_COMM_WORLD);
      MPI_Send(vel_x, Npart, MY_MPI_REAL, MasterPE, 104, MPI_COMM_WORLD);
      MPI_Send(vel_y, Npart, MY_MPI_REAL, MasterPE, 105, MPI_COMM_WORLD);
      MPI_Send(vel_z, Npart, MY_MPI_REAL, MasterPE, 106, MPI_COMM_WORLD);
      MPI_Send(id, Npart, MY_MPI_INTEGER, MasterPE, 107, MPI_COMM_WORLD);
   }
   
   return; 	 
}
   
   
void Output::write_hcosmo_parallel(real *pos_x, real *pos_y, real *pos_z, 	 
                                   real *vel_x, real *vel_y, real *vel_z, 	 
                                   integer *id, int Npart, string outBase) { 	 
   FILE *outFile;
   ostringstream outName;
   long i;
   int MyPE;
   
   MPI_Comm_rank(MPI_COMM_WORLD, &MyPE);
   
   outName << outBase << ".bin." << MyPE;
   outFile = fopen(outName.str().c_str(), "wb"); 	 
   for(i=0; i<Npart; i++) { 
      fwrite(&pos_x[i], sizeof(real), 1, outFile); 	 
      fwrite(&vel_x[i], sizeof(real), 1, outFile); 	 
      fwrite(&pos_y[i], sizeof(real), 1, outFile); 	 
      fwrite(&vel_y[i], sizeof(real), 1, outFile); 	 
      fwrite(&pos_z[i], sizeof(real), 1, outFile); 	 
      fwrite(&vel_z[i], sizeof(real), 1, outFile); 	 
      fwrite(&id[i], sizeof(integer), 1, outFile); 	 
   }
   fclose(outFile); 	 
   
   return; 	 
}


#ifdef H5PART
void Output::write_hcosmo_h5(real *pos_x, real *pos_y, real *pos_z,
                     real *vel_x, real *vel_y, real *vel_z,
                     integer *id, int Npart, string outBase, 
                     h5part_int64_t ng, h5part_int64_t ng2d, h5part_int64_t np, 
                     h5part_int64_t rL, string indatName) {
   long i;
   int status=0;
   
   int H5call = 0;
   H5PartFile *H5file;
   ostringstream fn;
   fn << outBase << ".h5";
   
#ifdef PARALLEL_IO
   H5file = H5PartOpenFileParallel(fn.str().c_str(),H5PART_APPEND,MPI_COMM_WORLD);
#else
   H5file = H5PartOpenFile(fn.str().c_str(),H5PART_APPEND);
#endif
   
   if(!H5file) {
      cout << "h5 file open failed: exiting!" << endl;
      exit(0);
   }
   
   H5PartWriteFileAttribString(H5file,"tUnit","s");
   H5PartWriteFileAttribString(H5file,"xUnit","Mpc/h");
   H5PartWriteFileAttribString(H5file,"yUnit","Mpc/h");
   H5PartWriteFileAttribString(H5file,"zUnit","Mpc/h");
   H5PartWriteFileAttribString(H5file,"pxUnit","km/s");
   H5PartWriteFileAttribString(H5file,"pyUnit","km/s");
   H5PartWriteFileAttribString(H5file,"pzUnit","km/s");
   H5PartWriteFileAttribString(H5file,"idUnit","1");
   
   H5PartWriteFileAttribString(H5file,"TIMEUnit","s");
   
   H5PartWriteFileAttrib(H5file, "ng", H5PART_INT64, &ng, 1);
   H5PartWriteFileAttrib(H5file, "ng2d", H5PART_INT64, &ng2d, 1);
   H5PartWriteFileAttrib(H5file, "np", H5PART_INT64, &np, 1);
   H5PartWriteFileAttrib(H5file, "rL", H5PART_INT64, &rL, 1);
   H5PartWriteFileAttribString(H5file, "input filename", indatName.c_str());
   
   void *varray = malloc(Npart*sizeof(double));
   double *farray = (double*)varray;
   h5part_int64_t *larray = (h5part_int64_t *)varray;
   
   /// Set current record/time step.
   H5PartSetStep(H5file, 0);
   H5PartSetNumParticles(H5file, Npart);
   
   for(size_t i=0; i<Npart; i++)
      farray[i] =  pos_x[i];
   H5PartWriteDataFloat64(H5file,"x",farray);
   for(size_t i=0; i<Npart; i++)
      farray[i] =  pos_y[i];
   H5PartWriteDataFloat64(H5file,"y",farray);
   for(size_t i=0; i<Npart; i++)
      farray[i] =  pos_z[i];
   H5PartWriteDataFloat64(H5file,"z",farray);
   
   for(size_t i=0; i<Npart; i++)
      farray[i] =  vel_x[i];
   H5PartWriteDataFloat64(H5file,"px",farray);
   for(size_t i=0; i<Npart; i++)
      farray[i] =  vel_y[i];
   H5PartWriteDataFloat64(H5file,"py",farray);
   for(size_t i=0; i<Npart; i++)
      farray[i] =  vel_z[i];
   H5PartWriteDataFloat64(H5file,"pz",farray);
   
   // Write particle id numbers.
   for (size_t i = 0; i < Npart; i++)
      larray[i] =  id[i];
   H5PartWriteDataInt64(H5file,"id",larray);
   
   H5Fflush(H5file->file,H5F_SCOPE_GLOBAL);
   
   if(varray)
      free(varray);
   
   H5PartCloseFile(H5file);
   
   return;
}
#endif


void Output::write_hcosmo_ascii(real *pos_x, real *pos_y, real *pos_z,
                       real *vel_x, real *vel_y, real *vel_z,
                       int Npart, string outBase) {
   long i;
   int MyPE, NumPEs, MasterPE, proc;
   MPI_Status status;
   ostringstream outName;
   ofstream of;
   
   MPI_Comm_size(MPI_COMM_WORLD, &NumPEs);
   MPI_Comm_rank(MPI_COMM_WORLD, &MyPE);
   MasterPE = 0;
   
   outName << outBase << ".txt";   
   
   if (MyPE == MasterPE){
      of.open(outName.str().c_str(), ios::out);
      for (i=0; i<Npart; ++i) {
         of  << pos_x[i] << ' ';
         of <<  vel_x[i] << ' ';
         of <<  pos_y[i] << ' ';
         of <<  vel_y[i] << ' ';
         of <<  pos_z[i] << ' ';
         of <<  vel_z[i] << endl;
      }
      for (proc=0; proc<NumPEs; ++proc){  // Get data from other processors:
         if (proc == MasterPE) continue;
         MPI_Recv(pos_x, Npart, MY_MPI_REAL, proc, 101, MPI_COMM_WORLD, &status);
         MPI_Recv(pos_y, Npart, MY_MPI_REAL, proc, 102, MPI_COMM_WORLD, &status);
         MPI_Recv(pos_z, Npart, MY_MPI_REAL, proc, 103, MPI_COMM_WORLD, &status);
         MPI_Recv(vel_x, Npart, MY_MPI_REAL, proc, 104, MPI_COMM_WORLD, &status);
         MPI_Recv(vel_y, Npart, MY_MPI_REAL, proc, 105, MPI_COMM_WORLD, &status);
         MPI_Recv(vel_z, Npart, MY_MPI_REAL, proc, 106, MPI_COMM_WORLD, &status);
         for (i=0; i<Npart; ++i) {
            of  << pos_x[i] << ' ';
            of <<  vel_x[i] << ' ';
            of <<  pos_y[i] << ' ';
            of <<  vel_y[i] << ' ';
            of <<  pos_z[i] << ' ';
            of <<  vel_z[i] << endl;
         }
      }
      of.close();
   }
   else {  // Send data to master processor:
      MPI_Send(pos_x, Npart, MY_MPI_REAL, MasterPE, 101, MPI_COMM_WORLD);
      MPI_Send(pos_y, Npart, MY_MPI_REAL, MasterPE, 102, MPI_COMM_WORLD);
      MPI_Send(pos_z, Npart, MY_MPI_REAL, MasterPE, 103, MPI_COMM_WORLD);
      MPI_Send(vel_x, Npart, MY_MPI_REAL, MasterPE, 104, MPI_COMM_WORLD);
      MPI_Send(vel_y, Npart, MY_MPI_REAL, MasterPE, 105, MPI_COMM_WORLD);
      MPI_Send(vel_z, Npart, MY_MPI_REAL, MasterPE, 106, MPI_COMM_WORLD);
   }
   
   return;
}


void Output::write_gadget(InputParser& par, real *pos_x, real *pos_y, real *pos_z,
                           real *vel_x, real *vel_y, real *vel_z, integer *id, 
                           int Npart, string outBase) {
   FILE *outFile;
   ostringstream outName;
   long i;
   int np;
   real rL, Om, h, redshift;
   gadget_header ghead;
   int MyPE, NumPEs, MasterPE, proc;
   MPI_Status status;
   float x, y, z, vx, vy, vz;
   unsigned int gID;
   
   MPI_Comm_size(MPI_COMM_WORLD, &NumPEs);
   MPI_Comm_rank(MPI_COMM_WORLD, &MyPE);
   MasterPE = 0;
   outName << outBase << ".gadget";
   
   // Gadget header -- dark matter only:
   ghead.flag_sfr = 0;
   ghead.flag_feedback = 0;
   ghead.flag_cooling = 0;
   ghead.num_files = 1;
   ghead.flag_stellarage = 0;
   ghead.flag_metals = 0;
   for (i=0; i<6; ++i) {ghead.npartTotalHighWord[i]=0;}
   ghead.flag_entropy_instead_u = 0;
   
   par.getByName("box_size", rL);
   par.getByName("Omega_m", Om);
   par.getByName("hubble", h);
   par.getByName("z_in", redshift);
   ghead.BoxSize = (double)rL; // In Mpc/h
   ghead.Omega0  = (double)Om;
   ghead.OmegaLambda = (double)(1.0-Om);  // Flat universe
   ghead.HubbleParam = (double)h;
   ghead.redshift = (double)redshift;
   ghead.time = (double)(1.0/(1.0+redshift));
   
   for (i=0; i<6; ++i) {
      ghead.npart[i] = 0;
      ghead.npartTotal[i] = 0;
      ghead.mass[i] = 0.0;
   }
   par.getByName("np", np);
   ghead.npart[1] = np*np*np;
   ghead.npartTotal[1] = np*np*np;
   ghead.mass[1] = (double)27.75197*Om*pow(rL/np, 3); // In 10^10 M_sun/h
   
   // Write it:
   int blksize;
#define SKIP {fwrite(&blksize, sizeof(int), 1, outFile);}
   blksize = sizeof(ghead);
   if (MyPE == MasterPE){
      outFile = fopen(outName.str().c_str(), "wb");
      SKIP
      fwrite(&ghead, sizeof(ghead), 1, outFile);
      SKIP
      blksize = np*np*np*3*sizeof(float);
      SKIP
      for(i=0; i<Npart; i++) { 
         x = (float)pos_x[i];  // In Mpc/h
         y = (float)pos_y[i];
         z = (float)pos_z[i];
         fwrite(&x, sizeof(float), 1, outFile); 	 
         fwrite(&y, sizeof(float), 1, outFile); 	 
         fwrite(&z, sizeof(float), 1, outFile);
      }
      for (proc=0; proc<NumPEs; ++proc){  // Get data from other processors:
         if (proc == MasterPE) continue;
         MPI_Recv(pos_x, Npart, MY_MPI_REAL, proc, 101, MPI_COMM_WORLD, &status);
         MPI_Recv(pos_y, Npart, MY_MPI_REAL, proc, 102, MPI_COMM_WORLD, &status);
         MPI_Recv(pos_z, Npart, MY_MPI_REAL, proc, 103, MPI_COMM_WORLD, &status);
         for (i=0; i<Npart; ++i) {
            x = (float)pos_x[i];  // In Mpc/h
            y = (float)pos_y[i];
            z = (float)pos_z[i];
            fwrite(&x, sizeof(float), 1, outFile); 	 
            fwrite(&y, sizeof(float), 1, outFile); 	 
            fwrite(&z, sizeof(float), 1, outFile);
         }
      }
      SKIP
   }
   else {  // Send data to master processor:
      MPI_Send(pos_x, Npart, MY_MPI_REAL, MasterPE, 101, MPI_COMM_WORLD);
      MPI_Send(pos_y, Npart, MY_MPI_REAL, MasterPE, 102, MPI_COMM_WORLD);
      MPI_Send(pos_z, Npart, MY_MPI_REAL, MasterPE, 103, MPI_COMM_WORLD);
   }
   if (MyPE == MasterPE){
      blksize = np*np*np*3*sizeof(float);
      SKIP
      for(i=0; i<Npart; i++) {
         vx = (float)vel_x[i];    // In km/s
         vy = (float)vel_y[i];
         vz = (float)vel_z[i];
         fwrite(&vx, sizeof(float), 1, outFile); 	 
         fwrite(&vy, sizeof(float), 1, outFile); 	 
         fwrite(&vz, sizeof(float), 1, outFile);
      }
      for (proc=0; proc<NumPEs; ++proc){  // Get data from other processors:
         if (proc == MasterPE) continue;
         MPI_Recv(vel_x, Npart, MY_MPI_REAL, proc, 104, MPI_COMM_WORLD, &status);
         MPI_Recv(vel_y, Npart, MY_MPI_REAL, proc, 105, MPI_COMM_WORLD, &status);
         MPI_Recv(vel_z, Npart, MY_MPI_REAL, proc, 106, MPI_COMM_WORLD, &status);
         for (i=0; i<Npart; ++i) {
            vx = (float)vel_x[i];    // In km/s
            vy = (float)vel_y[i];
            vz = (float)vel_z[i];
            fwrite(&vx, sizeof(float), 1, outFile); 	 
            fwrite(&vy, sizeof(float), 1, outFile); 	 
            fwrite(&vz, sizeof(float), 1, outFile);
         }
      }
      SKIP
   }
   else {  // Send data to master processor:
      MPI_Send(vel_x, Npart, MY_MPI_REAL, MasterPE, 104, MPI_COMM_WORLD);
      MPI_Send(vel_y, Npart, MY_MPI_REAL, MasterPE, 105, MPI_COMM_WORLD);
      MPI_Send(vel_z, Npart, MY_MPI_REAL, MasterPE, 106, MPI_COMM_WORLD);
   }
   if (MyPE == MasterPE){
      blksize = np*np*np*sizeof(int);
      SKIP
      for(i=0; i<Npart; i++) {
         gID = (int)id[i];
         fwrite(&gID, sizeof(int), 1, outFile);
      }
      for (proc=0; proc<NumPEs; ++proc){  // Get data from other processors:
         if (proc == MasterPE) continue;
         MPI_Recv(id, Npart, MY_MPI_INTEGER, proc, 107, MPI_COMM_WORLD, &status);
         for (i=0; i<Npart; ++i) {
            gID = (int)id[i];
            fwrite(&gID, sizeof(int), 1, outFile);
         }
      }
      SKIP
      fclose(outFile);
   }
   else {  // Send data to master processor:
      MPI_Send(id, Npart, MY_MPI_INTEGER, MasterPE, 107, MPI_COMM_WORLD);
   }
   
   return;
}
   
#ifdef USENAMESPACE
}
#endif
