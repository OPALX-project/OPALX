/* Example driver for initializer */

#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include "Basedata.h"
#include "Initializer.h"

#include <hdf5.h>
#include "H5Part.h"
extern "C" {
    #include "H5Block.h"
    #include "H5BlockTypes.h"
}

#include <string>

using namespace std;

#ifdef INITIALIZER_NAMESPACE
using namespace initializer; 
#endif

string create_outName(string outBase, int rank) {
  ostringstream outName;
  if (rank > 0)
      outName << outBase << "." << rank;
  else
      outName << outBase << ".dat" ;
  return outName.str();
}



void grid2phys(float *pos_x, float *pos_y, float *pos_z,
	       float *vel_x, float *vel_y, float *vel_z,
	       int Npart, int np, float rL) {
  int i;

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

void write_hcosmo_h5(float *pos_x, float *pos_y, float *pos_z,
		  float *vel_x, float *vel_y, float *vel_z,
		  int *id, int Npart, string outBase, 
          h5part_int64_t ng, h5part_int64_t ng2d, h5part_int64_t np, 
          h5part_int64_t rL, string indatName) {
  int i;
  int status=0;

  int H5call = 0;
  H5PartFile *H5file;
  ostringstream fn;
  fn << outBase << ".h5";

  //  H5PartSetVerbosityLevel(4);

#ifdef PARALLEL_IO
  H5file = H5PartOpenFileParallel(fn.str().c_str(), H5PART_WRITE, MPI_COMM_WORLD);
#else
  H5file = H5PartOpenFile(fn.str().c_str(), H5PART_WRITE);
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

  /// Write particle id numbers.
  for (size_t i = 0; i < Npart; i++)
      larray[i] =  id[i];
  H5PartWriteDataInt64(H5file,"id",larray);

  H5Fflush(H5file->file,H5F_SCOPE_GLOBAL);

  //  if(varray)
  //  free(varray);
    
  h5part_int64_t cres = H5PartCloseFile(H5file);
  
  if (cres < 0)
    cout << "Error code close h5 " << cres  << endl;



  return;
}

void write_hcosmo(float *pos_x, float *pos_y, float *pos_z,
		  float *vel_x, float *vel_y, float *vel_z,
		  int *id, int Npart, const char *outName) {
  int i;
  int status=0;
  FILE *outFile;
  outFile = fopen(outName, "wb");
  for(i=0; i<Npart; i++) {
    fwrite(&pos_x[i], sizeof(float), 1, outFile);
    fwrite(&vel_x[i], sizeof(float), 1, outFile);
    fwrite(&pos_y[i], sizeof(float), 1, outFile);
    fwrite(&vel_y[i], sizeof(float), 1, outFile);
    fwrite(&pos_z[i], sizeof(float), 1, outFile);
    fwrite(&vel_z[i], sizeof(float), 1, outFile);
    fwrite(&status, sizeof(int), 1, outFile);
    fwrite(&id[i], sizeof(int), 1, outFile);
  }
  fclose(outFile);
  
  return;
}


void write_hcosmo_fort(float *pos_x, float *pos_y, float *pos_z,
		       float *vel_x, float *vel_y, float *vel_z,
		       int Npart, const char *outName) {
    int i;
    ofstream of;
    
    of.open(outName,ios::out);

    for (size_t i = 0 ; i < Npart; ++i) {
	of  << pos_x[i] << ' ';
	of <<  vel_x[i] << ' ';
	of <<  pos_y[i] << ' ';
	of <<  vel_y[i] << ' ';
	of <<  pos_z[i] << ' ';
	of <<  vel_z[i] << endl;
    }
    of.close();  
   
    return;
}


int main(int argc, char * const argv[]) {
  
  if(argc < 4) {
    fprintf(stderr,"USAGE: init <indatName> <tfName> <outBase>\n");
    exit(-1);
  }
  string indatName = argv[1];
  string tfName = argv[2];
  string outBase = argv[3];
  
  Basedata bdata(indatName.c_str());

  float *pos_x, *pos_y, *pos_z, *vel_x, *vel_y, *vel_z;
  long Npart;
  int NumProcs, rank, *id, i;
  
  MPI::Init();

  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Npart = (bdata.np()/NumProcs)*bdata.np()*bdata.np();

  pos_x = (float *)malloc(Npart*sizeof(float));
  pos_y = (float *)malloc(Npart*sizeof(float));
  pos_z = (float *)malloc(Npart*sizeof(float));
  vel_x = (float *)malloc(Npart*sizeof(float));
  vel_y = (float *)malloc(Npart*sizeof(float));
  vel_z = (float *)malloc(Npart*sizeof(float));

  init_particles(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, 
		 bdata, tfName.c_str(), MPI_COMM_WORLD);

  if (bdata.pbchoice() == 2)
      grid2phys(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
		Npart, bdata.np(), bdata.rL() );

  id = (int *)malloc(Npart*sizeof(int));
  for(i=0; i<Npart; i++) {
      id[i] = Npart*rank + i;
  }
  
#ifdef H5PART
  h5part_int64_t ng = bdata.ng();
  h5part_int64_t ng2d = bdata.ng2d();
  h5part_int64_t np = bdata.np();
  h5part_int64_t rL = bdata.rL();
  np *= np * np;
  write_hcosmo_h5(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
		  id, Npart, outBase, ng, ng2d, np, rL, indatName);

  /* 
     Only for debug purposes, should be in the same format than the
     mc2 fort.66
   */
  if (NumProcs==1) 
    write_hcosmo_fort(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
		      Npart, create_outName(outBase,-1).c_str() );
#else
  write_hcosmo(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
	       id, Npart, create_outName(outBase,rank).c_str() );
#endif

  free(pos_x);
  free(pos_y);
  free(pos_z);
  free(vel_x);
  free(vel_y);
  free(vel_z);

  free(id);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return(0);
}

