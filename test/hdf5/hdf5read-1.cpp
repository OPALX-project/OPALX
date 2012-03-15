#include <hdf5.h>
#include "Ippl.h"
#include <hdf5.h>
#include "H5Part.hh"
#include "H5Block.hh"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
int main(int argc,char *argv[]){
  Ippl ippl(argc, argv);
  Inform msg("hdf5read ");  
  Inform msg2all("hdf5read ",INFORM_ALL_NODES);  
  Timer mytimer;

  H5PartSetVerbosityLevel(0);

  H5PartFile *file;

  /* parallel file read */
  file=H5PartOpenFileParallel("dataH5.dat",H5PART_READ,MPI_COMM_WORLD);
  if(!file) {
    ERRORMSG("File open failed:  exiting!" << endl);
    exit(0);
  }

  H5PartSetStep(file,0);      

  int nt,nds;
  char name[64];
  double magNumRx = 0.0;
  double magNumRy = 0.0;
  double magNumRz = 0.0;

  nt =H5PartGetNumSteps(file);     /* get number of steps in file */
  nds=H5PartGetNumDatasets(file);  /* get number of datasets in timestep 0 */ 

  for(int i=0;i<nds;i++){          /* and print out those names */
    H5PartGetDatasetName(file,i,name,64);
    msg << "Dataset " << i << " name= " << string(name) << endl;
  }
  msg << "Number of datasteps in file is " << nt << endl;
  
  long long  *globN = (long long*) malloc(Ippl::getNodes()*sizeof(long long));


  for(int i=0;i<nt;i++){ 
    H5PartSetStep(file,i);

    long long actStep = 0;
    long long starti = 0;
    long long endi   = 0;

    H5PartReadStepAttrib(file,"Step",&actStep);
    H5PartReadStepAttrib(file,"nloc",globN);
    
    if (Ippl::myNode() != 0) {
      for(int i=0; i< Ippl::myNode(); i++) 
	starti += globN[i];
    }
    endi = starti + globN[Ippl::myNode()];
    H5PartSetView(file,starti,endi);
    
    int N=(int)H5PartGetNumParticles(file);
    void *varray = malloc(N*8);
    double *farray = (double*)varray;
      
    mytimer.clear();    
    mytimer.start();
    H5PartReadDataFloat64(file,"x",farray);
    mytimer.stop();
    for(int k=0;k<N;k++)
      magNumRx += farray[k];
    reduce(magNumRx,magNumRx,OpAddAssign());
  
    mytimer.start();
    H5PartReadDataFloat64(file,"y",farray);
    mytimer.stop();
    for(int k=0;k<N;k++)
      magNumRy += farray[k];
    reduce(magNumRy,magNumRy,OpAddAssign());

    mytimer.start();
    H5PartReadDataFloat64(file,"z",farray);
    mytimer.stop();
    for(int k=0;k<N;k++)
      magNumRz += farray[k];
    reduce(magNumRz,magNumRz,OpAddAssign());      

    double rate = (3*N*sizeof(double)/1000000.0) / mytimer.clock_time();
  
    reduce(rate,rate,OpAddAssign());      
    reduce(N,N,OpAddAssign());      
    
    msg << "Number of particles " << N << " in file set " << i
	<< " magic number x = " << magNumRx << " magic number y = " << magNumRy << " magic number z = " << magNumRz 
	<< " Read from disk br= " << rate << " [MB/s] "<< endl;      

    magNumRx = magNumRy = magNumRz = 0.0;      

    // ada read field
    
    h5part_int64_t herr;
    h5part_float64_t *data;
    h5part_int64_t l[6];

    stringstream lstr;
    lstr << "layout" << Ippl::myNode();

    herr = H5BlockReadFieldAttrib (file,"EFmag", lstr.str().c_str(),l); 
    if (herr < 0)
      exit(1);
    
    h5part_int64_t dataLen = (l[1]-l[0]+1) * (l[3]-l[2]+1) * (l[5]-l[4]+1) ;
    
    data = (h5part_float64_t *) malloc ( dataLen * sizeof ( *data ) );
    
    herr = H5BlockDefine3DFieldLayout (file, l[0], l[1], l[2], l[3], l[4], l[5]);
    if (herr < 0)
      msg << "H5BlockDefine3DFieldLayout (file, l[0], l[1], l[2], l[3], l[4], l[5]) failed" << endl;
    
    herr = H5Block3dReadScalarField ( file, "EFmag", data );
    if (herr < 0)
      msg << "H5BlockReadScalarField  failed" << endl;
    
    double sum = 0.0;
    for( h5part_int64_t i = 0; i<dataLen; i++)
      sum+= data[i];
    reduce(sum,sum, OpAddAssign());

    msg << "sum= " << sum << endl;
    
    if(farray)
      free(farray);

    if(data)
      free(data);

  }

  if(globN)
    free(globN);
  
  H5PartCloseFile(file);
  
  Ippl::Comm->barrier();
  return 0;
}
