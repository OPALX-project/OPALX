// ------------------------------------------------------------------------
// $RCSfile: DataSink.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DataSink
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2004/06/02 19:38:54 $
// $Author: adelmann $
// $Log: DataSink.cpp,v $
// Revision 1.3  2004/06/02 19:38:54  adelmann
//
// ------------------------------------------------------------------------

#include <iomanip>
#include <string>
#include <hdf5.h>
#include "H5Part.h"
#include "DataSink.h"
#include "../ChargedParticles/ChargedParticles.hh"
extern "C" {
#include "H5Block.h"
#include "H5BlockTypes.h"
}

DataSink::DataSink(string fn, ChargedParticles<TT,3> *univ)
{
    univ_m = univ;

    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");

    /// Set file write flags to true. These will be set to false after first
    /// write operation.
    firstWriteH5part_m = true;

    int pos=fn.find(string("."), 0);
    fn.erase(pos, fn.size() - pos);

    fn += string(".h5");

    /// Open H5 file. Check that it opens correctly.
#ifdef PARALLEL_IO
    H5file_m=H5PartOpenFileParallel(fn.c_str(),H5PART_WRITE,MPI_COMM_WORLD);
#else
    H5file_m=H5PartOpenFile(fn.c_str(),H5PART_WRITE);
#endif

    if(!H5file_m) {
        ERRORMSG("h5 file open failed: exiting!" << endl);
        exit(0);
    }

    /// Write file attributes.
    writeH5FileAttributes();

    /// Set current record/time step to 0.
    H5call_m = 0;

}

DataSink::DataSink(string fn, ChargedParticles<TT,3> *univ, int restartStep)
{
    univ_m = univ;
    
    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");

    /// Set file write flags. Since this is a restart, assume H5 file is old.
    firstWriteH5part_m = false;

#ifdef PARALLEL_IO
    H5file_m = H5PartOpenFileParallel(fn.c_str(),H5PART_READ,MPI_COMM_WORLD);
#else
    H5file_m = H5PartOpenFile(fn.c_str(),H5PART_READ);
#endif

    if(!H5file_m) {
        ERRORMSG("h5 file open failed: exiting!" << endl);
        exit(0);
    }
    
    //currently always restart from last step
    int numStepsInFile = H5PartGetNumSteps(H5file_m);
    restartStep = numStepsInFile;
    
    // Use same dump frequency.
    //int dumpfreq = 0;
    //H5PartReadFileAttrib(H5file_m, "dump frequency", &dumpfreq);
    
    // Set current record/time step to restart step.
    H5call_m = restartStep-1;

    readPhaseSpace();
}

DataSink::~DataSink()
{
    if (H5file_m)
	H5PartCloseFile(H5file_m);
    Ippl::Comm->barrier();
}

void DataSink::writeH5FileAttributes() {

    /// Write file attributes to describe phase space to H5 file.
    H5PartWriteFileAttribString(H5file_m,"tUnit","s");
    H5PartWriteFileAttribString(H5file_m,"xUnit","Mpc/h");
    H5PartWriteFileAttribString(H5file_m,"yUnit","Mpc/h");
    H5PartWriteFileAttribString(H5file_m,"zUnit","Mpc/h");
    H5PartWriteFileAttribString(H5file_m,"pxUnit","km/s");
    H5PartWriteFileAttribString(H5file_m,"pyUnit","km/s");
    H5PartWriteFileAttribString(H5file_m,"pzUnit","km/s");
    H5PartWriteFileAttribString(H5file_m,"idUnit","1");

    H5PartWriteFileAttribString(H5file_m,"TIMEUnit","s");

    /// Write file dump frequency.
    //TODO: do we need this?
    //h5part_int64_t dumpfreq = Options::psDumpFreq;
    h5part_int64_t dumpfreq = 1;
    H5PartWriteFileAttrib(H5file_m, "dump frequency", H5PART_INT64, &dumpfreq, 1);

    /// Reset first write flag.
    firstWriteH5part_m = false;
}

void DataSink::writePhaseSpace(TT time, TT z, int step) {
    Inform msg("writePhaseSpace ");
    msg << " number " << H5call_m << endl;
    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);
 
    size_t nTot                   = univ_m->getTotalNum();
    size_t nLoc                   = univ_m->getLocalNum();

#ifdef IPPL_USE_SINGLE_PRECISION
    void *varray = malloc(nLoc*sizeof(float));
    float *farray = (float*)varray;
    void *lvarray = malloc(nLoc*sizeof(double));
    h5part_int64_t *larray = (h5part_int64_t *)lvarray;
#else
    void *varray = malloc(nLoc*sizeof(double));
    double *farray = (double*)varray;
    h5part_int64_t *larray = (h5part_int64_t *)varray;
#endif


    ///Get the particle decomposition from all the compute nodes.
    size_t *locN = (size_t *) malloc(Ippl::getNodes()*sizeof(size_t));
    size_t  *globN = (size_t*) malloc(Ippl::getNodes()*sizeof(size_t));

    for(int i=0; i<Ippl::getNodes(); i++) {
        globN[i] = locN[i]=0;
    }

    locN[Ippl::myNode()] = nLoc;
    reduce(locN, locN + Ippl::getNodes(), globN, OpAddAssign());

    /// Set current record/time step.
    H5PartSetStep(H5file_m, H5call_m);
    H5PartSetNumParticles(H5file_m, nLoc);

    H5PartWriteStepAttrib(H5file_m,"SPOS",     H5T_NATIVE_DOUBLE, &z, 1);
    H5PartWriteStepAttrib(H5file_m,"TIME",     H5T_NATIVE_DOUBLE, &time,1);

    /// Write number of compute nodes.
    H5PartWriteStepAttrib(H5file_m,"nloc",H5T_NATIVE_INT64, globN, Ippl::getNodes());

    /// Write univ phase space.
    for(size_t i=0; i<nLoc; i++)
        farray[i] =  univ_m->R[i](0);
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"x",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"x",farray);
#endif

    for(size_t i=0; i<nLoc; i++)
        farray[i] =  univ_m->R[i](1);
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"y",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"y",farray);
#endif

    for(size_t i=0; i<nLoc; i++)
        farray[i] =  univ_m->R[i](2);
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"z",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"z",farray);
#endif

    for(size_t i=0; i<nLoc; i++)
        farray[i] =  univ_m->V[i](0);
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"px",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"px",farray);
#endif

    for(size_t i=0; i<nLoc; i++)
        farray[i] =  univ_m->V[i](1);
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"py",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"py",farray);
#endif

    for(size_t i=0; i<nLoc; i++)
        farray[i] =  univ_m->V[i](2);
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"pz",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"pz",farray);
#endif

    /// Write particle id numbers.
    for (size_t i = 0; i < nLoc; i++)
        larray[i] =  univ_m->ID[i];
    H5PartWriteDataInt64(H5file_m,"id",larray);

    H5Fflush(H5file_m->file,H5F_SCOPE_GLOBAL);

    /// Step record/time step index.
    H5call_m++;

    if(varray)
        free(varray);
#ifdef IPPL_USE_SINGLE_PRECISION
    if(lvarray)
        free(lvarray);
#endif

    /// Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}


void DataSink::readPhaseSpace()
{

    H5PartSetStep(H5file_m, H5call_m);
    int N = (int)H5PartGetNumParticles(H5file_m);

    //TODO: do a more sophisticated distribution of particles?
    //my guess is that the end range index is EXCLUSIVE!
    int numberOfParticlesPerNode = (int) floor((double) N / Ippl::getNodes());
    long long starti = Ippl::myNode() * numberOfParticlesPerNode;
    long long endi = 0;
    // ensure that we don't miss any particle in the end
    if(Ippl::myNode() == Ippl::getNodes() - 1)
        endi = -1;
    else
        endi = starti + numberOfParticlesPerNode;

    H5PartSetView(H5file_m,starti,endi);
    N = (int)H5PartGetNumParticles(H5file_m);

    //double actualT;
    //H5PartReadStepAttrib(H5file,"TIME",&actualT);

#ifdef IPPL_USE_SINGLE_PRECISION
    void *varray = malloc(N*sizeof(float));
    float *farray = (float*)varray;
#else
    void *varray = malloc(N*sizeof(double));
    double *farray = (double*)varray;
#endif

    univ_m->create(N);

#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartReadDataFloat32(H5file_m,"x",farray);
#else
    H5PartReadDataFloat64(H5file_m,"x",farray);
#endif
    for (unsigned long int n=0; n < N; ++n) {
        univ_m->R[n](0)=farray[n];
	univ_m->M[n]=1.0;
    }

#ifdef IPPL_USE_SINGLE_PRECISION    
    H5PartReadDataFloat32(H5file_m,"y",farray);
#else
    H5PartReadDataFloat64(H5file_m,"y",farray);
#endif
    for (unsigned long int n=0; n < N; ++n)
        univ_m->R[n](1)=farray[n];

#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartReadDataFloat32(H5file_m,"z",farray);
#else
    H5PartReadDataFloat64(H5file_m,"z",farray);
#endif
    for (unsigned long int n=0; n < N; ++n)
        univ_m->R[n](2)=farray[n];

#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartReadDataFloat32(H5file_m,"px",farray);
#else
    H5PartReadDataFloat64(H5file_m,"px",farray);
#endif
    for (unsigned long int n=0; n < N; ++n)
        univ_m->V[n](0)=farray[n];

#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartReadDataFloat32(H5file_m,"py",farray);
#else
    H5PartReadDataFloat64(H5file_m,"py",farray);
#endif
    for (unsigned long int n=0; n < N; ++n)
        univ_m->V[n](1)=farray[n];

#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartReadDataFloat32(H5file_m,"pz",farray);
#else
    H5PartReadDataFloat64(H5file_m,"pz",farray);
#endif
    for (unsigned long int n=0; n < N; ++n)
        univ_m->V[n](2)=farray[n];

    if(farray)
        free(farray);

    H5PartCloseFile(H5file_m);

    return;
}


/***************************************************************************
 * $RCSfile: DataSink.cpp,v $   $Author: adelmann $
 * $Revision: 1.3 $   $Date: 2004/06/02 19:38:54 $
 ***************************************************************************/

