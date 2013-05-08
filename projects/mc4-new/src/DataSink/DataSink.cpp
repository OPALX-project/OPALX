#include <iomanip>
#include <string>
#include <hdf5.h>
#include "H5hut.h"
#include "DataSink.h"
#include "../ChargedParticles/ChargedParticles.hh"
#include "H5Block.h"


// backward compatibility with 1.99.5
#if defined(H5_O_FLUSHSTEP)
#define H5_FLUSH_STEP H5_O_FLUSHSTEP
#endif

DataSink::DataSink(std::string fn, ChargedParticles<TT,3> *univ)
{
    univ_m = univ;

    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");

    /// Set file write flags to true. These will be set to false after first
    /// write operation.
    firstWriteH5part_m = true;

    int pos=fn.find(std::string("."), 0);
    fn.erase(pos, fn.size() - pos);

    fn += std::string(".h5");

#ifdef PARALLEL_IO
    H5file_m = H5OpenFile(fn.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, MPI_COMM_WORLD);
#else
    H5file_m = H5OpenFile(fn.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, 0);
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

DataSink::DataSink(std::string fn, ChargedParticles<TT,3> *univ, int restartStep)
{
    univ_m = univ;

    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");

    /// Set file write flags. Since this is a restart, assume H5 file is old.
    firstWriteH5part_m = false;

#ifdef PARALLEL_IO
    H5file_m = H5OpenFile(fn.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, MPI_COMM_WORLD);
#else
    H5file_m = H5OpenFile(fn.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, 0);
#endif

    if(!H5file_m) {
        ERRORMSG("h5 file open failed: exiting!" << endl);
        exit(0);
    }

    //currently always restart from last step
    int numStepsInFile = H5GetNumSteps(H5file_m);
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
        H5CloseFile(H5file_m);
    Ippl::Comm->barrier();
}

void DataSink::writeH5FileAttributes() {

    /// Write file attributes to describe phase space to H5 file.
    H5WriteFileAttribString(H5file_m,"tUnit","s");
    H5WriteFileAttribString(H5file_m,"xUnit","Mpc/h");
    H5WriteFileAttribString(H5file_m,"yUnit","Mpc/h");
    H5WriteFileAttribString(H5file_m,"zUnit","Mpc/h");
    H5WriteFileAttribString(H5file_m,"pxUnit","km/s");
    H5WriteFileAttribString(H5file_m,"pyUnit","km/s");
    H5WriteFileAttribString(H5file_m,"pzUnit","km/s");
    H5WriteFileAttribString(H5file_m,"idUnit","1");

    H5WriteFileAttribString(H5file_m,"TIMEUnit","s");

    /// Write file dump frequency.
    //TODO: do we need this?
    //h5part_int64_t dumpfreq = Options::psDumpFreq;
    h5_int64_t dumpfreq = 1;
    H5WriteFileAttribInt64(H5file_m, "dump frequency", &dumpfreq, 1);

    /// Reset first write flag.
    firstWriteH5part_m = false;
}

void DataSink::writePhaseSpaceNeutrinos(TT time, TT z, int step, size_t nNeutr) {
    Inform msg("dumpUniverse ");
    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);

    size_t nLoc                   = nNeutr; 

#ifdef IPPL_USE_SINGLE_PRECISION
    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(float)]);
    float *farray = reinterpret_cast<float *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());
#else
    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());
#endif

    ///Get the particle decomposition from all the compute nodes.
    std::unique_ptr<size_t[]> locN(new size_t[Ippl::getNodes()]);
    std::unique_ptr<size_t[]> globN(new size_t[Ippl::getNodes()]);

    for(int i=0; i<Ippl::getNodes(); i++) {
        globN[i] = locN[i]=0;
    }

    locN[Ippl::myNode()] = nLoc;
    reduce(locN.get(), locN.get() + Ippl::getNodes(), globN.get(), OpAddAssign());

    /// Set current record/time step.
    H5SetStep(H5file_m, H5call_m);
    H5PartSetNumParticles(H5file_m, nLoc);

    H5WriteStepAttribFloat64(H5file_m,"SPOS", &z, 1);
    H5WriteStepAttribFloat64(H5file_m,"TIME", &time,1);

    /// Write number of compute nodes.
    H5WriteStepAttribInt64(H5file_m,"nloc", (h5_int64_t *)globN.get(), Ippl::getNodes());

    /// Write univ phase space.
    size_t k=0;
    for(size_t i=0; i<nLoc; i++) {
        if (univ_m->M[i] == -1.0) {
            farray[k] =  univ_m->R[i](0);
            k++;
        }
    }
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"x",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"x",farray);
#endif
    k=0;
    for(size_t i=0; i<nLoc; i++) {
        if (univ_m->M[i] == -1.0) {
            farray[k] =  univ_m->R[i](1);
            k++;
        }
    }
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"y",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"y",farray);
#endif

    k=0;
    for(size_t i=0; i<nLoc; i++) {
        if (univ_m->M[i] == -1.0) {
            farray[k] =  univ_m->R[i](2);
            k++;
        }
    }
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"z",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"z",farray);
#endif

    k=0;
    for(size_t i=0; i<nLoc; i++){
        if (univ_m->M[i] == -1.0) {
            farray[k] =  univ_m->V[i](0);
            k++;
        }
    }
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"px",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"px",farray);
#endif

    k=0;
    for(size_t i=0; i<nLoc; i++) {
        if (univ_m->M[i] == -1.0) {
            farray[k] =  univ_m->V[i](1);
            k++;
        }
    }
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"py",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"py",farray);
#endif
    k=0;
    for(size_t i=0; i<nLoc; i++) {
        if (univ_m->M[i] == -1.0) {
            farray[k] =  univ_m->V[i](2);
            k++;
        }
    }
#ifdef IPPL_USE_SINGLE_PRECISION
    H5PartWriteDataFloat32(H5file_m,"pz",farray);
#else
    H5PartWriteDataFloat64(H5file_m,"pz",farray);
#endif

    /// Write particle id numbers.
    k=0;
    for (size_t i = 0; i < nLoc; i++){
        if (univ_m->M[i] == -1.0) {
            larray[k] =  univ_m->ID[i];
            k++;
        }
    }
    H5PartWriteDataInt64(H5file_m,"id",larray);

    //    H5Fflush(H5file_m->file,H5F_SCOPE_GLOBAL);

    msg << " step number " << H5call_m << " np " << k << endl;

    /// Step record/time step index.
    H5call_m++;

    IpplTimings::stopTimer(H5PartTimer_m);

}

void DataSink::writePhaseSpace(TT time, TT z, int step) {
    Inform msg("dumpUniverse ");
    msg << " step number " << H5call_m << endl;
    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);

    //size_t nTot                   = univ_m->getTotalNum();
    size_t nLoc                   = univ_m->getLocalNum();

#ifdef IPPL_USE_SINGLE_PRECISION
    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(float)]);
    float *farray = reinterpret_cast<float *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());
#else
    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());
#endif

    ///Get the particle decomposition from all the compute nodes.
    std::unique_ptr<size_t[]> locN(new size_t[Ippl::getNodes()]);
    std::unique_ptr<size_t[]> globN(new size_t[Ippl::getNodes()]);


    for(int i=0; i<Ippl::getNodes(); i++) {
        globN[i] = locN[i]=0;
    }

    locN[Ippl::myNode()] = nLoc;
    reduce(locN.get(), locN.get() + Ippl::getNodes(), globN.get(), OpAddAssign());

    /// Set current record/time step.
    H5SetStep(H5file_m, H5call_m);
    H5PartSetNumParticles(H5file_m, nLoc);

    H5WriteStepAttribFloat64(H5file_m,"SPOS", &z, 1);
    H5WriteStepAttribFloat64(H5file_m,"TIME", &time,1);

    /// Write number of compute nodes.
    H5WriteStepAttribInt64(H5file_m,"nloc", (h5_int64_t *)globN.get(), Ippl::getNodes());

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

    //    H5Fflush(H5file_m->file,H5F_SCOPE_GLOBAL);

    /// Step record/time step index.
    H5call_m++;

    /// Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}


void DataSink::readPhaseSpace()
{

    H5SetStep(H5file_m, H5call_m);
    unsigned long int N = static_cast<unsigned long int>(H5PartGetNumParticles(H5file_m));

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
    N = static_cast<unsigned long int>(H5PartGetNumParticles(H5file_m));

    //double actualT;
    //H5PartReadStepAttrib(H5file,"TIME",&actualT);

#ifdef IPPL_USE_SINGLE_PRECISION
    std::unique_ptr<char[]> varray(new char[(N)*sizeof(float)]);
    float *farray = reinterpret_cast<float *>(varray.get());
#else
    std::unique_ptr<char[]> varray(new char[(N)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
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

    H5CloseFile(H5file_m);

    return;
}

