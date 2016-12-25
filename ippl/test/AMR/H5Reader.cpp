#include "H5Reader.h"

#include <vector>
#include <cassert>

#include "Ippl.h"

H5Reader::H5Reader(const std::string& filename)
    : filename_m(filename), file_m(0)
{ }

void H5Reader::open(int step) {
    close();
    
    
#if defined (USE_H5HUT2)
    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = Ippl::getComm();
    h5_err_t h5err = H5SetPropFileMPIOCollective (props, &comm);
#if defined (NDEBUG)
    (void)h5err;
#endif
    assert (h5err != H5_ERR);
    file_m = H5OpenFile (filename_m.c_str(), H5_O_RDONLY, props);
    assert (file_m != (h5_file_t)H5_ERR);
#else
    file_m = H5OpenFile(filename_m.c_str(), H5_FLUSH_STEP | H5_O_RDONLY, Ippl::getComm());
    assert (file_m != (void*)H5_ERR);
#endif
    
    H5SetStep(file_m, step);
}


void H5Reader::close() {
    if (file_m) {
        Ippl::Comm->barrier();
        
        H5CloseFile(file_m);
        file_m = 0;
    }
}


void H5Reader::read(Distribution::container_t& x,
                    Distribution::container_t& px,
                    Distribution::container_t& y,
                    Distribution::container_t& py,
                    Distribution::container_t& z,
                    Distribution::container_t& pz,
                    Distribution::container_t& q,
                    size_t firstParticle,
                    size_t lastParticle)
{
    // cyclotron: y <--> z
    
    h5_ssize_t numParticles = getNumParticles();
    
//     H5PartSetNumParticles(numParticles);
    
    H5PartSetView(file_m, firstParticle, lastParticle);

    numParticles = lastParticle - firstParticle + 1;

    std::vector<char> buffer(numParticles * sizeof(h5_float64_t));
    h5_float64_t *f64buffer = reinterpret_cast<h5_float64_t*>(&buffer[0]);

    READDATA(Float64, file_m, "x", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        x[n] = f64buffer[n];
        
        x[n] -= int(x[n]);
    }

    READDATA(Float64, file_m, "y", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        z[n] = f64buffer[n];
        
        // do shift due to box
        z[n] -= int(z[n]);
    }

    READDATA(Float64, file_m, "z", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        y[n] = f64buffer[n];
    }

    READDATA(Float64, file_m, "px", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        px[n] = f64buffer[n];
    }

    READDATA(Float64, file_m, "py", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        pz[n] = f64buffer[n];
    }

    READDATA(Float64, file_m, "pz", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        py[n] = f64buffer[n];
    }

//     READDATA(Float64, file_m, "q", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        q[n] = 1.0; //f64buffer[n];
    }
/*
    READDATA(Float64, file_m, "mass", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        bunch.M[n] = f64buffer[n];
    }*/
}

h5_ssize_t H5Reader::getNumParticles() {
    return H5PartGetNumParticles(file_m);
}