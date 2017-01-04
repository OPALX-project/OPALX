#include "H5Reader.h"

#include <vector>
#include <cassert>

#include "Ippl.h"

H5Reader::H5Reader(const std::string& filename)
    : filename_m(filename), file_m(0)
{ }

H5Reader::H5Reader()
    : filename_m(""), file_m(0)
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
	y[n] -= int(y[n]);
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


void H5Reader::writeScalarField(const container_t& scalfield,
                                const Array<Geometry>& geom)
{
    
    std::string fname = "test_scalfield.h5";
    h5_file_t file = H5OpenFile (fname.c_str(), H5_O_WRONLY, H5_PROP_DEFAULT);
    H5SetStep (file, 0);
    
    for (int l = 0; l < scalfield.size(); ++l) {
        int ii = 0;
        int gridnr = 0;
        for (MFIter mfi(scalfield[l]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const FArrayBox& field = (scalfield[l])[mfi];
            
            H5Block3dSetView(file,
                             bx.loVect()[0], bx.hiVect()[0],
                             bx.loVect()[1], bx.hiVect()[1],
                             bx.loVect()[2], bx.hiVect()[2]);
            
            std::cout << "#points: " << bx.numPts()
                      << " i_dim: " << bx.loVect()[0] << " " << bx.hiVect()[0]
                      << " j_dim: " << bx.loVect()[1] << " " << bx.hiVect()[1]
                      << " k_dim: " << bx.loVect()[2] << " " << bx.hiVect()[2]
                      << std::endl;
            
            std::unique_ptr<h5_float64_t[]> data(
                new h5_float64_t[bx.numPts()]
            );
            
            // Fortran storing convention, i.e. column-major
            for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                    for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                        IntVect ivec(i, j, k);
//                         std::cout << ii << " " << ivec << std::endl; std::cin.get();
                        data[ii] = field(ivec, 0 /*component*/);
                        ++ii;
                    }
                }
            }
            std::string group = "rho/level-" + std::to_string(l) + "/grid-" + std::to_string(gridnr);
//             std::string group = "rho";
            H5Block3dWriteScalarFieldFloat64(file, group.c_str(), data.get());
            
            RealBox rb = geom[l].ProbDomain();
            
            H5Block3dSetFieldOrigin(file, group.c_str(),
                                    (h5_float64_t)rb.lo(0),
                                    (h5_float64_t)rb.lo(1),
                                    (h5_float64_t)rb.lo(2));
            
            H5Block3dSetFieldSpacing(file, group.c_str(),
                                     (h5_float64_t)(geom[0].CellSize(0)),
                                     (h5_float64_t)(geom[0].CellSize(1)),
                                     (h5_float64_t)(geom[0].CellSize(2)));
            
            ++gridnr;
        }
    }
    
    H5CloseFile (file);
}


void H5Reader::writeVectorField(const container_t& vecfield) {
    
}