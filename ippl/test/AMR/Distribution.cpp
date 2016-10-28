#include "Distribution.h"

#include "H5Reader.h"


// ----------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// ----------------------------------------------------------------------------

void Distribution::uniform(double lower, double upper, size_t nloc, int seed) {
    
    nloc_m = nloc;
    
    std::mt19937_64 mt(seed);
    
    std::uniform_real_distribution<> dist(lower, upper);
    
    x_m.resize(nloc);
    y_m.resize(nloc);
    z_m.resize(nloc);
    
    px_m.resize(nloc);
    py_m.resize(nloc);
    pz_m.resize(nloc);
    
    for (size_t i = 0; i < nloc; ++i) {
        x_m[i] = dist(mt);
        y_m[i] = dist(mt);
        z_m[i] = dist(mt);
        
        px_m[i] = dist(mt);
        py_m[i] = dist(mt);
        pz_m[i] = dist(mt);
    }
}

void Distribution::gaussian(double mean, double stddev, size_t nloc, int seed) {
    
    nloc_m = nloc;
    
    std::mt19937_64 mt(seed);
    
    std::normal_distribution<double> dist(mean, stddev);
    
    x_m.resize(nloc);
    y_m.resize(nloc);
    z_m.resize(nloc);
    
    px_m.resize(nloc);
    py_m.resize(nloc);
    pz_m.resize(nloc);
    
    for (size_t i = 0; i < nloc; ++i) {
        x_m[i] = dist(mt);
        y_m[i] = dist(mt);
        z_m[i] = dist(mt);
        
        px_m[i] = dist(mt);
        py_m[i] = dist(mt);
        pz_m[i] = dist(mt);
    }
}

void Distribution::readH5(const std::string& filename, int step) {
    H5Reader h5(filename);
    
    h5.open(step);
    
    
    
    long nloc_m = h5.getNumParticles();
    size_t numParticlesPerNode = nloc_m / Ippl::getNodes();

    size_t firstParticle = numParticlesPerNode * Ippl::myNode();
    size_t lastParticle = firstParticle + numParticlesPerNode - 1;
    if (Ippl::myNode() == Ippl::getNodes() - 1)
        lastParticle = nloc_m - 1;

    nloc_m = lastParticle - firstParticle + 1;
    
    x_m.resize(nloc_m);
    y_m.resize(nloc_m);
    z_m.resize(nloc_m);
    
    px_m.resize(nloc_m);
    py_m.resize(nloc_m);
    pz_m.resize(nloc_m);
    
    h5.read(x_m, px_m,
            y_m, py_m,
            z_m, pz_m,
            firstParticle,
            lastParticle);
    
    h5.close();
}


void Distribution::injectBeam(PartBunchBase& bunch) {
    
    // create memory space
    bunch.create(nloc_m);
    
    for (int i = bunch.getLocalNum() - 1; i >= 0; --i) {
        bunch.setR(Vector_t(x_m[i], y_m[i], z_m[i]), i);
        bunch.setP(Vector_t(px_m[i], py_m[i], pz_m[i]), i);
        
        x_m.pop_back();
        y_m.pop_back();
        z_m.pop_back();
        px_m.pop_back();
        py_m.pop_back();
        pz_m.pop_back();
    }
}