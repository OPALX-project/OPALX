#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

/*!
 * @file Distribution.h
 * @details Generates distributions
 * for a particle bunch.
 * @authors Matthias Frey \n
 *          Andreas Adelmann \n
 *          Ann Almgren \n
 *          Weiqun Zhang
 * @date LBNL, October 2016
 */

#include <random>
#include <iostream>

#include "PartBunchBase.h"


/// Create particle distributions
class Distribution {
    
public:
    typedef std::vector<double> container_t;
    
    
public:
    
    /// Generate an uniform particle distribution
    /*!
     * @param lower boundary
     * @param upper boundary
     * @param nloc is the local number of particles
     * @param seed of the Mersenne-Twister
     */
    void uniform(double lower, double upper, size_t nloc, int seed);
    
    /// Generate a Gaussian particle distribution
    /*!
     * @param mean also called centroid
     * @param stddev is the standard deviation
     * @param nloc is the local number of particles
     * @param seed of the Mersenne-Twister
     */
    void gaussian(double mean, double stddev, size_t nloc, int seed);
    
    /// Transfer distribution to particle bunch object.
    /*! @param bunch is either an AmrPartBunch or an PartBunch object
     */
    void injectBeam(PartBunchBase& bunch);
    
private:
    container_t x_m;    ///< Horizontal particle positions [m]
    container_t y_m;    ///< Vertical particle positions [m]
    container_t z_m;    ///< Longitudinal particle positions [m]
    size_t nloc_m;      ///< Local number of particles
};



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
    
    for (size_t i = 0; i < nloc; ++i) {
        x_m[i] = dist(mt);
        y_m[i] = dist(mt);
        z_m[i] = dist(mt);
    }
}

void Distribution::gaussian(double mean, double stddev, size_t nloc, int seed) {
    
    nloc_m = nloc;
    
    std::mt19937_64 mt(seed);
    
    std::uniform_real_distribution<> dist(mean, stddev);
    
    x_m.resize(nloc);
    y_m.resize(nloc);
    z_m.resize(nloc);
    
    for (size_t i = 0; i < nloc; ++i) {
        x_m[i] = dist(mt);
        y_m[i] = dist(mt);
        z_m[i] = dist(mt);
    }
}


void Distribution::injectBeam(PartBunchBase& bunch) {
    
    // create memory space
    bunch.create(nloc_m);
    
    for (size_t i = 0; i < bunch.getLocalNum(); ++i)
        bunch.setR(Vector_t(x_m[i], y_m[i], z_m[i]), i);
}

#endif