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
 * @brief Generate different particle distributions
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
    
    /// Read in a distribution from an H5 file
    /*!
     * @param filename is the path and name of the H5 file.
     * @param step to be read in.
     */
    void readH5(const std::string& filename, int step);
    
    /// Transfer distribution to particle bunch object.
    /*! @param bunch is either an AmrPartBunch or an PartBunch object
     */
    void injectBeam(PartBunchBase& bunch);
    
private:
    container_t x_m;    ///< Horizontal particle positions [m]
    container_t y_m;    ///< Vertical particle positions [m]
    container_t z_m;    ///< Longitudinal particle positions [m]
    
    container_t px_m;   ///< Horizontal particle momentum
    container_t py_m;   ///< Vertical particle momentum
    container_t pz_m;   ///< Longitudinal particle momentum
    
    size_t nloc_m;      ///< Local number of particles
};

#endif