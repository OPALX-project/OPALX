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
#include <array>

#include <Array.H>
#include <Geometry.H>

#ifdef IPPL_AMR
    #include "ippl-amr/AmrParticleBase.h"
    #include "ippl-amr/ParticleAmrLayout.h"
    #include "ippl-amr/PartBunchAmr.h"
#else
    #include "boxlib-amr/PartBunchBase.h"
#endif

/// Create particle distributions
class Distribution {
    
public:
    typedef std::vector<double> container_t;
    typedef Vektor<double, BL_SPACEDIM> Vector_t;
    
public:
    
    Distribution();
    
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
    
    /// Generate a two stream instability distribution according to B.\ Ulmer
    /*!
     * @param lower boundary of domain
     * @param upper boundary of domain
     * @param nx is the number of grid points in cooordinate space
     * @param nv is the number of grid points in velocity space
     * @param vmax is the max. velocity
     * @param alpha is the amplitude of the initial disturbance
     */
    void twostream(const Vector_t& lower, const Vector_t& upper,
                   const Vektor<std::size_t, 3>& nx, const Vektor<std::size_t, 3>& nv,
                   const Vektor<double, 3>& vmax, double alpha = 0.5);
    
    /// Generate a uniform particle disitribution per cell
    /*!
     * Each cell of the domain receives the same number of particles.
     * @param geom is the geometry obtained from Amr (for each level).
     * @param ba are all boxes per level
     * @param nr is the number of gridpoints in each direction of coarsest level
     * @param nParticles is the number of particles per cell
     * @param seed of the Mersenne-Twister
     */
    void uniformPerCell(const Array<Geometry>& geom,
                        const Array<BoxArray>& ba,
                        const Vektor<std::size_t, 3>& nr,
                        std::size_t nParticles, int seed);
    
    /// Read in a distribution from an H5 file
    /*!
     * @param filename is the path and name of the H5 file.
     * @param step to be read in.
     */
    void readH5(const std::string& filename, int step);
    
    /// Transfer distribution to particle bunch object.
    /*! @param bunch is either an AmrPartBunch or an PartBunch object
     * @param doDelete removes all particles already in bunch before
     * injection.
     * @param shift all particles, each direction independently
     */
    void injectBeam(
#ifdef IPPL_AMR
        PartBunchAmr< ParticleAmrLayout<double, BL_SPACEDIM> > & bunch,
#else
        PartBunchBase& bunch,
#endif
        bool doDelete = true, std::array<double, 3> shift = {{0.0, 0.0, 0.0}});
    
    /// Update a distribution (only single-core)
    /*! @param bunch is either an AmrPartBunch or an PartBunch object
     * @param filename is the path and name of the H5 file
     * @param step to be read in from a H5 file
     */
    void setDistribution(
#ifdef IPPL_AMR
        PartBunchAmr< ParticleAmrLayout<double, BL_SPACEDIM> >& bunch,
#else
        PartBunchBase& bunch,
#endif
        const std::string& filename, int step);
    
    /// Write the particles to a text file that can be read by OPAL. (sec. 11.3 in OPAL manual)
    /*!
     * @param pathname where to store.
     */
    void print2file(std::string pathname);
    
private:
    container_t x_m;    ///< Horizontal particle positions [m]
    container_t y_m;    ///< Vertical particle positions [m]
    container_t z_m;    ///< Longitudinal particle positions [m]
    
    container_t px_m;   ///< Horizontal particle momentum
    container_t py_m;   ///< Vertical particle momentum
    container_t pz_m;   ///< Longitudinal particle momentum
    container_t q_m;    ///< Particle charge (always set to 1.0, except for Distribution::readH5)
    container_t mass_m; ///< Particl mass (always set to 1.0, except for Distribution::readH5 and Distribution::twostream)
    
    size_t nloc_m;      ///< Local number of particles
    size_t ntot_m;      ///< Total number of particles
};

#endif
