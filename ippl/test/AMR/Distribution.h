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

#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>

#ifdef IPPL_AMR
    #include "ippl-amr/AmrParticleBase1.h"
    #include "ippl-amr/ParticleAmrLayout.h"
    #include "ippl-amr/PartBunchAmr.h"
#else
    #include "boxlib-amr/PartBunchBase.h"
#endif

using amrex::Vector;
using amrex::BoxArray;
using amrex::Geometry;

/// Create particle distributions
class Distribution {

public:
    typedef std::vector<double> container_t;
    typedef Vektor<double, AMREX_SPACEDIM> Vector_t;
    typedef std::function<double(const Vector_t&,
                         const Vector_t&,
                         double,
                         double)> function_t;
    
    // Only used in the function "special"
    enum Type {
        kTwoStream,
        kLandauDamping,
        kRecurrence
    };

private:
    function_t twoStream = [&](const Vector_t& pos,
                               const Vector_t& vel,
                               double alpha, double k)
    {
#if AMREX_SPACEDIM == 3
        double factor = 1.0 / ( M_PI * 30.0 );
        double v2 = vel[0] * vel[0] +
                    vel[1] * vel[1] +
                    vel[2] * vel[2];
        
        double f = factor * std::exp(-0.5 * v2) *
                    (1.0 + alpha * std::cos(k * pos[2])) *
                    (1.0 + 5.0/*0.5*/ * vel[2] * vel[2]);
#elif AMREX_SPACEDIM == 2
        double factor = 1.0 / ( 12.0 * M_PI );
        double v2 = vel[0] * vel[0] +
                    vel[1] * vel[1];
        
        double f = factor * std::exp(-0.5 * v2) *
                    (1.0 + alpha * std::cos(k * pos[0])) *
                    (1.0 + 5.0 * vel[0] * vel[0]);
#else
        double f = 0.0;
        throw std::runtime_error("DIM > 1");
#endif
        return f;
    };
    
    function_t landauDamping = [&](const Vector_t& pos,
                                   const Vector_t& vel,
                                   double alpha, double k)
    {
#if AMREX_SPACEDIM == 3
        double factor = 1.0 / ( 2.0 * M_PI * std::sqrt(2.0 * M_PI) );
        double v2 = vel[0] * vel[0] +
                    vel[1] * vel[1] +
                    vel[2] * vel[2];
        
        double f = factor * std::exp( -0.5 * v2) *
                    (1.0 + alpha * ( std::cos( k * pos[2] ) +
                                     std::cos( k * pos[1] ) +
                                     std::cos( k * pos[0])
                                   )
                    );
#elif AMREX_SPACEDIM == 2
        double factor = 1.0 / ( 2.0 * M_PI );
        double v2 = vel[0] * vel[0] +
                    vel[1] * vel[1];
        
        double f = factor * std::exp(-0.5 * v2) *
                    (1.0 + alpha * ( std::cos(k * pos[0]) *
                                     std::cos(k * pos[1])
                                   )
                    );
#else
        double f = 0.0;
        throw std::runtime_error("DIM > 1");
#endif
        return f;
    };
    
    
    // free-streaming
    function_t recurrence = [&](const Vector_t& pos,
                               const Vector_t& vel,
                               double alpha, double k)
    {
#if AMREX_SPACEDIM == 3
        double factor = 1.0 / ( 2.0 * M_PI * std::sqrt(2.0 * M_PI) );
        double v2 = vel[0] * vel[0] +
                    vel[1] * vel[1] +
                    vel[2] * vel[2];
        
        double f = alpha * factor * std::exp(-0.5 * v2) *
                    (
                        std::cos( k * pos[2] ) +
                        std::cos( k * pos[1] ) +
                        std::cos( k * pos[0] )
                    );
#else
        double f = 0.0;
        throw std::runtime_error("DIM > 2");
#endif
        return f;
    };
    
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

    /// More general version of Gaussian particle distribution
    /*!
     * @param mean for each direction independently
     * @param stddev for each direction
     * @param nloc local number of particles
     * @param seed of Mersenne-Twister
     */
    void gaussian(const double* mean, const double* stddev,
		  size_t nloc, int seed);
    
    /// Generate particle distributions according to B.\ Ulmer
    /*!
     * Supported distributions:\n
     * - 2-stream instability
     * - Landau damping
     * - recurrence
     * @param lower boundary of domain
     * @param upper boundary of domain
     * @param nx is the number of grid points in cooordinate space
     * @param nv is the number of grid points in velocity space
     * @param vmax is the max. velocity
     * @param type is either kTwoStream, kLandauDamping or kRecurrence
     * @param alpha is the amplitude of the initial disturbance
     * @param kk is phase
     */
    void special(const Vector_t& lower, const Vector_t& upper,
                   const Vektor<std::size_t, AMREX_SPACEDIM>& nx,
                 const Vektor<std::size_t, AMREX_SPACEDIM>& nv,
                   const Vektor<double, AMREX_SPACEDIM>& vmax,
                 Type type, double alpha = 0.5, double kk = 0.5);
    
    /// Generate a uniform particle disitribution per cell
    /*!
     * Each cell of the domain receives the same number of particles.
     * @param geom is the geometry obtained from Amr (for each level).
     * @param ba are all boxes per level
     * @param nr is the number of gridpoints in each direction of coarsest level
     * @param nParticles is the number of particles per cell
     * @param seed of the Mersenne-Twister
     */
    void uniformPerCell(const Vector<Geometry>& geom,
                        const Vector<BoxArray>& ba,
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
        PartBunchAmr< ParticleAmrLayout<double, AMREX_SPACEDIM> > & bunch,
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
        PartBunchAmr< ParticleAmrLayout<double, AMREX_SPACEDIM> >& bunch,
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
