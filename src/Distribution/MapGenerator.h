/**
 * @file MapGenerator.h
 * This class is based on the paper "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons" (2012)
 * of Dr. Christian Baumgarten.
 * It has one template parameter that specifies the type of the variables and containers.
 *
 * @author Matthias Frey
 * @version 1.0
 */
#ifndef MAPGENERATOR_H
#define MAPGENERATOR_H

#include <cmath>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>

#include "matrix_vector_operation.h"
#include "physics.h"

/// @brief This class generates the matrices for the one turn matrix of a cyclotron.
template<typename Value_type, typename Size_type>
class MapGenerator
{
    public:
        /// Type of variables
        typedef Value_type value_type;
        /// Type for specifying sizes
        typedef Size_type size_type;
        /// Type for specifying matrices
        typedef boost::numeric::ublas::matrix<value_type> matrix_type;
        /// Type for specifying vectors
        typedef std::vector<value_type> vector_type;

        /// Constructor
        /*!
         * @param N is the number of angle splits
         * @param I is the current, \f$ \left[I\right] = A \f$
         * @param gamma is the relativistic factor, \f$ \left[\gamma\right] = 1 \f$
         * @param wo is the orbital frequency, \f$ \left[\omega_{o}\right] = \frac{1}{s} \f$
         * @param nh is the harmonic number, \f$ \left[N_{h}\right] = 1 \f$
         * @param m is the mass of the particles \f$ [m] = \frac{MeV}{c^{2}} \f$
         */
        MapGenerator(size_type, value_type, value_type, value_type, value_type, value_type);

        /// Computes the space charge map for a given angle
        /*!
         * @param sigma is the sigma matrix
         * @param ds is the step size (angle dependent)
         */
        matrix_type generateMsc(const matrix_type&, value_type);

        /// Computes the cyclotron map for a given angle
        /*!
         * @param h is the inverse bending radius (angle dependent), \f$ \left[h\right] = \frac{1}{m} \f$
         * @param fidx is the field index (angle dependent)
         * @param ds is the step size (angle dependent)
         * @param order specifies the order of the Taylor expansion for the cyclotron map
         */
        matrix_type generateMcyc(value_type,value_type,value_type,size_type);

        /// Combines the space charge maps (for each angle one) and the cyclotron maps (for each angle one) to the one turn map
        /*!
         * @param h is vector of inverse bending radiuses for each angle, \f$ \left[h\right] = \frac{1}{m} \f$
         * @param fidx is a vector of field indexes for each angle
         * @param ds is vector of path lengths for each angle
         * @param sigma is a list of sigma matrices for each angle
         */
        void combine(const vector_type&, const vector_type&, const vector_type&, /*const*/ std::vector<matrix_type>&, size_type);

        /// Combines the space charge maps (for each angle one) and the cyclotron maps (for each angle one) to the ont turn map, taking lists of maps
        /*!
         * @param Mscs is a list of space charge maps (the higher the index, the higher the angle)
         * @param Mcycs is a list of cyclotron maps (the higher the index, the higher the angle)
         */
        void combine(std::vector<matrix_type>&, std::vector<matrix_type>&);

        /// Returns the one turn map
        matrix_type getMap();

    private:
        /// Number of angle splits
        size_type N_m;
        /// Relativistic factor
        value_type gamma_m;
        /// Relativistic factor squared
        value_type gamma2_m;
        /// One-turn matrix
        matrix_type Mturn_m;
        /// Constant, \f$ \left[K_{3}\right] = m \f$
        value_type K3_m;
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
MapGenerator<Value_type, Size_type>::MapGenerator(size_type N, value_type I, value_type gamma, value_type wo, value_type nh, value_type m)
: N_m(N), gamma_m(gamma), gamma2_m(gamma*gamma), Mturn_m(6,6)
{
    value_type beta = std::sqrt(1.0 - 1.0 / gamma2_m);

    // convert m from MeV/c^2 to eV*s^{2}/m^{2}
    m *= 1.0e6 / (physics::c * physics::c);

    // formula (57)
    value_type lam = 2.0 * M_PI*physics::c / (wo * nh); // wavelength, [lam] = m
    K3_m = 3.0 * physics::q0 * I * lam / (20.0 * std::sqrt(5.0) * M_PI * physics::eps0 * m *
           physics::c * physics::c * physics::c * beta * beta * gamma2_m * gamma_m);            // [K3_m] = m
}

template<typename Value_type, typename Size_type>
typename MapGenerator<Value_type, Size_type>::matrix_type MapGenerator<Value_type, Size_type>::generateMsc(const matrix_type& sigma, value_type ds) {
    /*
     * Computes space charge matrix according to 
     * formula (20) and (21) from paper:
     * Transverse-Longitudinal Coupling by Space Charge in Cyclotrons
     * 
     * The functions "__M6SC" and "__M6sc" of mtx.c of Dr. Christian Baumgarten are used as implementation guide.
     */
    // space charge map for one angle
    matrix_type Msc = boost::numeric::ublas::identity_matrix<value_type>(6);

    value_type milli = 1.0e-3;

    // formula (30), (31)
    // [sigma(0,0)] = mm^{2} rad --> [sx] = [sy] = [sz] = mm
    // multiply with 0.001 to get meter --> [sx] = [sy] = [sz] = m 
    value_type sx = std::sqrt(std::fabs(sigma(0,0))) * milli;
    value_type sy = std::sqrt(std::fabs(sigma(2,2))) * milli;
    value_type sz = std::sqrt(std::fabs(sigma(4,4))) * milli;

    value_type tmp = sx * sy;                                           // [tmp] = m^{2}
    
    value_type f = std::sqrt(tmp) / (3.0 * gamma_m * sz);               // [f] = 1
    value_type kxy = K3_m * std::fabs(1.0 - f) * ds / ((sx + sy) * sz); // [kxy] = 1/m
    
    Msc(1,0) = kxy / sx;                                                // [Msc(1,0)] = 1/m^{2}
    Msc(3,2) = kxy / sy;                                                // [Msc(3,2)] = 1/m^{2}
    Msc(5,4) = K3_m * f * ds * gamma2_m / (tmp * sz);                   // [Msc(5,4)] = 1/m^{2}
    
    return Msc;
}

template<typename Value_type, typename Size_type>
typename MapGenerator<Value_type, Size_type>::matrix_type MapGenerator<Value_type, Size_type>::generateMcyc(value_type h,value_type fidx,
        value_type ds, size_type order) {

    /*
     * Computes the cyclotron matrix according to formula
     * (20) and (44) from the paper:
     * Transverse-Longitudinal Coupling by Space Charge in Cyclotrons
     * 
     * The implementation in the function "scMakeCycFromFieldData" of ring_ana5.c
     * of Dr. Christian Baumgarten is used as implementation guide
     */

    // force matrix for one angle
    matrix_type F = boost::numeric::ublas::zero_matrix<value_type>(6);
    F(0,1) = 1.0;
    F(1,0) = - h * h - fidx;    // p. 6 of paper: "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons"
    F(1,5) = h;
    F(2,3) = 1.0;
    F(3,2) = fidx;              // p. 6 of paper: "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons"
    F(4,0) = - h;
    F(4,5) = 1.0 / gamma2_m;

    // cyclotron map for one angle;
    matrix_type Mcyc = matt_boost::taylor_exp<value_type>(F,ds,order);
    
    return Mcyc;
}

template<typename Value_type, typename Size_type>
void MapGenerator<Value_type, Size_type>::combine(const vector_type& h, const vector_type& fidx, const vector_type& ds,
        /*const*/ std::vector<matrix_type>& sigma, size_type order) {

    Mturn_m = boost::numeric::ublas::identity_matrix<value_type>(6);

    for (size_type i = 0; i < N_m; ++i)
        Mturn_m = matt_boost::gemmm<matrix_type>(generateMsc(sigma[i],ds[i]),generateMcyc(h[i],fidx[i],ds[i],order),Mturn_m);
}

template<typename Value_type, typename Size_type>
void MapGenerator<Value_type, Size_type>::combine(std::vector<matrix_type>& Mscs, std::vector<matrix_type>& Mcycs) {
    
    if (N_m != Mscs.size() || N_m != Mcycs.size())
        Error::message("MapGenerator<Value_type, Size_type>:combine(std::vector<matrix_type>&, std::vector<matrix_type>&", Error::size);

    Mturn_m = boost::numeric::ublas::identity_matrix<value_type>(6);
    
    for (size_type i = 0; i < N_m; ++i)
        Mturn_m = matt_boost::gemmm<matrix_type>(Mscs[i],Mcycs[i],Mturn_m);
}

template<typename Value_type, typename Size_type>
typename MapGenerator<Value_type, Size_type>::matrix_type MapGenerator<Value_type, Size_type>::getMap() {
    return Mturn_m;
}

#endif
