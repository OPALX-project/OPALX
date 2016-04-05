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
#include <stdexcept>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>

#include "matrix_vector_operation.h"

#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FMatrix.h"

/// @brief This class generates the matrices for the one turn matrix of a cyclotron.
template<typename Value_type, typename Size_type, typename Series_type,
         typename Map_type, typename Hamiltonian_type, typename Space_charge_type>
class MapGenerator
{
    public:
        /// Type of variables
        typedef Value_type value_type;
        /// Type for specifying sizes
        typedef Size_type size_type;
        /// Type of truncated power series
        typedef Series_type series_type;
        /// Type of a map
        typedef Map_type map_type;
        /// Type of the Hamiltonian
        typedef Hamiltonian_type hamiltonian_type;
        /// Type of the Hamiltonian representing the space charge
        typedef Space_charge_type space_charge_type;
        
        /// Type for specifying matrices
        typedef boost::numeric::ublas::matrix<value_type> matrix_type;
        /// Type for specifying vectors
        typedef std::vector<value_type> vector_type;
        
        /// Constructor
        /*!
         * @param nMaps is the number of maps
         */
        MapGenerator(size_type);

        /// Generates a map based on the Hamiltonian for a given angle
        /*!
         * @param H represents the Hamiltonian
         * @param ds is the step size (angle dependent)
         * @param order is the truncation order of the Taylor series of the exponential function
         */
        matrix_type generateMap(const series_type&, value_type, size_type);
        
        /// Combines the space charge maps (for each angle one) and the cyclotron maps (for each angle one) to the ont turn map, taking lists of maps
        /*!
         * @param Mscs is a list of space charge maps (the higher the index, the higher the angle)
         * @param Mcycs is a list of cyclotron maps (the higher the index, the higher the angle)
         */
        void combine(std::vector<matrix_type>&, std::vector<matrix_type>&);
        
        /// Returns the one turn map
        matrix_type getMap();
        
    private:
        /// Number of maps
        size_type nMaps_m;
        /// One-turn matrix
        matrix_type Mturn_m;
//         /// Stores the one turn transfer map
//         map_type M_m;
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------


template<typename Value_type, typename Size_type, typename Series_type, typename Map_type, typename Hamiltonian_type, typename Space_charge_type>
MapGenerator<Value_type, Size_type, Series_type, Map_type, Hamiltonian_type, Space_charge_type>::MapGenerator(size_type nMaps)
: nMaps_m(nMaps)
{}

template<typename Value_type, typename Size_type, typename Series_type, typename Map_type, typename Hamiltonian_type, typename Space_charge_type>
typename MapGenerator<Value_type, Size_type, Series_type, Map_type, Hamiltonian_type, Space_charge_type>::matrix_type MapGenerator<Value_type, Size_type, Series_type, Map_type, Hamiltonian_type, Space_charge_type>::generateMap/*cyc*/(const series_type& H, value_type ds, size_type order) {
    
    // expand
    map_type M = ExpMap(-H * ds, order);
    
    // get linear part
    FMatrix<value_type,6,6> matrix = M.linearTerms();
    
    matrix_type map = boost::numeric::ublas::zero_matrix<value_type>(6);
    
    for (size_type i = 0; i < 6; ++i) {
        for (size_type j = 0; j < 6; ++j) {
            map(i,j) = matrix[i][j];
        }
    }
    
    return map;
}

template<typename Value_type, typename Size_type, typename Series_type, typename Map_type, typename Hamiltonian_type, typename Space_charge_type>
void MapGenerator<Value_type, Size_type, Series_type, Map_type, Hamiltonian_type, Space_charge_type>::combine(std::vector<matrix_type>& Mscs, std::vector<matrix_type>& Mcycs) {
    
    if (nMaps_m != Mscs.size() || nMaps_m != Mcycs.size())
        throw std::length_error("Error in MapGenerator::combine: Wrong vector dimensions.");

    Mturn_m = boost::numeric::ublas::identity_matrix<value_type>(6);
    
    for (size_type i = 0; i < nMaps_m; ++i)
        Mturn_m = matt_boost::gemmm<matrix_type>(Mscs[i],Mcycs[i],Mturn_m);
}

template<typename Value_type, typename Size_type, typename Series_type, typename Map_type, typename Hamiltonian_type, typename Space_charge_type>
typename MapGenerator<Value_type, Size_type, Series_type, Map_type, Hamiltonian_type, Space_charge_type>::matrix_type MapGenerator<Value_type, Size_type, Series_type, Map_type, Hamiltonian_type, Space_charge_type>::getMap() {
    return Mturn_m;
}

#endif
