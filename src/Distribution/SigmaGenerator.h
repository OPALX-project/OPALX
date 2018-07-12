/**
 * @file SigmaGenerator.h
 * The SigmaGenerator class uses the class <b>ClosedOrbitFinder</b> to get the parameters (inverse bending radius, path length
 * field index and tunes) to initialize the sigma matrix. 
 * The main function of this class is <b>match(value_type, size_type)</b>, where it iteratively
 * tries to find a matched distribution for given emittances, energy and current. 
 * The computation stops when the L2-norm is smaller than a user-defined tolerance. \n
 * In default mode it prints all space charge maps, cyclotron maps and second moment matrices. The orbit properties, i.e.
 * tunes, average radius, orbit radius, inverse bending radius, path length, field index and frequency error, are printed
 * as well.
 *
 * @author Matthias Frey
 * @version 1.0
 * @version 2.1 (implemented by Cristopher Cortes)
 * The original version used a private class <b>RDM</b> to calculate the coupled part of the sigma matrix
 * in version 2.1 this was replaced by an Eigenvalue Solver
 * @version 2.2 (implemented by Cristopher Cortes)
 * In version 2.2 for every permutation of the Eigenvalue-order the algorithm searches for a matched distribution
 * @version 2.3 (implemented by Cristopher Cortes)
 * The algorithm puts the vertical plane on the right position, and the other eigenvectors can be permutated
 */
#ifndef SIGMAGENERATOR_H
#define SIGMAGENERATOR_H

#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <limits>
#include <list>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <complex>
#include <cstdlib>

#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include "Utilities/OpalException.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#if BOOST_VERSION >= 106000
#include <boost/numeric/odeint/integrate/check_adapter.hpp>
#endif

#include <boost/filesystem.hpp>

#include "matrix_vector_operation.h"
#include "ClosedOrbitFinder.h"
#include "MapGenerator.h"
#include "Harmonics.h"

extern Inform *gmsg;

/// @brief This class computes the matched distribution
template<typename Value_type, typename Size_type>
class SigmaGenerator
{
    public:
        /// Type of variables
        typedef Value_type value_type;
        /// Type of constant variables
        typedef const value_type const_value_type;
        /// Type for specifying sizes
        typedef Size_type size_type;
	/// Type for storing complex numbers
	typedef std::complex<value_type> complex_number_type;
        /// Type for storing maps
        typedef boost::numeric::ublas::matrix<value_type> matrix_type;
        /// Type for storing the sparse maps
        typedef boost::numeric::ublas::compressed_matrix<value_type,boost::numeric::ublas::row_major> sparse_matrix_type;
        /// Type for storing vectors
        typedef boost::numeric::ublas::vector<value_type> vector_type;
	/// Type for storing complex matrices
	typedef boost::numeric::ublas::matrix<complex_number_type> complex_matrix_type;
	/// Type for storing complex vectors
	typedef boost::numeric::ublas::vector<complex_number_type> complex_vector_type;

	/// Type for storing gsl real Matrices
	typedef gsl_matrix* gsl_matrix_type;
	/// Type for storing gsl complex Vectors
	typedef gsl_vector_complex* gsl_vector_complex_type;
        /// Type for storing gsl complex Matrices
	typedef gsl_matrix_complex* gsl_matrix_complex_type;
	/// Type for storing gsl workspaces
	typedef gsl_eigen_nonsymmv_workspace* gsl_eigen_nonsymmv_workspace_type;
	/// Type for storing gsl complex numbers
	typedef gsl_complex gsl_complex_type;
	/// Type for storing a permutation matrix
	typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix_type;
	
	/// Container for storing the properties for each angle
        typedef std::vector<value_type> container_type;
        /// Type of the truncated powere series
        typedef FTps<value_type,2*3> Series;
        /// Type of a map
        typedef FVps<value_type,2*3> Map;
        /// Type of the Hamiltonian for the cyclotron
        typedef std::function<Series(value_type,value_type,value_type)> Hamiltonian;
        /// Type of the Hamiltonian for the space charge
        typedef std::function<Series(value_type,value_type,value_type)> SpaceCharge;

        /// Constructs an object of type SigmaGenerator
        /*!
         * @param I specifies the current for which a matched distribution should be found, \f$ [I] = A \f$
         * @param ex is the emittance in x-direction (horizontal), \f$ \left[\varepsilon_{x}\right] = \pi\ mm\ mrad  \f$
         * @param ey is the emittance in y-direction (longitudinal), \f$ \left[\varepsilon_{y}\right] = \pi\ mm\ mrad  \f$
         * @param ez is the emittance in z-direction (vertical), \f$ \left[\varepsilon_{z}\right] = \pi\ mm\ mrad  \f$
         * @param wo is the orbital frequency, \f$ \left[\omega_{o}\right] = \frac{1}{s} \f$
         * @param E is the energy, \f$ \left[E\right] = MeV \f$
         * @param nh is the harmonic number
         * @param m is the mass of the particles \f$ \left[m\right] = \frac{MeV}{c^{2}} \f$
         * @param Emin is the minimum energy [MeV] needed in cyclotron, \f$ \left[E_{min}\right] = MeV \f$
         * @param Emax is the maximum energy [MeV] reached in cyclotron, \f$ \left[E_{max}\right] = MeV \f$
         * @param nSector is the number of sectors (symmetry assumption)
         * @param N is the number of integration steps (closed orbit computation). That's why its also the number
         *        of maps (for each integration step a map)
         * @param fieldmap is the location of the file that specifies the magnetic field
         * @param truncOrder is the truncation order for power series of the Hamiltonian
         * @param scaleFactor for the magnetic field (default: 1.0)
         * @param write is a boolean (default: true). If true all maps of all iterations are stored, otherwise not.
         */
        SigmaGenerator(value_type I, value_type ex, value_type ey, value_type ez,
                       value_type wo, value_type E, value_type nh, value_type m,
                       value_type Emin, value_type Emax, size_type nSector,
                       size_type N, const std::string& fieldmap,
                       size_type truncOrder, value_type scaleFactor = 1.0, bool write = true);

        /// Searches for a matched distribution, for every permutation of the Eigenvalues
        /*!
         * Returns "true" if at least one matched distribution was found
         * @param accuracy is used for computing the equilibrium orbit (to a certain accuracy)
         * @param maxit is the maximum number of iterations (the program stops if within this value no stationary
         *        distribution was found)
         * @param maxitOrbit is the maximum number of iterations for finding closed orbit
         * @param angle defines the start of the sector (one can choose any angle between 0° and 360°)
	 * @param guess value of radius for closed orbit finder
         * @param type specifies the magnetic field format (e.g. CARBONCYCL)
         * @param harmonic is a boolean. If "true" the harmonics are used instead of the closed orbit finder.
         */
        bool match(value_type accuracy, size_type maxit, size_type maxitOrbit, value_type angle,
                   value_type guess, const std::string& type, bool harmonic);

	/// Returns True if an Inverse Matrix was found
	/*
	 * Uses lu_factorize and lu_substitute in uBLAS to invert the matrix
	 * Source: Numerical Recipies in C, 2nd ed., by Press, Teukolsky, Vetterling & Flannery
	 * (It stores the inverse matrix in invR) 
	 */
	bool InvertMatrix(const complex_matrix_type& R, complex_matrix_type& invR);
	
        /// Checks if the sigma-matrix is an eigenellipse and returns the L2 error.
        /*!
         * The idea of this function is taken from Dr. Christian Baumgarten's program.
         * @param Mturn is the one turn transfer matrix
         * @param sigma is the sigma matrix to be tested
         */
        value_type isEigenEllipse(const matrix_type& Mturn, const matrix_type& sigma);

        /// Returns the converged stationary distribution
        matrix_type& getSigma();

        /// Returns the number of iterations needed for convergence (if not converged, it returns zero)
        size_type getIterations() const;

        /// Returns the error (if the program didn't converged, one can call this function to check the error)
        value_type getError() const;

        /// Returns the emittances (ex,ey,ez) in \f$ \pi\ mm\ mrad \f$ for which the code converged 
	///(since the whole simulation is done on normalized emittances)
        std::array<value_type,3> getEmittances() const;

    private:
        /// Stores the value of the current, \f$ \left[I\right] = A \f$
        value_type I_m;
        /// Stores the desired emittances, \f$ \left[\varepsilon_{x}\right] = \left[\varepsilon_{y}\right] = \left[\varepsilon_{z}\right] = mm \ mrad \f$
        std::array<value_type,3> emittance_m;
        /// Is the orbital frequency, \f$ \left[\omega_{o}\right] = \frac{1}{s} \f$
        value_type wo_m;
        /// Stores the user-define energy, \f$ \left[E\right] = MeV \f$
        value_type E_m;
        /// Relativistic factor (which can be computed out ot the kinetic energy and rest mass (potential energy)), \f$ \left[\gamma\right] = 1 \f$
        value_type gamma_m;
        /// Relativistic factor squared
        value_type gamma2_m;
        /// Harmonic number, \f$ \left[N_{h}\right] = 1 \f$
        value_type nh_m;
        /// Velocity (c/v), \f$ \left[\beta\right] = 1 \f$
        value_type beta_m;
        /// Is the mass of the particles, \f$ \left[m\right] = \frac{MeV}{c^{2}} \f$
        value_type m_m;
        /// Is the number of iterations needed for convergence
        size_type niterations_m;
        /// Is true if converged, false otherwise
        bool converged_m;
        /// Minimum energy needed in cyclotron, \f$ \left[E_{min}\right] = MeV \f$
        value_type Emin_m;
        /// Maximum energy reached in cyclotron, \f$ \left[E_{max}\right] = MeV \f$
        value_type Emax_m;
        /// Number of (symmetric) sectors
        size_type nSector_m;
        /// Number of integration steps for closed orbit computation
        size_type N_m;
        /// Number of integration steps per sector (--> also: number of maps per sector)
        size_type nStepsPerSector_m;
        /// Error of computation
        value_type error_m;

        /// Location of magnetic fieldmap
        std::string fieldmap_m;

        /// Truncation order of the power series for the Hamiltonian (cyclotron and space charge)
        size_type truncOrder_m;

        /// Decides for writing output (default: true)
        bool write_m;
        
        /// Scale the magnetic field values (default: 1.0)
        value_type scaleFactor_m;

        /// Stores the stationary distribution (sigma matrix)
        matrix_type sigma_m;
	/// Stores the found matched distributions
	std::vector<matrix_type> matchedSigmas_m;
        /// Vector storing the second moment matrix for each angle
        std::vector<matrix_type> sigmas_m;

        /// Vector containing the possible permutations
	std::vector<std::vector<size_type>> permutations_m;

        /// Stores the Hamiltonian of the cyclotron
        Hamiltonian H_m;

        /// Stores the Hamiltonian for the space charge
        SpaceCharge Hsc_m;

        /// All variables x, px, y, py, z, delta
        Series x_m, px_m, y_m, py_m, z_m, delta_m;
	
	/// Stores the dimension of the matrix 
	unsigned short const matrixDimension{6};
	

        /// Initializes a first guess of the sigma matrix with the assumption of a spherical symmetric beam (ex = ey = ez). For each angle split the same initial guess is taken.
        /*!
         * @param nuz is the vertical tune
         * @param ravg is the average radius of the closed orbit
         */
        void initialize(value_type, value_type);
	
	/// Searches for a matched distribution for a given permutation
	/*
	 * Returns true if a matched distribution was found within the accuracy constrain, returns false otherwise
	 * @param accuracy is used for computing the equilibrium orbit (to a certain accuracy)
	 * @param h contains the inverse radius for each step
	 * @param fidx contains the angle for each step
	 * @param ds contains each step ds
	 * @param tunes contains the computed tunes
	 * @param ravg contains the average radius for the computed close orbit
	 * @param maxit is the maximum number of iterations (the program stops if within this value no stationary
         *        distribution was found)
	 * @param harmonic is a boolean. If "true" the harmonics are used instead of the closed orbit finder.(not implemented yet)
	 * @param permutation is the permutation of the order of the Eigenvalues (default: 0) 
	 *
	 */
	
	bool findMatchDistribution(value_type accuracy,
				   container_type& h,
				   container_type& fidx,
				   container_type& ds,
				   std::pair<value_type,value_type> tunes,
				   value_type ravg,
				   size_type maxit,
				   std::vector<matrix_type>& Mcycs,
				   bool harmonic,
				   size_type permutation = 0);
	
		
	/// Sets the Eigenvectors with a given order in R and its inverse in invR
	/*
	 * 
	 * @param Mturn is the one turn transfer matrix
	 * @param R is the matrix containing the Eigenvectors of Mturn
	 * @param invR is the inverse of the R matrix
	 */
	bool setEigenVectors(const matrix_type& Mturn, complex_matrix_type& R, complex_matrix_type& invR, size_type perm);

	/// Computes the new initial sigma matrix (Version 2.0)
	/*!
         * @param R is the transformation matrix containig the Eigenvalues of M
         * @param invR is the inverse transformation matrix
         */
	matrix_type updateInitialSigma(const complex_matrix_type&, const complex_matrix_type&);

        /// Computes new sigma matrices (one for each angle)
        /*!
         * Mscs is a vector of all space charge maps
         * Mcycs is a vector of all cyclotron maps
         */
        void updateSigma(const std::vector<matrix_type>&, const std::vector<matrix_type>&);
	
	/// Writes the Matched Distributions
	/*
	 * @param match is the matched distribution sigma matrix
	 * @param permutation is the permutation used
	 * @param Mcycs is the container with all the cyclotron maps for each angle
	 * @param Mscs is the container with the space charge maps for each angle
	 */
	void writeMatched(const matrix_type& match, 
			  const size_type permutation,
			  const std::vector<matrix_type>& Mcycs,
			  const std::vector<matrix_type>& Mscs,
			  const matrix_type& Mcyc);
        /// Returns the L2-error norm between the old and new sigma-matrix
        /*!
         * @param oldS is the old sigma matrix
         * @param newS is the new sigma matrix
         */
        value_type L2ErrorNorm(const matrix_type&, const matrix_type&);
        
        /// Returns the L1-error norm between the old and new sigma-matrix
        /*!
         * @param oldS is the old sigma matrix
         * @param newS is the new sigma matrix
         */
        value_type L1ErrorNorm(const matrix_type&, const matrix_type&);

        /// Transforms a floating point value to a string
        /*!
         * @param val is the floating point value which is transformed to a string
         */
        std::string float2string(value_type val);
        
        /// Called within SigmaGenerator::match().
        /*!
         * @param tunes
         * @param ravg is the average radius [m]
         * @param r_turn is the radius [m]
         * @param peo is the momentum
         * @param h_turn is the inverse bending radius
         * @param fidx_turn is the field index
         * @param ds_turn is the path length element
         */
        void writeOrbitOutput_m(const std::pair<value_type,value_type>& tunes,
                                const value_type& ravg,
                                const value_type& freqError,
                                const container_type& r_turn,
                                const container_type& peo,
                                const container_type& h_turn,
                                const container_type&  fidx_turn,
                                const container_type&  ds_turn);
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
SigmaGenerator<Value_type, Size_type>::SigmaGenerator(value_type I, value_type ex, value_type ey, value_type ez, value_type wo,
        value_type E, value_type nh, value_type m, value_type Emin, value_type Emax, size_type nSector,
        size_type N, const std::string& fieldmap, size_type truncOrder, value_type scaleFactor, bool write)
: I_m(I), wo_m(wo), E_m(E),
  gamma_m(E/m+1.0), gamma2_m(gamma_m*gamma_m),
  nh_m(nh), beta_m(std::sqrt(1.0-1.0/gamma2_m)), m_m(m), niterations_m(0), converged_m(false),
  Emin_m(Emin), Emax_m(Emax), nSector_m(nSector), N_m(N), nStepsPerSector_m(N/nSector),
  error_m(std::numeric_limits<value_type>::max()),
  fieldmap_m(fieldmap), truncOrder_m(truncOrder), write_m(write),
  scaleFactor_m(scaleFactor), sigmas_m(nStepsPerSector_m),permutations_m{{0,1,2,3},{1,2,3,0},{0,1,3,2},{1,0,2,3},{0,2,3,1},{0,3,2,1},
  {1,0,3,2},{1,3,2,0},{2,1,0,3},{2,0,1,3},{2,3,0,1},{2,3,1,0},{3,1,0,2},{3,2,0,1}, {3,2,1,0},{3,0,1,2}}
  // The permutations up converge to a matched distribution, but the first four give the 4 unique solutions
  // Following permutations threw an OPAL error, whose reason I still don't know. 
  //{0,2,1,3},{0,3,1,2},{1,2,0,3},{1,3,0,2},{2,1,3,0},{2,0,3,1},{3,1,2,0},{3,0,2,1}
{
    // set emittances (initialization like that due to old compiler version)
    // [ex] = [ey] = [ez] = pi*mm*mrad --> [emittance] = mm mrad
    emittance_m[0] = ex * Physics::pi;
    emittance_m[1] = ey * Physics::pi;
    emittance_m[2] = ez * Physics::pi;

    // minimum beta*gamma
    value_type minGamma = Emin_m / m_m + 1.0;
    value_type bgam = std::sqrt(minGamma * minGamma - 1.0);

    // normalized emittance (--> multiply with beta*gamma)
    // [emittance] = mm mrad
    emittance_m[0] *= bgam;
    emittance_m[1] *= bgam;
    emittance_m[2] *= bgam;

    // Define the Hamiltonian
    Series::setGlobalTruncOrder(truncOrder_m);

    // infinitesimal elements
    x_m = Series::makeVariable(0);
    px_m = Series::makeVariable(1);
    y_m = Series::makeVariable(2);
    py_m = Series::makeVariable(3);
    z_m = Series::makeVariable(4);
    delta_m = Series::makeVariable(5);

    H_m = [&](value_type h, value_type kx, value_type ky) {
        return 0.5*px_m*px_m + 0.5*kx*x_m*x_m - h*x_m*delta_m +
               0.5*py_m*py_m + 0.5*ky*y_m*y_m +
               0.5*delta_m*delta_m/gamma2_m;
    };

    Hsc_m = [&](value_type sigx, value_type sigy, value_type sigz) {
        // convert m from MeV/c^2 to eV*s^{2}/m^{2}
        value_type m = m_m * 1.0e6 / (Physics::c * Physics::c);

        // formula (57)
        value_type lam = 2.0 * Physics::pi*Physics::c / (wo_m * nh_m); // wavelength, [lam] = m
        value_type K3 = 3.0 * /* physics::q0 */ 1.0 * I_m * lam / (20.0 * std::sqrt(5.0) * Physics::pi * Physics::epsilon_0 * m *
                        Physics::c * Physics::c * Physics::c * beta_m * beta_m * gamma_m * gamma2_m);            // [K3] = m

        value_type milli = 1.0e-3;

        // formula (30), (31)
        // [sigma(0,0)] = mm^{2} rad --> [sx] = [sy] = [sz] = mm
        // multiply with 0.001 to get meter --> [sx] = [sy] = [sz] = m
        value_type sx = std::sqrt(std::fabs(sigx)) * milli;
        value_type sy = std::sqrt(std::fabs(sigy)) * milli;
        value_type sz = std::sqrt(std::fabs(sigz)) * milli;

        value_type tmp = sx * sy;                                           // [tmp] = m^{2}

        value_type f = std::sqrt(tmp) / (3.0 * gamma_m * sz);               // [f] = 1
        value_type kxy = K3 * std::fabs(1.0 - f) / ((sx + sy) * sz); // [kxy] = 1/m

        value_type Kx = kxy / sx;
        value_type Ky = kxy / sy;
        value_type Kz = K3 * f / (tmp * sz);

        return -0.5 * Kx * x_m * x_m
               -0.5 * Ky * y_m * y_m
               -0.5 * Kz * z_m * z_m * gamma2_m;
     };
}

template<typename Value_type, typename Size_type>
  bool SigmaGenerator<Value_type, Size_type>::match(value_type accuracy,
                                                    size_type maxit,
                                                    size_type maxitOrbit,
                                                    value_type angle,
                                                    value_type rguess,
                                                    const std::string& type,
                                                    bool harmonic)
{   
    /* compute the equilibrium orbit for energy E
     * and get the the following properties:
     * - inverse bending radius h
     * - step sizes of path ds
     * - tune nuz
     */
    // object for space cyclotron map
    try{
        MapGenerator<value_type,size_type,Series,Map,Hamiltonian,SpaceCharge> mapgen(nStepsPerSector_m);

	container_type h(nStepsPerSector_m), r(nStepsPerSector_m), ds(nStepsPerSector_m), fidx(nStepsPerSector_m);
	value_type ravg = 0.0;
	std::pair<value_type,value_type> tunes;

	std::string resultDir = "data";
	if ( !boost::filesystem::exists(resultDir) ) {
	    boost::filesystem::create_directory(resultDir);
	} 
	if (!harmonic) {
	    ClosedOrbitFinder<value_type, size_type,
	        boost::numeric::odeint::runge_kutta4<container_type> > cof(E_m, m_m, wo_m, N_m, accuracy,
									   maxitOrbit, Emin_m, Emax_m,
									   nSector_m, fieldmap_m, rguess,
									   type, scaleFactor_m, false);

	    // properties of one turn
	    container_type h_turn = cof.getInverseBendingRadius();
	    container_type r_turn = cof.getOrbit(angle);
	    container_type ds_turn = cof.getPathLength();
	    container_type fidx_turn = cof.getFieldIndex();
	    tunes = cof.getTunes();
	    ravg = cof.getAverageRadius();                   // average radius

	    container_type peo = cof.getMomentum(angle);

	    // write to terminal
	    *gmsg << "* ----------------------------" << endl
		  << "* Closed orbit info (Gordon units):" << endl
		  << "*" << endl
		  << "* average radius: " << ravg << " [m]" << endl
		  << "* initial radius: " << r_turn[0] << " [m]" << endl
		  << "* initial momentum: " << peo[0] << " [Beta Gamma]" << endl
		  << "* frequency error: " << cof.getFrequencyError() << endl
		  << "* horizontal tune: " << tunes.first << endl
		  << "* vertical tune: " << tunes.second << endl
		  << "* ----------------------------" << endl << endl;

	    // compute the number of steps per degree
	    value_type deg_step = N_m / 360.0;
	    // compute starting point of computation
	    size_type start = deg_step * angle;

	    // copy properties of the length of one sector (--> nStepsPerSector_m)
	    std::copy_n(r_turn.begin()+start,nStepsPerSector_m, r.begin());
	    std::copy_n(h_turn.begin()+start,nStepsPerSector_m, h.begin());
	    std::copy_n(fidx_turn.begin()+start,nStepsPerSector_m, fidx.begin());
	    std::copy_n(ds_turn.begin()+start,nStepsPerSector_m, ds.begin());

	    if(write_m){
	        writeOrbitOutput_m(tunes,ravg, cof.getFrequencyError(),r_turn,peo,h_turn,fidx_turn,ds_turn);
	    }
            
	} else {
	    *gmsg << "Not yet supported." << endl;
	    return false;
	}

	if(Options::cloTuneOnly)
	  throw OpalException("Do only CLO and tune calculation","... ");

	// computes the cyclotron maps for each angle and stores them into a vector
	std::vector<matrix_type> Mcycs(nStepsPerSector_m);
    
	// calculate only for a single sector (a nSector_-th) of the whole cyclotron
	for (size_type i = 0; i < nStepsPerSector_m; ++i) 
	    Mcycs[i] = mapgen.generateMap(H_m(h[i],h[i]*h[i]+fidx[i],-fidx[i]),ds[i],truncOrder_m);
    
	// search for a match distribution for each permutation 
	for(size_type permutation = 0; permutation < 1; permutation++){
	    findMatchDistribution(accuracy,h,fidx,ds,tunes,ravg,maxit,Mcycs,harmonic,2);
	}
	std::cout << "OPAL> * Number of matched distributions: " <<matchedSigmas_m.size() << std::endl;
	std::cout << "OPAL> * Last matched distribution information: "<<std::endl;
	
    }catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    //Returns true if at least one matched distribution was found.
    if ( 0 < matchedSigmas_m.size()) return (0 < matchedSigmas_m.size());
    //Returns an identity matrix, if  no matched distribution was found
    else{
        *gmsg << "No matched distribution found, returning unity matrix instead.";
        matrix_type One(6,6);
	for(size_type i = 0; i < 6; i++){
	  for(size_type j = 0; j < 6; j++) One(i,j) = 0;
	}
	for(size_type i = 0 ; i < 6; i++) One(i,i) = 1;
        sigma_m.swap(One);
	return true;
    }
}

template<typename Value_type, typename Size_type> bool SigmaGenerator<Value_type, Size_type>::InvertMatrix(const complex_matrix_type& R,complex_matrix_type& invR){
    //creates a working copy of R
    complex_matrix_type A(R);
    //permutation matrix for the LU-factorization
    pmatrix_type pm(A.size1());
    //LU-factorization
    int res = lu_factorize(A,pm);
    if( res != 0) return false;
    // create identity matrix of invR
    invR.assign(boost::numeric::ublas::identity_matrix<complex_number_type>(A.size1()));
    // backsubstitute to get the inverse
    boost::numeric::ublas::lu_substitute(A,pm,invR);

    return true;
}
template<typename Value_type, typename Size_type>
  typename SigmaGenerator<Value_type, Size_type>::value_type SigmaGenerator<Value_type, Size_type>::isEigenEllipse(const matrix_type& Mturn, const matrix_type& sigma) {
    // compute sigma matrix after one turn
    matrix_type newSigma = matt_boost::gemmm<matrix_type>(Mturn,sigma,boost::numeric::ublas::trans(Mturn));

    // return L2 error
    return L2ErrorNorm(sigma,newSigma);
}

template<typename Value_type, typename Size_type>
inline typename SigmaGenerator<Value_type, Size_type>::matrix_type& SigmaGenerator<Value_type, Size_type>::getSigma() {
    return sigma_m;
}

template<typename Value_type, typename Size_type>
inline typename SigmaGenerator<Value_type, Size_type>::size_type SigmaGenerator<Value_type, Size_type>::getIterations() const {
    return (converged_m) ? niterations_m : size_type(0);
}

template<typename Value_type, typename Size_type>
inline typename SigmaGenerator<Value_type, Size_type>::value_type SigmaGenerator<Value_type, Size_type>::getError() const {
    return error_m;
}

template<typename Value_type, typename Size_type>
inline std::array<Value_type,3> SigmaGenerator<Value_type, Size_type>::getEmittances() const {
    value_type bgam = gamma_m*beta_m;
    return std::array<value_type,3>{{
		    emittance_m[0]/Physics::pi/bgam, emittance_m[1]/Physics::pi/bgam, emittance_m[2]/Physics::pi/bgam}};
}

// -----------------------------------------------------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type>
void SigmaGenerator<Value_type, Size_type>::initialize(value_type nuz, value_type ravg) {
    /*
     * The initialization is based on the analytical solution of using a spherical symmetric beam in the paper:
     * Transverse-longitudinal coupling by space charge in cyclotrons
     * by Dr. Christian Baumgarten
     * (formulas: (46), (56), (57))
     */


    /* Units:
     * ----------------------------------------------
     * [wo] = 1/s
     * [nh] = 1
     * [q0] = e
     * [I] = A
     * [eps0] = (A*s)^{2}/(N*m^{2})
     * [E0] = MeV/(c^{2}) (with speed of light c)
     * [beta] = 1
     * [gamma] = 1
     * [m] = kg
     *
     * [lam] = m
     * [K3] = m
     * [alpha] = 10^{3}/(pi*mrad)
     * ----------------------------------------------
     */

    // helper constants
    value_type invbg = 1.0 / (beta_m * gamma_m);
    value_type micro = 1.0e-6;
    value_type mega = 1.0e6;

    // convert mass m_m from MeV/c^2 to eV*s^{2}/m^{2}
    value_type m = m_m * mega/(Physics::c * Physics::c);        // [m] = eV*s^{2}/m^{2}, [m_m] = MeV/c^2

    /* Emittance [ex] = [ey] = [ez] = mm mrad (emittance_m are normalized emittances
     * (i.e. emittance multiplied with beta*gamma)
     */
    value_type ex = emittance_m[0] * invbg;                        // [ex] = mm mrad
    value_type ey = emittance_m[1] * invbg;                        // [ey] = mm mrad
    value_type ez = emittance_m[2] * invbg;                        // [ez] = mm mrad

    // convert normalized emittances: mm mrad --> m rad (mm mrad: millimeter milliradian)
    ex *= micro;
    ey *= micro;
    ez *= micro;

    // initial guess of emittance, [e] = m rad
    value_type e = std::cbrt(ex * ey * ez);             // cbrt computes cubic root (C++11) <cmath>

    // cyclotron radius [rcyc] = m
    value_type rcyc = ravg / beta_m;

    // "average" inverse bending radius
    value_type h = 1.0 / ravg;            // [h] = 1/m

    // formula (57)
    value_type lam = 2.0 * Physics::pi * Physics::c / (wo_m * nh_m); // wavelength, [lam] = m
    value_type K3 = 3.0 * /* physics::q0 */ 1.0 * I_m * lam / (20.0 * std::sqrt(5.0) * Physics::pi * Physics::epsilon_0 * m *
                    Physics::c * Physics::c * Physics::c * beta_m * beta_m * gamma2_m * gamma_m);               // [K3] = m

    value_type alpha = /* physics::q0 */ 1.0 * Physics::mu_0 * I_m / (5.0 * std::sqrt(10.0) * m * Physics::c *
                       gamma_m * nh_m) * std::sqrt(rcyc * rcyc * rcyc / (e * e * e));                           // [alpha] = 1/rad --> [alpha] = 1

    value_type sig0 = std::sqrt(2.0 * rcyc * e) / gamma_m;                                                      // [sig0] = m*sqrt(rad) --> [sig0] = m

    // formula (56)
    value_type sig;                                     // [sig] = m
    if (alpha >= 2.5) {
        sig = sig0 * std::cbrt(1.0 + alpha);            // cbrt computes cubic root (C++11) <cmath>
    } else if (alpha >= 0) {
        sig = sig0 * (1 + alpha * (0.25 - 0.03125 * alpha));
    } else {
        throw OpalException("SigmaGenerator::initialize()", "Negative alpha value: " + std::to_string(alpha) + " < 0");
    }

    // K = Kx = Ky = Kz
    value_type K = K3 * gamma_m / (3.0 * sig * sig * sig);   // formula (46), [K] = 1/m^{2}
    value_type kx = h * h * gamma2_m;                        // formula (46) (assumption of an isochronous cyclotron), [kx] = 1/m^{2}

    value_type a = 0.5 * kx - K;    // formula (22) (with K = Kx = Kz), [a] = 1/m^{2}
    value_type b = K * K;           // formula (22) (with K = Kx = Kz and kx = h^2*gamma^2), [b] = 1/m^{4}


    // b must be positive, otherwise no real-valued frequency
    if (b < 0)
        throw OpalException("SigmaGenerator::initialize()", "Negative value --> No real-valued frequency.");

    value_type tmp = a * a - b;           // [tmp] = 1/m^{4}
    if (tmp < 0)
        throw OpalException("SigmaGenerator::initialize()", "Square root of negative number.");

    tmp = std::sqrt(tmp);               // [tmp] = 1/m^{2}

    if (a < tmp)
        throw OpalException("Error in SigmaGenerator::initialize()", "Square root of negative number.");

    if (h * h * nuz * nuz <= K)
        throw OpalException("SigmaGenerator::initialize()", "h^{2} * nu_{z}^{2} <= K (i.e. square root of negative number)");

    value_type Omega = std::sqrt(a + tmp);                // formula (22), [Omega] = 1/m
    value_type omega = std::sqrt(a - tmp);                // formula (22), [omega] = 1/m

    value_type A = h / (Omega * Omega + K);           // formula (26), [A] = m
    value_type B = h / (omega * omega + K);           // formula (26), [B] = m
    value_type invAB = 1.0 / (B - A);                 // [invAB] = 1/m

    // construct initial sigma-matrix (formula (29, 30, 31)
    // Remark: We multiply with 10^{6} (= mega) to convert emittances back.
    // 1 m^{2} = 10^{6} mm^{2}
    matrix_type sigma = boost::numeric::ublas::zero_matrix<value_type>(6);
    // formula (30), [sigma(0,0)] = m^2 rad * 10^{6} = mm^{2} rad = mm mrad
    sigma(0,0) = invAB * (B * ex / Omega + A * ez / omega) * mega;
    // [sigma(0,5)] = [sigma(5,0)] = m rad * 10^{6} = mm mrad	// 1000: m --> mm and 1000 to promille
    sigma(0,5) = sigma(5,0) = invAB * (ex / Omega + ez / omega) * mega;
    // [sigma(1,1)] = rad * 10^{6} = mrad (and promille)
    sigma(1,1) = invAB * (B * ex * Omega + A * ez * omega) * mega;
    // [sigma(1,4)] = [sigma(4,1)] = m rad * 10^{6} = mm mrad
    sigma(1,4) = sigma(4,1) = invAB * (ex * Omega+ez * omega) / (K * gamma2_m) * mega;
    // formula (31), [sigma(2,2)] = m rad * 10^{6} = mm mrad
    sigma(2,2) = ey / (std::sqrt(h * h * nuz * nuz - K)) * mega;
     
    sigma(3,3) = (ey * mega) * (ey * mega) / sigma(2,2);
    // [sigma(4,4)] = m^{2} rad * 10^{6} = mm^{2} rad = mm mrad
    sigma(4,4) = invAB * (A * ex * Omega + B * ez * omega) / (K * gamma2_m) * mega;
    // formula (30), [sigma(5,5)] = rad * 10^{6} = mrad (and promille)
    sigma(5,5) = invAB * (ex / (B * Omega) + ez / (A * omega)) * mega;

    // fill in initial guess of the sigma matrix (for each angle the same guess)
    sigmas_m.resize(nStepsPerSector_m);
    for (typename std::vector<matrix_type>::iterator it = sigmas_m.begin(); it != sigmas_m.end(); ++it) {
        *it = sigma;
    }

    if (write_m) {
        std::string energy = float2string(E_m);
        std::ofstream writeInit("data/maps/InitialSigmaPerAngleForEnergy"+energy+"MeV.dat",std::ios::app);
        writeInit << sigma << std::endl;
        writeInit.close();
    }
}

template<typename Value_type, typename Size_type> 
  bool SigmaGenerator<Value_type, Size_type>::findMatchDistribution(value_type accuracy,
								    container_type& h,
								    container_type& fidx,
								    container_type& ds,
								    std::pair<value_type,value_type> tunes,
								    value_type ravg,
								    size_type maxit,
								    std::vector<matrix_type>& Mcycs,
								    bool harmonic,
								    size_type permutation){
    // object for space charge map
    MapGenerator<value_type,size_type,Series,Map,Hamiltonian,SpaceCharge> mapgen(nStepsPerSector_m);
    value_type const_ds = 0.0;

    try {
        // stores computed space charge map for each angle
        std::vector<matrix_type> Mscs(nStepsPerSector_m);

	// (inverse) transformation matrix
        complex_matrix_type R(6,6), invR(6,6);

        // new initial sigma matrix
        matrix_type newSigma(6,6);
        value_type weight = 0.65;

	// for exiting the loop if maxit is reached
        bool stop = false;

	// initialize sigma matrices (for each angle one) (first guess)
        initialize(tunes.second,ravg);

        // calculate only for a single sector (a nSector_-th) of the whole cyclotron
        for (size_type i = 0; i < nStepsPerSector_m; ++i) {
            if (!harmonic) {
                Mscs[i]  = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),sigmas_m[i](2,2),sigmas_m[i](4,4)),ds[i],truncOrder_m);
             } else {
                Mscs[i] = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),sigmas_m[i](2,2),sigmas_m[i](4,4)),const_ds,truncOrder_m);
             }
         }
	 
	 // Twiss-Matrix (One turn matrix)
         mapgen.combine(Mscs,Mcycs);
         matrix_type Mturn = mapgen.getMap();

         while (error_m > accuracy && !stop) {
	     /*
	      * Version 2.0 decouple was substituded by setEigenVectors
	      * and updateInitialSigma was changed according to Eq. 18 Sec. 2.4 
	      * Frey Semester thesis
	      */
	     // Returns the Eigenvalues of the Twiss-Matrix,
	     // orders them and computes R and invR
	     setEigenVectors(Mturn,R,invR,permutation);
	 	
	     // computes the new sigma matrix
	     newSigma = updateInitialSigma(R,invR);
		
	     // computes new sigma matrices for all angles (except for initial sigma)
	     updateSigma(Mscs,Mcycs);
		
	     // computes error
	     error_m = L2ErrorNorm(sigmas_m[0],newSigma);
		
	     // write new initial sigma-matrix into vector
	     sigmas_m[0] = weight*newSigma + (1.0-weight)*sigmas_m[0];

	     // compute new space charge maps
	     for (size_type i = 0; i < nStepsPerSector_m; ++i) {
	         if (!harmonic) {
		     Mscs[i] = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),sigmas_m[i](2,2),sigmas_m[i](4,4)),ds[i],truncOrder_m);
		 } else {
		     Mscs[i] = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),sigmas_m[i](2,2),sigmas_m[i](4,4)),const_ds,truncOrder_m);
		 }
	     }   

	     // construct new one turn transfer matrix M
	     mapgen.combine(Mscs,Mcycs);
	     Mturn = mapgen.getMap();

	     // check if number of iterations has maxit exceeded.
	     stop = (niterations_m++ > maxit);

	     // check for convergence
	     converged_m = error_m < accuracy;
	     // store converged sigma-matrix
	     if(converged_m){
	         for(size_type i = 0; i < 6; i++){
	             for(size_type j = 0; j < 6; j++){
	                 if(std::abs(newSigma(i,j)) < 1e-12) newSigma(i,j) = 0.;
	             }
		 }
		 mapgen.combine(Mcycs);
		 matrix_type Mcyc = mapgen.getMap();
		 if(write_m) writeMatched(newSigma,permutation,Mcycs,Mscs,Mcyc);
		 sigma_m.swap(newSigma);
		 matchedSigmas_m.push_back(newSigma);
	     }
	 }
	 // reset for next permutation
	 error_m = 1;
         niterations_m = 0;
    } catch(const std::exception& e) {
        std::cout << "No convergence for this permutation"<<std::endl;
	std::cout << "Probably no invert Matrix found"<<std::endl;
        std::cerr << e.what() << std::endl;
    }

    return converged_m;
}

template<typename Value_type, typename Size_type> 
  bool SigmaGenerator<Value_type, Size_type>::setEigenVectors(const matrix_type& MTurn, complex_matrix_type& R, complex_matrix_type& invR, size_type perm ){
    
    gsl_vector_complex_type eigenValues;
    gsl_matrix_complex_type eigenVectors;
    gsl_eigen_nonsymmv_workspace_type w;
    gsl_matrix_type MTurngsl;
    gsl_complex zGsl;
    //Initialization of the gsl Environment variables 
    MTurngsl = gsl_matrix_alloc(matrixDimension,matrixDimension);
    eigenValues = gsl_vector_complex_alloc(matrixDimension);
    eigenVectors = gsl_matrix_complex_alloc(matrixDimension,matrixDimension);
    w = gsl_eigen_nonsymmv_alloc(matrixDimension);
    //turns the Twiss-Matrix into gsl-Matrix
    for(size_type i = 0; i < matrixDimension; i++){
        for(size_type j = 0; j < matrixDimension; j++){
            gsl_matrix_set(MTurngsl,i,j,MTurn(i,j));
        }
    }
    //calculates the eigenValues and eigenVectors of MTurngsl
    //and stores them in R
    gsl_eigen_nonsymmv(MTurngsl,eigenValues,eigenVectors,w);
    gsl_eigen_nonsymmv_free(w);
    //Sets the Eigenvectors in R
    for(size_type i = 0; i < matrixDimension; i++){
        gsl_vector_complex_view evec_i = gsl_matrix_complex_column(eigenVectors,i);
        for(size_type j = 0;j < matrixDimension; j++){
	    zGsl = gsl_vector_complex_get(&evec_i.vector,j);
	    complex_number_type z(GSL_REAL(zGsl),GSL_IMAG(zGsl));
	    R(i,j) = z;
	}
    }

    complex_vector_type eigen(6);
    for(size_type i = 0; i < matrixDimension; i++){
        zGsl = gsl_vector_complex_get(eigenValues,i);
        complex_number_type z(GSL_REAL(zGsl),GSL_IMAG(zGsl));
	eigen(i) = z;
    }
    gsl_vector_complex_free(eigenValues);
    gsl_matrix_complex_free(eigenVectors);
    gsl_matrix_free(MTurngsl);

    // Sorts the Eigenvectors such that
    // Re(eig_x) < Re(eig_y) < Re(eig_z)

    std::vector<size_type> permutation = permutations_m[perm];

    bool isZdirection = false;
    std::vector<complex_vector_type> zVectors;
    std::vector<complex_vector_type> xyVectors; 

    for(size_type i = 0; i < 6; i++){
        complex_number_type z = R(i,0);
	if(std::abs(z) < 1e-13) z = 0.;
	if(z == 0.) isZdirection = true;
	
	complex_vector_type v(6);
	if(isZdirection){
	    for(size_type j = 0;j < matrixDimension; j++){
	        complex_number_type z = R(i,j);
		v(j) = z;
	    }
	    zVectors.push_back(v);
	}
	else{
	    for(size_type j = 0;j < matrixDimension; j++){
	        complex_number_type z = R(i,j);
		v(j) = z;
	    }
	    xyVectors.push_back(v);
	}
        isZdirection = false;
    }
    // Norms the Eigenvectors
     for(size_type i = 0; i < 4; i++){
        value_type norm{0};
        for(size_type j = 0; j < 6; j++) norm += std::norm(xyVectors[i](j));
	for(size_type j = 0; j < 6; j++) xyVectors[i](j) /= std::sqrt(norm);
    }
    for(size_type i = 0; i < 2; i++){
        value_type norm{0};
        for(size_type j = 0; j < 6; j++) norm += std::norm(zVectors[i](j));
	for(size_type j = 0; j < 6; j++) zVectors[i](j) /= std::sqrt(norm);
    }
    
    for(value_type i = 0; i < 6; i++){
        R(i,0) = xyVectors[permutation[0]](i);
        R(i,1) = xyVectors[permutation[1]](i);
        R(i,2) = zVectors[0](i);
        R(i,3) = zVectors[1](i);
        R(i,4) = xyVectors[permutation[2]](i);
        R(i,5) = xyVectors[permutation[3]](i); 
    }

    //Sets the invert matrix of R in invR
    return InvertMatrix(R,invR);
}

template<typename Value_type, typename Size_type>
typename SigmaGenerator<Value_type, Size_type>::matrix_type SigmaGenerator<Value_type, Size_type>::updateInitialSigma(const complex_matrix_type& R, const complex_matrix_type& invR) {
    
    matrix_type newSigmareal(6,6);
    complex_matrix_type D(6,6),newSigma(6,6),S(6,6);
    
    // Initialization of D and S Matrix
    for(size_type i = 0; i < matrixDimension; i++){
        for(size_type j = 0; j < matrixDimension; j++){
	    S(i,j) = D(i,j) = 0.;
        }
    }
    // Build Skew-Matrix
    S(0,1) = S(2,3) = S(4,5) = 1.;
    S(1,0) = S(3,2) = S(5,4) = -1.;

    // Build new D-Matrix
    // Section 2.4 Eq. 18 in M. Frey Semester thesis
    // D = diag(i*emx,-i*emx,i*emy,-i*emy,i*emz, -i*emz) 
    
    value_type invbg = 1.0 / (beta_m * gamma_m);
    complex_number_type im(0,1);    

    for(size_type i = 0; i < 3; i++){
            D(2*i,2*i) = emittance_m[i]*invbg*im;
            D(2*i+1,2*i+1) = -emittance_m[i]*invbg*im;
    }

    // Computing of new Sigma
    // sigma = -R*D*R^{-1}*S

    newSigma = boost::numeric::ublas::prod(invR,S);
    newSigma = boost::numeric::ublas::prod(D,newSigma);
    newSigma = boost::numeric::ublas::prod(-R,newSigma);
    
    for(value_type i = 0; i < 6; i++){
        for(value_type j = 0; j < 6; j++){
	  newSigmareal(i,j) = newSigma(i,j).real(); 
        }
    }

    for(value_type i = 0; i < 6; i++){
        if(newSigmareal(i,i) < 0.) newSigmareal(i,i) *= -1.0;
    }
    
    return newSigmareal;
}
 
template<typename Value_type, typename Size_type>
void SigmaGenerator<Value_type, Size_type>::updateSigma(const std::vector<matrix_type>& Mscs, const std::vector<matrix_type>& Mcycs) {
    matrix_type M = boost::numeric::ublas::matrix<value_type>(6,6);

    std::ofstream writeSigma;

    if (write_m) {
        std::string energy = float2string(E_m);
        writeSigma.open("data/maps/SigmaPerAngleForEnergy"+energy+"MeV.dat",std::ios::app);
    }

    // initial sigma is already computed
    for (size_type i = 1; i < nStepsPerSector_m; ++i) {
        // transfer matrix for one angle
        M = boost::numeric::ublas::prod(Mscs[i - 1],Mcycs[i - 1]);
        // transfer the matrix sigma
        sigmas_m[i] = matt_boost::gemmm<matrix_type>(M,sigmas_m[i - 1],boost::numeric::ublas::trans(M));

        if (write_m)
            writeSigma << sigmas_m[i] << std::endl;
    }

    if (write_m) {
        writeSigma << std::endl;
        writeSigma.close();
    }
}

template<typename Value_type, typename Size_type>
typename SigmaGenerator<Value_type, Size_type>::value_type SigmaGenerator<Value_type, Size_type>::L2ErrorNorm(const matrix_type& oldS, const matrix_type& newS) {
    // compute difference
    matrix_type diff = newS - oldS;

    // sum squared error up and take square root
    return std::sqrt(std::inner_product(diff.data().begin(),diff.data().end(),diff.data().begin(),0.0));
}

template<typename Value_type, typename Size_type> void SigmaGenerator<Value_type, Size_type>::writeMatched(const matrix_type& match, const size_type permutation,const std::vector<matrix_type>& Mcycs, const std::vector<matrix_type>& Mscs, const matrix_type& Mcyc){
    
    std::string energy = float2string(E_m);
    std::ofstream writeSigmas("data/Matches"+energy+"MeV.txt", std::ios::app);
    writeSigmas<<std::left<<std::setprecision(5)
               <<std::setw(12)<<E_m
               <<std::setw(12)<<I_m
	       <<std::setw(12)<<permutation
               <<std::setw(12)<<niterations_m
               <<std::setw(12)<<emittance_m[0]
               <<std::setw(12)<<emittance_m[1]
	       <<std::endl;
    for(size_type i = 0; i < 6; i++){
        for(size_type j = 0; j < 6; j++){
	    writeSigmas<<std::left<<std::setprecision(5)
	               <<std::setw(12)<<match(i,j);
        }
        writeSigmas<<std::endl;
    }
    writeSigmas<<std::endl;
    writeSigmas.close();
    
    std::ofstream writeOneTurn("data/OneTurnMatrix"+energy+"MeV.txt", std::ios::app);
    for(size_type i = 0; i < 6; i++){
        for(size_type j = 0; j < 6; j++){
	    writeOneTurn<<std::left<<std::setprecision(5)
	               <<std::setw(12)<<Mcyc(i,j);
        }
        writeOneTurn<<std::endl;
    }
    writeOneTurn<<std::endl;
    writeOneTurn.close();

    std::ofstream writeSigma;
    std::stringstream ss;
    ss << permutation;
    std::string p = ss.str();
    writeSigma.open("data/SigmaPerAngleForEnergy"+p+"_"+energy+"MeV.dat",std::ios::app);
    matrix_type M = boost::numeric::ublas::matrix<value_type>(6,6);
    // initial sigma is already computed
    for (size_type i = 1; i < nStepsPerSector_m; ++i) {
        // transfer matrix for one angle
        M = boost::numeric::ublas::prod(Mscs[i - 1],Mcycs[i - 1]);
        // transfer the matrix sigma
        sigmas_m[i] = matt_boost::gemmm<matrix_type>(M,sigmas_m[i - 1],boost::numeric::ublas::trans(M));
	
	for(size_type j = 0; j < 6; j++){
	    for(size_type k = 0; k < 6; k++){
	      if(std::abs(sigmas_m[i](j,k)) < 10e-13) sigmas_m[i](j,k) = 0.; 
	        writeSigma<<std::left<<std::setprecision(5)
			   <<std::setw(12)<<sigmas_m[i](j,k);
	    }
	    writeSigma<<std::endl;
	}
    }
    writeSigma << std::endl;
    writeSigma.close();
}
 
template<typename Value_type, typename Size_type>
typename SigmaGenerator<Value_type, Size_type>::value_type
    SigmaGenerator<Value_type, Size_type>::L1ErrorNorm(const matrix_type& oldS,
                                                       const matrix_type& newS)
{
    // compute difference
    matrix_type diff = newS - oldS;
    
    std::for_each(diff.data().begin(), diff.data().end(),
                  [](value_type& val) {
                      return std::abs(val);
                  });

    // sum squared error up and take square root
    return std::accumulate(diff.data().begin(), diff.data().end(), 0.0);
}

template<typename Value_type, typename Size_type>
std::string SigmaGenerator<Value_type, Size_type>::float2string(value_type val) {
    std::ostringstream out;
    out << std::setprecision(6) << val;
    return out.str();
}

template<typename Value_type, typename Size_type>
void SigmaGenerator<Value_type, Size_type>::writeOrbitOutput_m(
    const std::pair<value_type,value_type>& tunes,
    const value_type& ravg,
    const value_type& freqError,
    const container_type& r_turn,
    const container_type& peo,
    const container_type& h_turn,
    const container_type&  fidx_turn,
    const container_type&  ds_turn)
{
    // write tunes
    std::ofstream writeTunes("data/Tunes.dat", std::ios::app);

    if(writeTunes.tellp() == 0) // if nothing yet written --> write description
        writeTunes << "energy [MeV]" << std::setw(15)
                   << "nur" << std::setw(25)
                   << "nuz" << std::endl;
    
    writeTunes << E_m << std::setw(30) << std::setprecision(10)
               << tunes.first << std::setw(25)
               << tunes.second << std::endl;
    
    // write average radius
    std::ofstream writeAvgRadius("data/AverageValues.dat", std::ios::app);
    
    if (writeAvgRadius.tellp() == 0) // if nothing yet written --> write description
        writeAvgRadius << "energy [MeV]" << std::setw(15)
                       << "avg. radius [m]" << std::setw(15)
                       << "r [m]" << std::setw(15)
                       << "pr [m]" << std::endl;
    
    writeAvgRadius << E_m << std::setw(25) << std::setprecision(10)
                   << ravg << std::setw(25) << std::setprecision(10)
                   << r_turn[0] << std::setw(25) << std::setprecision(10)
                   << peo[0] << std::endl;

    // write frequency error
    std::ofstream writePhase("data/FrequencyError.dat",std::ios::app);

    if(writePhase.tellp() == 0) // if nothing yet written --> write description
        writePhase << "energy [MeV]" << std::setw(15)
                   << "freq. error" << std::endl;
    
    writePhase << E_m << std::setw(30) << std::setprecision(10)
               << freqError << std::endl;

    // write other properties
    std::string energy = float2string(E_m);
    std::ofstream writeProperties("data/PropertiesForEnergy"+energy+"MeV.dat", std::ios::out);
    writeProperties << std::left
                    << std::setw(25) << "orbit radius"
                    << std::setw(25) << "orbit momentum"
                    << std::setw(25) << "inverse bending radius"
                    << std::setw(25) << "field index"
                    << std::setw(25) << "path length" << std::endl;
        
    for (size_type i = 0; i < r_turn.size(); ++i) {
        writeProperties << std::setprecision(10) << std::left
                        << std::setw(25) << r_turn[i]
                        << std::setw(25) << peo[i]
                        << std::setw(25) << h_turn[i]
                        << std::setw(25) << fidx_turn[i]
                        << std::setw(25) << ds_turn[i] << std::endl;
    }
        
    // close all files within this if-statement
    writeTunes.close();
    writeAvgRadius.close();
    writePhase.close();
    writeProperties.close();
}

#endif
