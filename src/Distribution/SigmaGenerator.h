/**
 * @file SigmaGenerator.h
 * The SigmaGenerator class uses the class <b>ClosedOrbitFinder</b> to get the parameters (inverse bending radius, path length
 * field index and tunes) to initialize the sigma matrix. It further has a private class member of type
 * <b>RDM</b> for decoupling. \n
 * The main function of this class is <b>match(value_type, size_type)</b>, where it iteratively tries to find a matched distribution for given
 * emittances, energy and current. The computation stops when the L2-norm is smaller than a user-defined tolerance. \n
 * In default mode it prints all space charge maps, cyclotron maps and second moment matrices. The orbit properties, i.e.
 * tunes, average radius, orbit radius, inverse bending radius, path length, field index and frequency error, are printed
 * as well.
 *
 * @author Matthias Frey
 * @version 1.0
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

#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include "Utilities/OpalException.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#if BOOST_VERSION >= 106000
#include <boost/numeric/odeint/integrate/check_adapter.hpp>
#endif

#include "rdm.h"

#include "matrix_vector_operation.h"
#include "ClosedOrbitFinder.h"
#include "MapGenerator.h"
#include "Harmonics.h"

#include <boost/numeric/ublas/io.hpp>

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
        /// Type for storing maps
        typedef boost::numeric::ublas::matrix<value_type> matrix_type;
        /// Type for storing the sparse maps
        typedef boost::numeric::ublas::compressed_matrix<value_type,boost::numeric::ublas::row_major> sparse_matrix_type;
        /// Type for storing vectors
        typedef boost::numeric::ublas::vector<value_type> vector_type;
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

        /// Searches for a matched distribution.
        /*!
         * Returns "true" if a matched distribution within the accuracy could be found, returns "false" otherwise.
         * @param accuracy is used for computing the equilibrium orbit (to a certain accuracy)
         * @param maxit is the maximum number of iterations (the program stops if within this value no stationary
         *        distribution was found)
         * @param maxitOrbit is the maximum number of iterations for finding closed orbit
         * @param angle defines the start of the sector (one can choose any angle between 0° and 360°)
	 * @param guess value of radius for closed orbit finder
         * @param harmonic is a boolean. If "true" the harmonics are used instead of the closed orbit finder.
         */
        bool match(value_type, size_type, size_type, value_type, value_type, bool);

        /// Block diagonalizes the symplex part of the one turn transfer matrix
        /** It computes the transformation matrix <b>R</b> and its inverse <b>invR</b>.
         * It returns a vector containing the four eigenvalues
         * (alpha,beta,gamma,delta) (see formula (38), paper: Geometrical method of decoupling)
         */
        /*!
         * @param Mturn is a 6x6 dimensional one turn transfer matrix
         * @param R is the 4x4 dimensional transformation matrix (gets computed)
         * @param invR is the 4x4 dimensional inverse transformation (gets computed)
         */
        vector_type decouple(const matrix_type&, sparse_matrix_type&, sparse_matrix_type&);

        /// Checks if the sigma-matrix is an eigenellipse and returns the L2 error.
        /*!
         * The idea of this function is taken from Dr. Christian Baumgarten's program.
         * @param Mturn is the one turn transfer matrix
         * @param sigma is the sigma matrix to be tested
         */
        value_type isEigenEllipse(const matrix_type&, const matrix_type&);

        /// Returns the converged stationary distribution
        matrix_type& getSigma();

        /// Returns the number of iterations needed for convergence (if not converged, it returns zero)
        size_type getIterations() const;

        /// Returns the error (if the program didn't converged, one can call this function to check the error)
        value_type getError() const;

        /// Returns the emittances (ex,ey,ez) in \f$ \pi\ mm\ mrad \f$ for which the code converged (since the whole simulation is done on normalized emittances)
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

        /// Vector storing the second moment matrix for each angle
        std::vector<matrix_type> sigmas_m;

        /// <b>RDM</b>-class member used for decoupling
        RDM<value_type, size_type> rdm_m;

        /// Stores the Hamiltonian of the cyclotron
        Hamiltonian H_m;

        /// Stores the Hamiltonian for the space charge
        SpaceCharge Hsc_m;

        /// All variables x, px, y, py, z, delta
        Series x_m, px_m, y_m, py_m, z_m, delta_m;

        /// Initializes a first guess of the sigma matrix with the assumption of a spherical symmetric beam (ex = ey = ez). For each angle split the same initial guess is taken.
        /*!
         * @param nuz is the vertical tune
         * @param ravg is the average radius of the closed orbit
         */
        void initialize(value_type, value_type);

        /// Reduces the 6x6 matrix to a 4x4 matrix (for more explanations consider the source code comments)
        template<class matrix>
            void reduce(matrix&);
        /// Expands the 4x4 matrix back to dimension 6x6 (for more explanations consider the source code comments)
        template<class matrix>
            void expand(matrix&);

        /// Computes the new initial sigma matrix
        /*!
         * @param M is the 6x6 one turn transfer matrix
         * @param eigen stores the 4 eigenvalues: alpha, beta, gamma, delta (in this order)
         * @param R is the transformation matrix
         * @param invR is the inverse transformation matrix
         */
        matrix_type updateInitialSigma(const matrix_type&, const vector_type&, sparse_matrix_type&, sparse_matrix_type&);

        /// Computes new sigma matrices (one for each angle)
        /*!
         * Mscs is a vector of all space charge maps
         * Mcycs is a vector of all cyclotron maps
         */
        void updateSigma(const std::vector<matrix_type>&, const std::vector<matrix_type>&);

        /// Returns the L2-error norm between the old and new sigma-matrix
        /*!
         * @param oldS is the old sigma matrix
         * @param newS is the new sigma matrix
         */
        value_type L2ErrorNorm(const matrix_type&, const matrix_type&);

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
  scaleFactor_m(scaleFactor), sigmas_m(nStepsPerSector_m)
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
  bool SigmaGenerator<Value_type, Size_type>::match(value_type accuracy, size_type maxit, size_type maxitOrbit, value_type angle, value_type rguess, bool harmonic) {
    /* compute the equilibrium orbit for energy E_
     * and get the the following properties:
     * - inverse bending radius h
     * - step sizes of path ds
     * - tune nuz
     */

    try {

        // object for space charge map and cyclotron map
        MapGenerator<value_type, size_type, Series, Map, Hamiltonian, SpaceCharge> mapgen(nStepsPerSector_m);

        // compute cyclotron map and space charge map for each angle and store them into a vector
        std::vector<matrix_type> Mcycs(nStepsPerSector_m), Mscs(nStepsPerSector_m);

        container_type h(nStepsPerSector_m), r(nStepsPerSector_m), ds(nStepsPerSector_m), fidx(nStepsPerSector_m);
        value_type ravg = 0.0, const_ds = 0.0;
        std::pair<value_type,value_type> tunes;

        if (!harmonic) {
            ClosedOrbitFinder<value_type, size_type,
                boost::numeric::odeint::runge_kutta4<container_type> > cof(E_m, m_m, wo_m, N_m, accuracy,
                                                                           maxitOrbit, Emin_m, Emax_m,
                                                                           nSector_m, fieldmap_m, rguess,
                                                                           scaleFactor_m, false);

            // properties of one turn
            container_type h_turn = cof.getInverseBendingRadius();
            container_type r_turn = cof.getOrbit(angle);
            container_type ds_turn = cof.getPathLength();
            container_type fidx_turn = cof.getFieldIndex();
            tunes = cof.getTunes();
            ravg = cof.getAverageRadius();                   // average radius

	    container_type peo = cof.getMomentum(angle);
            
            // write properties to file
            if (write_m)
                writeOrbitOutput_m(tunes, ravg, cof.getFrequencyError(),
                                   r_turn, peo, h_turn, fidx_turn, ds_turn);
            
//             // tune calculation from map
//             std::vector<matrix_type> map(nStepsPerSector_m);
//             
            // compute the number of steps per degree
            value_type deg_step = N_m / 360.0;
            // compute starting point of computation
            size_type start = deg_step * angle;

            // copy properties of the length of one sector (--> nStepsPerSector_m)
            std::copy_n(r_turn.begin()+start,nStepsPerSector_m, r.begin());
            std::copy_n(h_turn.begin()+start,nStepsPerSector_m, h.begin());
            std::copy_n(fidx_turn.begin()+start,nStepsPerSector_m, fidx.begin());
            std::copy_n(ds_turn.begin()+start,nStepsPerSector_m, ds.begin());
//             
//             for (std::size_t i = 0; i < map.size(); ++i)
//                 map[i] = mapgen.generateMap(H_m(h[i],
//                                                 h[i]*h[i]+fidx[i],
//                                                 -fidx[i]),
//                                             ds[i],
//                                             truncOrder_m);
//             
//             matrix_type mtot = mapgen.combine(map);
//             map.clear();
//             tunes = mapgen.computeTunes(mtot);
//             tunes.first = tunes.first * nSector_m;
//             tunes.second = tunes.second * nSector_m;


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

        } else {
            *gmsg << "Not yet supported." << endl;
            return false;
//             MagneticField<double> bField(fieldmap_m);
//             bField.read(scaleFactor_m);
            
//             int nradial = bField.getNradSteps();
//             int ntheta = bField.getNtetSteps();
            
//             Harmonics<value_type, size_type> H(wo_m,Emin_m, Emax_m, ntheta, nradial, nSector_m, E_m, m_m);
//             Mcycs = H.computeMap("data/inj2sym_mainharms.4","data/inj2sym_mainharms.8",4);
//             ravg = H.getRadius();
//             tunes = H.getTunes();
//             const_ds = H.getPathLength();
        }


        if(Options::cloTuneOnly)
            throw OpalException("Do only CLO and tune calculation","... ");

        // initialize sigma matrices (for each angle one) (first guess)
        initialize(tunes.second,ravg);

        // for writing
        std::ofstream writeMturn, writeMcyc, writeMsc;

        if (write_m) {

            std::string energy = float2string(E_m);

            writeMturn.open("data/maps/OneTurnMapForEnergy"+energy+"MeV.dat",std::ios::app);
            writeMsc.open("data/maps/SpaceChargeMapPerAngleForEnergy"+energy+"MeV.dat",std::ios::app);
            writeMcyc.open("data/maps/CyclotronMapPerAngleForEnergy"+energy+"MeV.dat",std::ios::app);

            writeMturn << "--------------------------------" << std::endl;
            writeMturn << "Iteration: 0 " << std::endl;
            writeMturn << "--------------------------------" << std::endl;

            writeMsc << "--------------------------------" << std::endl;
            writeMsc << "Iteration: 0 " << std::endl;
            writeMsc << "--------------------------------" << std::endl;
        }

        // calculate only for a single sector (a nSector_-th) of the whole cyclotron
        for (size_type i = 0; i < nStepsPerSector_m; ++i) {
            if (!harmonic) {
                Mcycs[i] = mapgen.generateMap(H_m(h[i],h[i]*h[i]+fidx[i],-fidx[i]),ds[i],truncOrder_m);
                Mscs[i] = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),sigmas_m[i](2,2),sigmas_m[i](4,4)),ds[i],truncOrder_m);
            } else {
                Mscs[i] = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),sigmas_m[i](2,2),sigmas_m[i](4,4)),const_ds,truncOrder_m);
            }

            if (write_m) {
                writeMcyc << Mcycs[i] << std::endl;
                writeMsc << Mscs[i] << std::endl;
            }
        }

        // one turn matrix
        mapgen.combine(Mscs,Mcycs);
        matrix_type Mturn = mapgen.getMap();

        if (write_m)
            writeMturn << Mturn << std::endl;

        // (inverse) transformation matrix
        sparse_matrix_type R, invR;

        // eigenvalues
        vector_type eigen(4);

        // new initial sigma matrix
        matrix_type newSigma(6,6);

        // for exiting loop
        bool stop = false;

        value_type weight = 0.05;

        while (error_m > accuracy && !stop) {
            // decouple transfer matrix and compute (inverse) tranformation matrix
            eigen = decouple(Mturn,R,invR);

            // construct new initial sigma-matrix
            newSigma = updateInitialSigma(Mturn,eigen,R,invR);

            // compute new sigma matrices for all angles (except for initial sigma)
            updateSigma(Mscs,Mcycs);

            // compute error
            error_m = L2ErrorNorm(sigmas_m[0],newSigma);

            // write new initial sigma-matrix into vector
            sigmas_m[0] = weight*newSigma + (1.0-weight)*sigmas_m[0];

            if (write_m) {
                writeMsc << "--------------------------------" << std::endl;
                writeMsc << "Iteration: " << niterations_m + 1 << std::endl;
                writeMsc << "--------------------------------" << std::endl;
            }

            // compute new space charge maps
            for (size_type i = 0; i < nStepsPerSector_m; ++i) {
                if (!harmonic) {
                    Mscs[i] = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),sigmas_m[i](2,2),sigmas_m[i](4,4)),ds[i],truncOrder_m);
                } else {
                    Mscs[i] = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),sigmas_m[i](2,2),sigmas_m[i](4,4)),const_ds,truncOrder_m);
                }

                if (write_m)
                    writeMsc << Mscs[i] << std::endl;
            }

            // construct new one turn transfer matrix M
            mapgen.combine(Mscs,Mcycs);
            Mturn = mapgen.getMap();

            if (write_m) {
                writeMturn << "--------------------------------" << std::endl;
                writeMturn << "Iteration: " << niterations_m + 1 << std::endl;
                writeMturn << "--------------------------------" << std::endl;
                writeMturn << Mturn << std::endl;
            }

            // check if number of iterations has maxit exceeded.
            stop = (niterations_m++ > maxit);
        }

        // store converged sigma-matrix
        sigma_m.resize(6,6,false);
        sigma_m.swap(newSigma);

        // returns if the sigma matrix has converged
        converged_m = error_m < accuracy;

        // Close files
        if (write_m) {
            writeMturn.close();
            writeMsc.close();
            writeMcyc.close();
        }

    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    
    if ( converged_m && write_m ) {
        // write tunes
        std::ofstream writeSigmaMatched("data/MatchedDistributions.dat", std::ios::app);
        
        std::array<double,3> emit = this->getEmittances();
        
        writeSigmaMatched << "* Converged (Ex, Ey, Ez) = (" << emit[0]
                          << ", " << emit[1] << ", " << emit[2]
                          << ") pi mm mrad for E= " << E_m << " (MeV)"
                          << std::endl << "* Sigma-Matrix " << std::endl;

        for(unsigned int i = 0; i < sigma_m.size1(); ++ i) {
            writeSigmaMatched << std::setprecision(4)  << std::setw(11)
                              << sigma_m(i,0);
            for(unsigned int j = 1; j < sigma_m.size2(); ++ j) {
                writeSigmaMatched << " & " <<  std::setprecision(4)
                                  << std::setw(11) << sigma_m(i,j);
            }
            writeSigmaMatched << " \\\\" << std::endl;
        }
        
        writeSigmaMatched << std::endl;
        writeSigmaMatched.close();
    }

    return converged_m;
}

template<typename Value_type, typename Size_type>
typename SigmaGenerator<Value_type, Size_type>::vector_type SigmaGenerator<Value_type, Size_type>::decouple(const matrix_type& Mturn, sparse_matrix_type& R,
        sparse_matrix_type& invR) {
    // copy one turn matrix
    matrix_type M(Mturn);

    // reduce 6x6 matrix to 4x4 matrix
    reduce<matrix_type>(M);

    // compute symplex part
    matrix_type Ms = rdm_m.symplex(M);

    // diagonalize and compute transformation matrices
    rdm_m.diagonalize(Ms,R,invR);

    /*
     * formula (38) in paper of Dr. Christian Baumgarten:
     * Geometrical method of decoupling
     *
     * 		[0, 	alpha, 	0, 	0;
     * F_{d} =	-beta, 	0, 	0, 	0;
     * 		0, 	0, 	0, 	gamma;
     * 		0, 	0, 	-delta,	0]
     *
     *
     */
    vector_type eigen(4);
    eigen(0) =   Ms(0,1);       // alpha
    eigen(1) = - Ms(1,0);       // beta
    eigen(2) =   Ms(2,3);       // gamma
    eigen(3) = - Ms(3,2);       // delta
    return eigen;
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
    return std::array<value_type,3>({emittance_m[0]/Physics::pi/bgam, emittance_m[1]/Physics::pi/bgam, emittance_m[2]/Physics::pi/bgam});
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
    value_type kilo = 1.0e3;

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
    sigma(0,0) = invAB * (B * ex / Omega + A * ez / omega) * mega;                      // formula (30), [sigma(0,0)] = m^2 rad * 10^{6} = mm^{2} rad = mm mrad
    sigma(0,5) = sigma(5,0) = invAB * (ex / Omega + ez / omega) * mega;                 // [sigma(0,5)] = [sigma(5,0)] = m rad * 10^{6} = mm mrad	// 1000: m --> mm and 1000 to promille
    sigma(1,1) = invAB * (B * ex * Omega + A * ez * omega) * mega;                      // [sigma(1,1)] = rad * 10^{6} = mrad (and promille)
    sigma(1,4) = sigma(4,1) = invAB * (ex * Omega+ez * omega) / (K * gamma2_m) * mega;  // [sigma(1,4)] = [sigma(4,1)] = m rad * 10^{6} = mm mrad
    sigma(2,2) = ey / (std::sqrt(h * h * nuz * nuz - K)) * mega;                        // formula (31), [sigma(2,2)] = m rad * 10^{6} = mm mrad
    sigma(3,3) = (ey * mega) * (ey * mega) / sigma(2,2);
    sigma(4,4) = invAB * (A * ex * Omega + B * ez * omega) / (K * gamma2_m) * mega;     // [sigma(4,4)] = m^{2} rad * 10^{6} = mm^{2} rad = mm mrad
    sigma(5,5) = invAB * (ex / (B * Omega) + ez / (A * omega)) * mega;                  // formula (30), [sigma(5,5)] = rad * 10^{6} = mrad (and promille)

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
template<class matrix>
void SigmaGenerator<Value_type, Size_type>::reduce(matrix& M) {
    /* The 6x6 matrix gets reduced to a 4x4 matrix in the following way:
     *
     * a11 a12 a13 a14 a15 a16
     * a21 a22 a23 a24 a25 a26          a11 a12 a15 a16
     * a31 a32 a33 a34 a35 a36  -->     a21 a22 a25 a26
     * a41 a42 a43 a44 a45 a46          a51 a52 a55 a56
     * a51 a52 a53 a54 a55 a56          a61 a62 a65 a66
     * a61 a62 a63 a64 a65 a66
     */

    // copy x- and z-direction to a 4x4 matrix_type
    matrix_type M4x4(4,4);
    for (size_type i = 0; i < 2; ++i) {
        // upper left 2x2 [a11,a12;a21,a22]
        M4x4(i,0) = M(i,0);
        M4x4(i,1) = M(i,1);
        // lower left 2x2 [a51,a52;a61,a62]
        M4x4(i + 2,0) = M(i + 4,0);
        M4x4(i + 2,1) = M(i + 4,1);
        // upper right 2x2 [a15,a16;a25,a26]
        M4x4(i,2) = M(i,4);
        M4x4(i,3) = M(i,5);
        // lower right 2x2 [a55,a56;a65,a66]
        M4x4(i + 2,2) = M(i + 4,4);
        M4x4(i + 2,3) = M(i + 4,5);
    }

    M.resize(4,4,false);
    M.swap(M4x4);
}

template<typename Value_type, typename Size_type>
template<class matrix>
void SigmaGenerator<Value_type, Size_type>::expand(matrix& M) {
    /* The 4x4 matrix gets expanded to a 6x6 matrix in the following way:
     *
     *                          a11 a12 0 0 a13 a14
     * a11 a12 a13 a14          a21 a22 0 0 a23 a24
     * a21 a22 a23 a24  -->     0   0   1 0 0   0
     * a31 a32 a33 a34          0   0   0 1 0   0
     * a41 a42 a43 a44          a31 a32 0 0 a33 a34
     *                          a41 a42 0 0 a43 a44
     */

    matrix M6x6 = boost::numeric::ublas::identity_matrix<value_type>(6,6);

    for (size_type i = 0; i < 2; ++i) {
        // upper left 2x2 [a11,a12;a21,a22]
        M6x6(i,0) = M(i,0);
        M6x6(i,1) = M(i,1);
        // lower left 2x2 [a31,a32;a41,a42]
        M6x6(i + 4,0) = M(i + 2,0);
        M6x6(i + 4,1) = M(i + 2,1);
        // upper right 2x2 [a13,a14;a23,a24]
        M6x6(i,4) = M(i,2);
        M6x6(i,5) = M(i,3);
        // lower right 2x2 [a22,a34;a43,a44]
        M6x6(i + 4,4) = M(i + 2,2);
        M6x6(i + 4,5) = M(i + 2,3);
    }

    // exchange
    M.swap(M6x6);
}

template<typename Value_type, typename Size_type>
typename SigmaGenerator<Value_type, Size_type>::matrix_type SigmaGenerator<Value_type, Size_type>::updateInitialSigma(const matrix_type& M, const vector_type& eigen,
        sparse_matrix_type& R, sparse_matrix_type& invR) {

    /*
     * This function is based on the paper of Dr. Christian Baumgarten:
     * Transverse-Longitudinal Coupling by Space Charge in Cyclotrons (2012)
     */

    /*
     * Function input:
     * - M: one turn transfer matrix
     * - eigen = {alpha, beta, gamma, delta}
     * - R: transformation matrix (in paper: E)
     * - invR: inverse transformation matrix (in paper: E^{-1}
     */

    /* formula (18):
     * sigma = -E*D*E^{-1}*S
     * with diagonal matrix D (stores eigenvalues of sigma*S (emittances apart from +- i),
     * skew-symmetric matrix (formula (13)), and tranformation matrices E, E^{-1}
     */

    // normalize emittances
    value_type invbg = 1.0 / (beta_m * gamma_m);
    value_type ex = emittance_m[0] * invbg;
    value_type ey = emittance_m[1] * invbg;
    value_type ez = emittance_m[2] * invbg;


    // alpha^2-beta*gamma = 1

    /* 0        eigen(0) 0        0
     * eigen(1) 0        0        0
     * 0        0        0        eigen(2)
     * 0        0        eigen(3) 0
     *
     * M = cos(mux)*[1, 0; 0, 1] + sin(mux)*[alpha, beta; -gamma, -alpha], Book, p. 242
     *
     * -----------------------------------------------------------------------------------
     * X-DIRECTION and Z-DIRECTION
     * -----------------------------------------------------------------------------------
     * --> eigen(0) = sin(mux)*betax
     * --> eigen(1) = -gammax*sin(mux)
     *
     * What is sin(mux)?   --> alphax = 0 --> -alphax^2+betax*gammax = betax*gammax = 1
     *
     * Convention: betax > 0
     *
     * 1) betax = 1/gammax
     * 2) eig0 = sin(mux)*betax
     * 3) eig1 = -gammax*sin(mux)
     *
     * eig0 = sin(mux)/gammax
     * eig1 = -gammax*sin(mux) <--> 1/gammax = -sin(mux)/eig1
     *
     * eig0 = -sin(mux)^2/eig1 --> -sin(mux)^2 = eig0*eig1      --> sin(mux) = sqrt(-eig0*eig1)
     *                                                          --> gammax = -eig1/sin(mux)
     *                                                          --> betax = eig0/sin(mux)
     */


    // x-direction
    value_type alphax = 0.0;
    value_type betax  = std::sqrt(std::fabs(eigen(0) / eigen(1)));
    value_type gammax = 1.0 / betax;

    // z-direction
    value_type alphaz = 0.0;
    value_type betaz  = std::sqrt(std::fabs(eigen(2) / eigen(3)));
    value_type gammaz = 1.0 / betaz;

    /*
     * -----------------------------------------------------------------------------------
     * Y-DIRECTION
     * -----------------------------------------------------------------------------------
     *
     * m22 m23
     * m32 m33
     *
     * m22 = cos(muy) + alpha*sin(muy)
     * m33 = cos(muy) - alpha*sin(muy)
     *
     * --> cos(muy) = 0.5*(m22 + m33)
     *     sin(muy) = sign(m32)*sqrt(1-cos(muy)^2)
     *
     * m22-m33 = 2*alpha*sin(muy) --> alpha = 0.5*(m22-m33)/sin(muy)
     *
     * m23 = betay*sin(muy)     --> betay = m23/sin(muy)
     * m32 = -gammay*sin(muy)   --> gammay = -m32/sin(muy)
     */

    value_type cosy = 0.5 * (M(2,2) + M(3,3));

    value_type sign = (std::signbit(M(2,3))) ? value_type(-1) : value_type(1);

    value_type invsiny = sign / std::sqrt(std::fabs( 1.0 - cosy * cosy));

    value_type alphay = 0.5 * (M(2,2) - M(3,3)) * invsiny;
    value_type betay  =   M(2,3) * invsiny;
    value_type gammay = - M(3,2) * invsiny;

    // Convention beta>0
    if (std::signbit(betay))    // singbit = true if beta<0, else false
        betay  *= -1.0;

    // diagonal matrix with eigenvalues
    matrix_type D = boost::numeric::ublas::zero_matrix<value_type>(6,6);
    // x-direction
    D(0,1) =   betax  * ex;
    D(1,0) = - gammax * ex;
    // y-direction
    D(2,2) =   alphay * ey;
    D(3,3) = - alphay * ey;
    D(2,3) =   betay  * ey;
    D(3,2) = - gammay * ey;
    // z-direction
    D(4,5) =   betaz  * ez;
    D(5,4) = - gammaz * ez;

    // expand 4x4 transformation matrices to 6x6
    expand<sparse_matrix_type>(R);
    expand<sparse_matrix_type>(invR);

    // symplectic matrix
    sparse_matrix_type S(6,6,6);
    S(0,1) = S(2,3) = S(4,5) = 1;
    S(1,0) = S(3,2) = S(5,4) = -1;

    // --> get D (formula (18), paper: Transverse-Longitudinal Coupling by Space Charge in Cyclotrons
    // sigma = -R*D*R^{-1}*S
    matrix_type sigma = matt_boost::gemmm<matrix_type>(-invR,D,R);
    sigma = boost::numeric::ublas::prod(sigma,S);

    if (write_m) {
        std::string energy = float2string(E_m);
        std::ofstream writeSigma("data/maps/SigmaPerAngleForEnergy"+energy+"MeV.dat",std::ios::app);

        writeSigma << "--------------------------------" << std::endl;
        writeSigma << "Iteration: " << niterations_m + 1 << std::endl;
        writeSigma << "--------------------------------" << std::endl;

        writeSigma << sigma << std::endl;
        writeSigma.close();
    }

    return sigma;
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