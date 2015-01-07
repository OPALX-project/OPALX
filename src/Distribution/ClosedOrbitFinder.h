/**
 * @file ClosedOrbitFinder.h
 * The algorithm is based on the paper of M. M. Gordon: "Computation of closed orbits and basic focusing properties for
 * sector-focused cyclotrons and the design of 'cyclops'" (1983)
 * As template arguments one chooses the type of the variables and the integrator for the ODEs. The supported steppers can be found on
 * http://www.boost.org/ where it is part of the library Odeint.
 *
 * @author Matthias Frey
 * @version 1.0
 */

#ifndef CLOSEDORBITFINDER_H
#define CLOSEDORBITFINDER_H

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "error.h"
#include "physics.h"
#include "physical_error.h"

#include "MagneticField.h"


#include <fstream>

// include headers for integration
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>

/// @brief Finds a closed orbit of a cyclotron for a given energy
template<typename Value_type, typename Size_type, class Stepper>
class ClosedOrbitFinder
{
    public:
        /// Type of variables
        typedef Value_type value_type;
        /// Type for specifying sizes
        typedef Size_type size_type;
        /// Type of container for storing quantities (path length, orbit, etc.)
        typedef std::vector<value_type> container_type;
        /// Type for holding state of ODE values
        typedef std::vector<value_type> state_type;

        /// Sets the initial values for the integration and calls findOrbit().
        /*!
         * @param E is the energy [MeV] to which the closed orbit should be found
         * @param wo is the nominal orbital frequency (see paper of Dr. C. Baumgarten: "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons" (2012), formula (1))
         * @param N specifies the number of splits (2pi/N)
         * @param accuracy specifies the accuracy of the closed orbit
         * @param maxit is the maximal number of iterations done. Program stops if closed orbit not found within this time.
         * @param Emin is the minimum energy [MeV] needed in cyclotron
         * @param Emax is the maximum energy [MeV] reached in cyclotron
         * @param fieldmap is the location of the file that specifies the magnetic field
         */
        ClosedOrbitFinder(value_type, value_type, size_type, value_type, size_type, value_type, value_type, const std::string&);

        /// Returns the inverse bending radius (size of container N+1)
        container_type& getInverseBendingRadius();

        /// Returns the step lengths of the path (size of container N+1)
        container_type& getPathLength();

        /// Returns the field index (size of container N+1)
        container_type& getFieldIndex();

        /// Returns the radial and vertical tunes (in that order)
        std::pair<value_type,value_type> getTunes();

        /// Returns the closed orbit (size of container N+1)
        container_type& getOrbit();

        /// Returns the relativistic factor gamma
        value_type getGamma();

        /// Returns the average orbit radius
        value_type getAverageRadius();

        /// Returns the phase
        value_type getPhase();

        /// Returns true if a closed orbit could be found
        bool isConverged();

    private:
        /// Computes the closed orbit
        /*!
         * @param accuracy specifies the accuracy of the closed orbit
         * @param maxit is the maximal number of iterations done for finding the closed orbit
         */
        bool findOrbit(value_type, size_type);

        /// Fills in the values of h_m, ds_m, fidx_m. It gets called by in by constructor.
        void computeOrbitProperties();

        /// This function is called by the function getTunes().
        /*! Transfer matrix Y = [y11, y12; y21, y22] (see Gordon paper for more details).
         * @param y are the positions (elements y11 and y12 of Y)
         * @param py2 is the momentum of the second solution (element y22 of Y)
         * @param ncross is the number of sign changes (\#crossings of zero-line)
         */
        value_type computeTune(const std::array<value_type,2>&, value_type, size_type);

        value_type christian_computeTune(const std::array<value_type,2>&, value_type, size_type);

        /// This function computes nzcross_ which is used to compute the tune in z-direction and the phase error
        void computeVerticalOscillations();

        /// Stores current position in horizontal direction for the solutions of the ODE with different initial values
        std::array<value_type,2> x_m; // x_m = [x1, x2]
        /// Stores current momenta in horizontal direction for the solutions of the ODE with different initial values
        std::array<value_type,2> px_m; // px_m = [px1, px2]
        /// Stores current position in longitudinal direction for the solutions of the ODE with different initial values
        std::array<value_type,2> z_m; // z_m = [z1, z2]
        /// Stores current momenta in longitudinal direction for the solutions of the ODE with different initial values
        std::array<value_type,2> pz_m; // pz_m = [pz1, pz2]

        /// Stores the inverse bending radius
        container_type h_m;
        /// Stores the step length
        container_type ds_m;
        /// Stores the radial orbit (size: N_m+1)
        container_type r_m;
        /// Stores the radial momentum
        container_type pr_m;
        /// Stores the field index
        container_type fidx_m;

        /// Counts the number of zero-line crossings in horizontal direction (used for computing horizontal tune)
        size_type nxcross_m;
        /// Counts the number of zero-line crossings in vertical direction (used for computing vertical tune)
        size_type nzcross_m; //#crossings of zero-line in x- and z-direction

        /// Is the energy for which the closed orbit should be found
        value_type E_m;
        /// Is the nominal orbital frequency
        value_type wo_m;
        /// Is the number of angle splits (2*pi/N) and the number of integration steps
        size_type N_m;
        /// Is the angle step size
        value_type dtheta_m;

        /// Is the relativistic factor
        value_type gamma_m;

        /// Is the average radius
        value_type ravg_m;

        /// Is the phase
        value_type phase_m;

        /// Boolean which tells if a closed orbit for this configuration could be found (get set by the function findOrbit)
        bool converged_m;

        /// Minimum energy needed in cyclotron
        value_type Emin_m;

        /// Maximum energy reached in cyclotron
        value_type Emax_m;

        /// Stores the last orbit value (since we have to return to the beginning to check the convergence in the findOrbit() function. This last value is then deleted from the array but is stored in lastOrbitVal_m to compute the vertical oscillations)
        value_type lastOrbitVal_m;

        /// Stores the last momentum value (since we have to return to the beginning to check the convergence in the findOrbit() function. This last value is then deleted from the array but is stored in lastMomentumVal_m to compute the vertical oscillations)
        value_type lastMomentumVal_m;
        
        /// Boolean which is true if computeVerticalOscillations() executed, otherwise false. This is used for checking in getTunes() and getPhase().
        bool vertOscDone_m;

        /// Location of magnetic field
        std::string fieldmap_m;

	/// some proerties of the magnetic field map
	int nr_m, nth_m, nsc_m;
	double rmin_m, dr_m, dth_m;

        /// Defines the stepper for integration of the ODE's
        Stepper stepper_m;

        /// ONLY FOR STAND-ALONE PROGRAM
        float** bmag_m;

        /// Christian's function for finding equilibrium orbit
        void christian_findOrbit(value_type, size_type);
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

    template<typename Value_type, typename Size_type, class Stepper>
ClosedOrbitFinder<Value_type, Size_type, Stepper>::ClosedOrbitFinder(value_type E, value_type wo, size_type N, value_type accuracy,
        size_type maxit,value_type Emin, value_type Emax, const std::string& fieldmap)
    : nxcross_m(0), nzcross_m(0), E_m(E), wo_m(wo), N_m(N), dtheta_m(2.0*M_PI/double(N)),
    gamma_m(E/physics::E0+1.0), ravg_m(0), phase_m(0), converged_m(false), Emin_m(Emin), Emax_m(Emax),
    lastOrbitVal_m(0.0), lastMomentumVal_m(0.0), vertOscDone_m(false), fieldmap_m(fieldmap), stepper_m()
{

    if (Emin_m > Emax_m || E_m < Emin_m || E > Emax_m)
        Error::message("ClosedOrbitFinder",Error::invalid);
    
    // reserve storage for the orbit and momentum (--> size = 0, capacity = N_m+1)
    /*
     * we need N+1 storage, since dtheta = 2pi/N (and not 2pi/(N-1)) that's why we need N+1 integration steps
     * to return to the origin (but the return size is N_m)
     */
    r_m.reserve(N + 1);
    pr_m.reserve(N + 1);
    
    // reserve memory of N_m for the properties (--> size = 0, capacity = N_m)
    h_m.reserve(N);
    ds_m.reserve(N);
    fidx_m.reserve(N);
    
    // compute closed orbit
    converged_m = findOrbit(accuracy, maxit);
    
    // compute h, ds, fidx, rav (average radius)
    computeOrbitProperties();
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type& ClosedOrbitFinder<Value_type, Size_type, Stepper>::getInverseBendingRadius() {
    return h_m;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type& ClosedOrbitFinder<Value_type, Size_type, Stepper>::getPathLength() {
    return ds_m;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type& ClosedOrbitFinder<Value_type, Size_type, Stepper>::getFieldIndex() {
    return fidx_m;
}

template<typename Value_type, typename Size_type, class Stepper>
std::pair<Value_type,Value_type> ClosedOrbitFinder<Value_type, Size_type, Stepper>::getTunes() {
    // compute radial tune
    value_type nur = computeTune(x_m,px_m[1],nxcross_m);
    
    // compute nzcross_m
    if (!vertOscDone_m)
        computeVerticalOscillations();
    
    // compute vertical tune
    value_type nuz = computeTune(z_m,pz_m[1],nzcross_m);

    return std::make_pair(nur,nuz);
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type& ClosedOrbitFinder<Value_type, Size_type, Stepper>::getOrbit() {
    return r_m;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::getGamma() {
    return gamma_m;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::getAverageRadius() {
    return ravg_m;
}

template<typename Value_type, typename Size_type, class Stepper>
typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::getPhase() {
    
    // if the vertical oscillations aren't computed, we have to, since there we also compuote the frequency error.
    if(!vertOscDone_m)
        computeVerticalOscillations();
    
    return phase_m;
}

template<typename Value_type, typename Size_type, class Stepper>
inline bool ClosedOrbitFinder<Value_type, Size_type, Stepper>::isConverged() {
    return converged_m;
} 

// -----------------------------------------------------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type, class Stepper>
bool ClosedOrbitFinder<Value_type, Size_type, Stepper>::findOrbit(value_type accuracy, size_type maxit) {
    /* REMARK TO GORDON
     * q' = 1/b = 1/bcon
     * a' = a = acon
     */

    // ada Injector 2  int nsc = 4, nr = 179, Nth = 144, nth = 144 / 4; value_type r0 = 0.24, dr = 0.02;

    MagneticField::ReadHeader(&nr_m, &nth_m, &rmin_m, &dr_m, &dth_m, &nsc_m, fieldmap_m);
    INFOMSG("Magnetic fiel map properties nr= " << nr_m << " nth= " << nth_m << " rmin= " << rmin_m << " dr= " << dr_m << " dth= " << dth_m << " nsc= " << nsc_m << endl);

    bmag_m = MagneticField::malloc2df(nth_m,nr_m);
    MagneticField::ReadSectorMap(bmag_m,nr_m,nth_m,1,fieldmap_m,0.0);

    MagneticField::MakeNFoldSymmetric(bmag_m,nth_m,nr_m,nth_m,nsc_m);
    value_type bint, brint, btint;

    // velocity: beta = v/c = sqrt(1-1/(gamma*gamma))
    if (gamma_m == 0)
        PhysicalError::message("ClosedOrbitFinder::findOrbit",PhysicalError::undefined);
    
    // resize vectors (--> size = N_m+1, capacity = N_m+1), note: we do N_m+1 integration steps
    r_m.resize(N_m+1);
    pr_m.resize(N_m+1);
    
    // store acon and bcon locally
    value_type acon = physics::acon(wo_m);               // [acon] = m
    value_type invbcon = 1.0 / physics::bcon(wo_m);        // [bcon] = MeV*s/(C*m^2) = 10^6 T = 10^7 kG (kilo Gauss)

    // helper constants
    value_type p2;                                      // p^2 = p*p
    value_type pr2;                                     // squared radial momentum (pr^2 = pr*pr)
    value_type ptheta, invptheta;                       // Gordon, formula (5c)
    value_type invdenom;                                // denominator for computing dr,dpr
    value_type xold = 0.0;                              // for counting nxcross

    // index for reaching next element of the arrays r and pr (no nicer way found yet)
    size_type idx = 0;
    // observer for storing the current value after each ODE step (e.g. Runge-Kutta step) into the containers of r and pr
    auto store = [&](state_type& y, const value_type t)
    {        
        r_m[idx] = y[0];
        pr_m[idx] = y[1];

        // count number of crossings (excluding starting point --> idx>0)
        nxcross_m += (idx > 0) * (y[4] * xold < 0);
        xold = y[4];
        
        ++idx;
    };

    // define the six ODEs (using lambda function)
    std::function<void(const state_type&, state_type&, const double)> orbit_integration = [&](const state_type &y, state_type &dydt, const double theta){
        pr2 = y[1] * y[1];
        if (p2 < pr2)
            PhysicalError::message("ClosedOrbitFinder::findOrbit",PhysicalError::negative);
            
        // Gordon, formula (5c)
        ptheta = std::sqrt(p2 - pr2);
        invptheta = 1.0 / ptheta;

        // intepolate values of magnetic field
        MagneticField::interpolate(&bint,&brint,&btint,theta * 180 / M_PI,nr_m,nth_m,y[0],rmin_m,dr_m,bmag_m);
        bint *= invbcon;
        brint *= invbcon;

        // Gordon, formula (5a)
        dydt[0] = y[0] * y[1] * invptheta;
        // Gordon, formula (5b)
        dydt[1] = ptheta - y[0] * bint;
        // Gordon, formulas (9a) and (9b)
        for (size_type i = 2; i < 5; i += 2) {
            dydt[i] = (y[1] * y[i] + y[0] * p2 * y[i+1] * invptheta * invptheta) * invptheta;
            dydt[i+1] = - y[1] * y[i+1] * invptheta - (bint + y[0] * brint) * y[i];
        }
    };

    // define initial state container for integration: y = {r, pr, x1, px1, x2, px2}
    state_type y(6);

    container_type err(2);                                      // difference of last and first value of r (1. element) and pr (2. element)
    container_type delta = {0.0, 0.0};                          // correction term for initial values: r = r + dr, pr = pr + dpr; Gordon, formula (17)
    value_type error = std::numeric_limits<value_type>::max();  // amplitude of error; Gordon, formula (18) (a = a')
    size_type niterations = 0;                                  // if niterations > maxit --> stop iteration

    /*
     * Christian:
     * N = 1440 ---> N = 720 ---> dtheta = 2PI/720 --> nsteps = 721
     * 
     * 0, 2, 4, ... ---> jeden zweiten berechnene: 1, 3, 5, ... interpolieren --> 1440 Werte
     * 
     * Matthias:
     * N = 1440 --> dtheta = 2PI/1440 --> nsteps = 1441
     * 
     * 0, 1, 2, 3, 4, 5, ... --> 1440 Werte
     * 
     */
    
    // iterate until suggested energy (start with minimum energy)
    value_type E = Emin_m;

    // energy dependent values
    value_type en = E / physics::E0;                      // en = E/E0 = E/(mc^2) (E0 is kinetic energy)
    value_type p = acon * std::sqrt(en * (2.0 + en));     // momentum [p] = m; Gordon, formula (3)
    value_type gamma2 = (1.0 + en) * (1.0 + en);          // = gamma^2
    value_type beta = std::sqrt(1.0 - 1.0 / gamma2);
    p2 = p * p;                                           // p^2 = p*p
    value_type invgamma4 = 1.0 / (gamma2 * gamma2);       // = 1/gamma^4

    // set initial values for radius and radial momentum for lowest energy Emin
    // orbit, [r] = m; Gordon, formula (20)
    // radial momentum; Gordon, formula (20)
    container_type init = {beta * acon, 0.0};

    // store initial values for updating values for higher energies
    container_type previous_init = {0.0, 0.0};
    
    while (E <= E_m) {

        // (re-)set inital values for r and pr
        r_m[0] = init[0]; 
        pr_m[0] = init[1];
        
        // integrate until error smaller than user-define accuracy
        do {
            // (re-)set inital values
            x_m[0]  = 1.0;               // x1; Gordon, formula (10)
            px_m[0] = 0.0;               // px1; Gordon, formula (10)
            x_m[1]  = 0.0;               // x2; Gordon, formula (10)
            px_m[1] = 1.0;               // px2; Gordon, formula (10)
            nxcross_m = 0;               // counts the number of crossings of x-axis (excluding first step)
            idx = 0;                     // index for looping over r and pr arrays

            // fill container with initial states
            y = {init[0],init[1], x_m[0], px_m[0], x_m[1], px_m[1]};
            
            // integrate from 0 to 2*pi (one has to get back to the "origin")
            boost::numeric::odeint::integrate_n_steps(stepper_m,orbit_integration,y,0.0,dtheta_m,N_m,store);
            
            // write new state
            x_m[0] = y[2];
            px_m[0] = y[3];
            x_m[1] = y[4];
            px_m[1] = y[5];
            
            // compute error (compare values of orbit and momentum for theta = 0 and theta = 2*pi)
            // (Note: size = N_m+1 --> last entry is N_m)
            err[0] = r_m[N_m] - r_m[0];      // Gordon, formula (14)
            err[1] = pr_m[N_m] - pr_m[0];    // Gordon, formula (14)
            
            // correct inital values of r and pr
            invdenom = 1.0 / (x_m[0] + px_m[1] - 2.0);
            delta[0] = ((px_m[1] - 1.0) * err[0] - x_m[1] * err[1]) * invdenom; // dr; Gordon, formula (16a)
            delta[1] = ((x_m[0] - 1.0) * err[1] - px_m[0] * err[0]) * invdenom; // dpr; Gordon, formula (16b)
            
            // improved initial values; Gordon, formula (17) (here it's used for higher energies)
            init[0] += delta[0];
            init[1] += delta[1];
            
            // compute amplitude of the error
            error = std::sqrt(delta[0] * delta[0] + delta[1] * delta[1] * invgamma4) / r_m[0];
            
        } while (error > accuracy && niterations++ < maxit);
        
        // reset iteration counter
        niterations = 0;
        
        // reset correction term
        delta[0] = delta[1] = 0.0;
        
        // increase energy by one
        E += 1.0;
        
        // set constants for new energy E
        en = E / physics::E0;                     // en = E/E0 = E/(mc^2) (E0 is kinetic energy)
        p = acon * std::sqrt(en * (2.0 + en));    // momentum [p] = m; Gordon, formula (3)
        p2 = p * p;                               // p^2 = p*p
        gamma2 = (1.0 + en) * (1.0 + en);
        invgamma4 = 1.0 / (gamma2 * gamma2);
    }
    
    /* store last entry, since it is needed in computeVerticalOscillations(), because we have to do the same
     * number of integrations steps there.
     */
    lastOrbitVal_m = r_m[N_m];           // needed in computeVerticalOscillations()
    lastMomentumVal_m = pr_m[N_m];       // needed in computeVerticalOscillations()
    
    // remove last entry (since we don't have to store [0,2pi], but [0,2pi[)  --> size = N_m, capacity = N_m+1
    r_m.pop_back();
    pr_m.pop_back();
    

    // returns true if converged, otherwise false
    return error < accuracy;
}

template<typename Value_type, typename Size_type, class Stepper>
Value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::computeTune(const std::array<value_type,2>& y, value_type py2, size_type ncross) {
    // Y = [y1, y2; py1, py2]
    
    // cos(mu)
    value_type cos = 0.5 * (y[0] + py2);
    
    value_type twopi = 2.0 * M_PI;
    value_type mu;
    
    // sign of sin(mu) has to be equal to y2
    bool negative = std::signbit(y[1]);
    
    bool uneven = (ncross % 2);
    
    if (std::fabs(cos) > 1.0) {
        // store the number of crossings
        value_type n = ncross;
        
        if (uneven)
            n = ncross - 1;
        
        // Gordon, formula (36b)
        value_type muPrime = -std::acosh(std::fabs(cos));      // mu'
        mu = n * M_PI + muPrime;
        
    } else {
        value_type muPrime = (uneven) ? std::acos(-cos) : std::acos(cos);    // mu'
        /* It has to be fulfilled: 0<= mu' <= pi
        * But since |cos(mu)| <= 1, we have
        * -1 <= cos(mu) <= 1 --> 0 <= mu <= pi (using above programmed line), such
        * that condition is already fulfilled.
        */
        
        // Gordon, formula (36)
        mu = ncross * M_PI + muPrime;
        
        // if sign(y[1]) > 0 && sign(sin(mu)) < 0
        if (!negative && std::signbit(std::sin(mu))) {
            mu = ncross * M_PI - muPrime;
        } else if (negative && !std::signbit(std::sin(mu))) {
            mu = ncross * M_PI - muPrime + twopi;
        }
    }
    
    // nu = mu/theta, where theta = integration domain
    return mu / twopi;
}

template<typename Value_type, typename Size_type, class Stepper>
void ClosedOrbitFinder<Value_type, Size_type, Stepper>::computeOrbitProperties() {
    /* 
     * The formulas for h, fidx and ds are from the paper:
     * "Tranverse-Longitudinal Coupling by Space Charge in Cyclotrons"
     * written by Dr. Christian Baumgarten (2012, PSI)
     * p. 6
     */
    
    value_type bint, brint, btint; // B, dB/dr, dB/dtheta

    value_type invbcon = 1.0 / physics::bcon(wo_m);
    value_type en = E_m / physics::E0;                                  // en = E/E0 = E/(mc^2) (E0 is kinetic energy)
    value_type p = physics::acon(wo_m) * std::sqrt(en * (2.0 + en));    // momentum [p] = m; Gordon, formula (3)
    value_type p2 = p * p;
    value_type theta = 0.0;                                             // angle for interpolating
    value_type ptheta;

    // resize of container (--> size = N_m, capacity = N_m)
    h_m.resize(N_m);
    fidx_m.resize(N_m);
    ds_m.resize(N_m);

    for (size_type i = 0; i < N_m; ++i) {
        // interpolate magnetic field
        MagneticField::interpolate(&bint,&brint,&btint,theta * 180.0 / M_PI,nr_m,nth_m,r_m[i],rmin_m,dr_m,bmag_m);
        bint *= invbcon;
        brint *= invbcon;
        btint *= invbcon;

        // inverse bending radius
        h_m[i] = bint / p;

        // local field index
        ptheta = std::sqrt(p2 - pr_m[i] * pr_m[i]);
        fidx_m[i] = (brint * ptheta - btint * pr_m[i] / r_m[i]) / p2; //(bint*bint);

        // path length element
        ds_m[i] = std::hypot(r_m[i] * pr_m[i] / ptheta,r_m[i]) * dtheta_m; // C++11 function

        // increase angle
        theta += dtheta_m;
    }

    // compute average radius
    ravg_m = std::accumulate(r_m.begin(),r_m.end(),0.0) / value_type(r_m.size());
}

template<typename Value_type, typename Size_type, class Stepper>
void ClosedOrbitFinder<Value_type, Size_type, Stepper>::computeVerticalOscillations() {
    
    vertOscDone_m = true;

    // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
   
    value_type bint, brint, btint; // B, dB/dr, dB/dtheta

    value_type en = E_m / physics::E0;                                  // en = E/E0 = E/(mc^2) with kinetic energy E0
    value_type p = physics::acon(wo_m) * std::sqrt(en *(en + 2.0));     // Gordon, formula (3)
    value_type p2 = p * p;                                              // p^2 = p*p
    size_type idx = 0;                                                  // index for going through container
    value_type pr2;                                                     // pr^2 = pr*pr
    value_type ptheta, invptheta;                                       // Gordon, formula (5c)
    value_type zold = 0.0;                                              // for counting nzcross

    // store bcon locally
    value_type invbcon = 1.0 / physics::bcon(wo_m);                     // [bcon] = MeV*s/(C*m^2) = 10^6 T = 10^7 kG (kilo Gauss)

    // define the ODEs (using lambda function)
    std::function<void(const state_type&, state_type&, const double)> vertical = [&](const state_type &y, state_type &dydt, const double theta){
        pr2 = y[1] * y[1];
        if (p2 < pr2)
            PhysicalError::message("ClosedOrbitFinder::findOrbit",PhysicalError::negative);
        
        // Gordon, formula (5c)
        ptheta = std::sqrt(p2 - pr2);
        invptheta = 1.0 / ptheta;

        // intepolate values of magnetic field
        MagneticField::interpolate(&bint,&brint,&btint,theta * 180 / M_PI,nr_m,nth_m,y[0],rmin_m,dr_m,bmag_m);
        bint *= invbcon;
        brint *= invbcon;
        btint *= invbcon;
        
        // We have to integrate r and pr again, otherwise we don't have the Runge-Kutta of the B-field
        // Gordon, formula (5a)
        dydt[0] = y[0] * y[1] * invptheta;
        // Gordon, formula (5b)
        dydt[1] = ptheta - y[0] * bint;

        // Gordon, formulas (22a) and (22b)
        for (size_type i = 2; i < 5; i += 2) {
            dydt[i] = y[0] * y[i+1] * invptheta;
            dydt[i+1] = (y[0] * brint - y[1] * invptheta * btint) * y[i];
        }

        // integrate phase
        dydt[6] = y[0] * invptheta * gamma_m - 1;
    };

    // to get next index for r and pr (to iterate over container)
    auto next = [&](state_type& y, const value_type t) {
        // number of times z2 changes sign
        nzcross_m += (idx > 0) * (y[4] * zold < 0);
        zold = y[4];
        ++idx;
    };

    // set initial state container for integration: y = {r, pr, z1, pz1, z2, pz2, phase}
    state_type y = {r_m[0], pr_m[0], 1.0, 0.0, 0.0, 1.0, 0.0};

    // add last element for integration (since we have to return to the initial point (--> size = N_m+1, capacity = N_m+1)
    r_m.push_back(lastOrbitVal_m);
    pr_m.push_back(lastMomentumVal_m);
    
    // integrate: assume no imperfections --> only integrate over a single sector (dtheta_m = 2pi/N_m)
    boost::numeric::odeint::integrate_n_steps(stepper_m,vertical,y,0.0,dtheta_m,N_m,next);
    
    // remove last element again (--> size = N_m, capacity = N_m+1)
    r_m.pop_back();
    pr_m.pop_back();
    
    // write new state
    z_m[0] = y[2];
    pz_m[0] = y[3];
    z_m[1] = y[4];
    pz_m[1] = y[5];
    phase_m = y[6] / (2.0 * M_PI);
}

#endif
