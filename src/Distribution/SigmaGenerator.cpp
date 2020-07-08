//
// Class: SigmaGenerator.h
// The SigmaGenerator class uses the class <b>ClosedOrbitFinder</b> to get the parameters(inverse bending radius, path length
// field index and tunes) to initialize the sigma matrix.
// The main function of this class is <b>match(double, unsigned int)</b>, where it iteratively tries to find a matched
// distribution for given
// emittances, energy and current. The computation stops when the L2-norm is smaller than a user-defined tolerance. \n
// In default mode it prints all space charge maps, cyclotron maps and second moment matrices. The orbit properties, i.e.
// tunes, average radius, orbit radius, inverse bending radius, path length, field index and frequency error, are printed
// as well.
//
// Copyright (c) 2014, 2018, Matthias Frey, Cristopher Cortes, ETH Zürich
// All rights reserved
//
// Implemented as part of the Semester thesis by Matthias Frey
// "Matched Distributions in Cyclotrons"
//
// Some adaptations done as part of the Bachelor thesis by Cristopher Cortes
// "Limitations of a linear transfer map method for finding matched distributions in high intensity cyclotrons"
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "SigmaGenerator.h"

#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include "matrix_vector_operation.h"
#include "ClosedOrbitFinder.h"
#include "MapGenerator.h"

#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <limits>
#include <list>
#include <numeric>
#include <sstream>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

extern Inform *gmsg;

SigmaGenerator::SigmaGenerator(double I,
                               double ex,
                               double ey,
                               double ez,
                               double E,
                               double m,
                               const Cyclotron* cycl,
                               unsigned int N,
                               unsigned int Nsectors,
                               unsigned int truncOrder,
                               bool write)
    : I_m(I)
    , wo_m(cycl->getRfFrequ(0)*1E6/cycl->getCyclHarm()*2.0*Physics::pi)
    , E_m(E)
    , gamma_m(E/m+1.0)
    , gamma2_m(gamma_m*gamma_m)
    , nh_m(cycl->getCyclHarm())
    , beta_m(std::sqrt(1.0-1.0/gamma2_m))
    , m_m(m)
    , niterations_m(0)
    , converged_m(false)
    , Emin_m(cycl->getFMLowE())
    , Emax_m(cycl->getFMHighE())
    , N_m(N)
    , nSectors_m(Nsectors)
    , nStepsPerSector_m(N/cycl->getSymmetry())
    , nSteps_m(N)
    , error_m(std::numeric_limits<double>::max())
    , truncOrder_m(truncOrder)
    , write_m(write)
    , sigmas_m(nStepsPerSector_m)
    , rinit_m(0.0)
    , prinit_m(0.0)
{
    // minimum beta*gamma
    double minGamma = Emin_m / m_m + 1.0;
    double bgam = std::sqrt(minGamma * minGamma - 1.0);

    // set emittances (initialization like that due to old compiler version)
    // [ex] = [ey] = [ez] = pi*mm*mrad --> [emittance] = m rad
    // normalized emittance (--> multiply with beta*gamma)
    emittance_m[0] = ex * Physics::pi * 1.0e-6 * bgam;
    emittance_m[1] = ey * Physics::pi * 1.0e-6 * bgam;
    emittance_m[2] = ez * Physics::pi * 1.0e-6 * bgam;

    // Define the Hamiltonian
    Series::setGlobalTruncOrder(truncOrder_m);

    // infinitesimal elements
    x_m = Series::makeVariable(0);
    px_m = Series::makeVariable(1);
    y_m = Series::makeVariable(2);
    py_m = Series::makeVariable(3);
    z_m = Series::makeVariable(4);
    delta_m = Series::makeVariable(5);

    H_m = [&](double h, double kx, double ky) {
        return 0.5*px_m*px_m + 0.5*kx*x_m*x_m - h*x_m*delta_m +
               0.5*py_m*py_m + 0.5*ky*y_m*y_m +
               0.5*delta_m*delta_m/gamma2_m;
    };

    Hsc_m = [&](double sigx, double sigy, double sigz) {
        // convert m from MeV/c^2 to eV*s^{2}/m^{2}
        double m = m_m * 1.0e6 / (Physics::c * Physics::c);

        // formula (57)
        double lam = Physics::two_pi*Physics::c / (wo_m * nh_m); // wavelength, [lam] = m
        double K3 = 3.0 * /* physics::q0 */ 1.0 * I_m * lam / (20.0 * std::sqrt(5.0) * Physics::pi * Physics::epsilon_0 * m *
                        Physics::c * Physics::c * Physics::c * beta_m * beta_m * gamma_m * gamma2_m);            // [K3] = m

        // formula (30), (31)
        // [sigma(0,0)] = m^{2} --> [sx] = [sy] = [sz] = m
        // multiply with 0.001 to get meter --> [sx] = [sy] = [sz] = m
        double sx = std::sqrt(std::abs(sigx));
        double sy = std::sqrt(std::abs(sigy));
        double sz = std::sqrt(std::abs(sigz));

        double tmp = sx * sy;                                           // [tmp] = m^{2}

        double f = std::sqrt(tmp) / (3.0 * gamma_m * sz);               // [f] = 1
        double kxy = K3 * std::abs(1.0 - f) / ((sx + sy) * sz); // [kxy] = 1/m

        double Kx = kxy / sx;
        double Ky = kxy / sy;
        double Kz = K3 * f / (tmp * sz);

        return -0.5 * Kx * x_m * x_m
               -0.5 * Ky * y_m * y_m
               -0.5 * Kz * z_m * z_m * gamma2_m;
     };
}


bool SigmaGenerator::match(double accuracy,
                           unsigned int maxit,
                           unsigned int maxitOrbit,
                           Cyclotron* cycl,
                           double denergy,
                           double rguess,
                           bool full)
{
    /* compute the equilibrium orbit for energy E_
     * and get the the following properties:
     * - inverse bending radius h
     * - step sizes of path ds
     * - tune nuz
     */

    try {
        if ( !full )
            nSteps_m = nStepsPerSector_m;

        // object for space charge map and cyclotron map
        MapGenerator<double,
                     unsigned int,
                     Series,
                     Map,
                     Hamiltonian,
                     SpaceCharge> mapgen(nSteps_m);

        // compute cyclotron map and space charge map for each angle and store them into a vector
        std::vector<matrix_type> Mcycs(nSteps_m), Mscs(nSteps_m);

        container_type h(nSteps_m), r(nSteps_m), ds(nSteps_m), fidx(nSteps_m);

        ClosedOrbitFinder<double, unsigned int,
            boost::numeric::odeint::runge_kutta4<container_type> > cof(m_m, N_m, cycl, false, nSectors_m);

        if ( !cof.findOrbit(accuracy, maxitOrbit, E_m, denergy, rguess) ) {
            throw OpalException("SigmaGenerator::match()",
                                "Closed orbit finder didn't converge.");
        }

        cof.computeOrbitProperties(E_m);

        // properties of one turn
        std::pair<double,double> tunes = cof.getTunes();
        double ravg = cof.getAverageRadius();                   // average radius

        double angle = cycl->getPHIinit();
        container_type h_turn = cof.getInverseBendingRadius(angle);
        container_type r_turn = cof.getOrbit(angle);
        container_type ds_turn = cof.getPathLength(angle);
        container_type fidx_turn = cof.getFieldIndex(angle);

        container_type peo = cof.getMomentum(angle);

        // write properties to file
        if (write_m)
            writeOrbitOutput_m(tunes, ravg, cof.getFrequencyError(),
                                r_turn, peo, h_turn, fidx_turn, ds_turn);

        // write to terminal
        *gmsg << "* ----------------------------" << endl
                << "* Closed orbit info:" << endl
                << "*" << endl
                << "* average radius: " << ravg << " [m]" << endl
                << "* initial radius: " << r_turn[0] << " [m]" << endl
                << "* initial momentum: " << peo[0] << " [Beta Gamma]" << endl
                << "* frequency error: " << cof.getFrequencyError()*100 <<" [ % ] "<< endl
                << "* horizontal tune: " << tunes.first << endl
                << "* vertical tune: " << tunes.second << endl
                << "* ----------------------------" << endl << endl;

        // copy properties
        std::copy_n(r_turn.begin(), nSteps_m, r.begin());
        std::copy_n(h_turn.begin(), nSteps_m, h.begin());
        std::copy_n(fidx_turn.begin(), nSteps_m, fidx.begin());
        std::copy_n(ds_turn.begin(), nSteps_m, ds.begin());

        rinit_m = r[0];
        prinit_m = peo[0];

        std::string fpath = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "maps"
        });
        if (!boost::filesystem::exists(fpath)) {
            boost::filesystem::create_directory(fpath);
        }

        // initialize sigma matrices (for each angle one) (first guess)
        initialize(tunes.second,ravg);

        // for writing
        std::ofstream writeMturn, writeMcyc, writeMsc;

        if (write_m) {

            std::string energy = float2string(E_m);

            std::string fname = Util::combineFilePath({
                OpalData::getInstance()->getAuxiliaryOutputDirectory(),
                "maps",
                "OneTurnMapsForEnergy" + energy + "MeV.dat"
            });

            writeMturn.open(fname, std::ios::app);

            fname = Util::combineFilePath({
                OpalData::getInstance()->getAuxiliaryOutputDirectory(),
                "maps",
                "SpaceChargeMapPerAngleForEnergy" + energy + "MeV_iteration_0.dat"
            });

            writeMsc.open(fname, std::ios::app);

            fname = Util::combineFilePath({
                OpalData::getInstance()->getAuxiliaryOutputDirectory(),
                "maps",
                "CyclotronMapPerAngleForEnergy" + energy + "MeV.dat"
            });

            writeMcyc.open(fname, std::ios::app);
        }

        // calculate only for a single sector (a nSector_-th) of the whole cyclotron
        for (unsigned int i = 0; i < nSteps_m; ++i) {
            Mcycs[i] = mapgen.generateMap(H_m(h[i],
                                              h[i]*h[i]+fidx[i],
                                              -fidx[i]),
                                          ds[i],truncOrder_m);


            Mscs[i]  = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),
                                                sigmas_m[i](2,2),
                                                sigmas_m[i](4,4)),
                                          ds[i],truncOrder_m);

            writeMatrix(writeMcyc, Mcycs[i]);
            writeMatrix(writeMsc, Mscs[i]);
        }

        if (write_m)
            writeMsc.close();

        // one turn matrix
        mapgen.combine(Mscs,Mcycs);
        matrix_type Mturn = mapgen.getMap();

        writeMatrix(writeMturn, Mturn);

        // (inverse) transformation matrix
        sparse_matrix_type R(6, 6), invR(6, 6);

        // new initial sigma matrix
        matrix_type newSigma(6,6);

        // for exiting loop
        bool stop = false;

        double weight = 0.5;

        while (error_m > accuracy && !stop) {
            // decouple transfer matrix and compute (inverse) tranformation matrix
            decouple(Mturn,R,invR);

            // construct new initial sigma-matrix
            newSigma = updateInitialSigma(Mturn, R, invR);

            // compute new sigma matrices for all angles (except for initial sigma)
            updateSigma(Mscs,Mcycs);

            // compute error
            error_m = L2ErrorNorm(sigmas_m[0],newSigma);

            // write new initial sigma-matrix into vector
            sigmas_m[0] = weight*newSigma + (1.0-weight)*sigmas_m[0];

            // compute new space charge maps
            for (unsigned int i = 0; i < nSteps_m; ++i) {
                Mscs[i] = mapgen.generateMap(Hsc_m(sigmas_m[i](0,0),
                                                    sigmas_m[i](2,2),
                                                    sigmas_m[i](4,4)),
                                                ds[i],truncOrder_m);

                if (write_m) {

                    std::string energy = float2string(E_m);

                    std::string fname = Util::combineFilePath({
                            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
                            "maps",
                            "SpaceChargeMapPerAngleForEnergy" + energy + "MeV_iteration_"
                            + std::to_string(niterations_m + 1)
                            + ".dat"
                    });

                    writeMsc.open(fname, std::ios::out);
                }

                writeMatrix(writeMsc, Mscs[i]);

                if (write_m) {
                    writeMsc.close();
                }
            }

            // construct new one turn transfer matrix M
            mapgen.combine(Mscs,Mcycs);
            Mturn = mapgen.getMap();

            writeMatrix(writeMturn, Mturn);

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
            writeMcyc.close();
        }

    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    if ( converged_m && write_m ) {
        // write tunes
        std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "MatchedDistributions.dat"
        });

        std::ofstream writeSigmaMatched(fname, std::ios::app);

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


void SigmaGenerator::eigsolve_m(const matrix_type& Mturn,
                                sparse_matrix_type& R)
{
    typedef gsl_matrix*                   gsl_matrix_t;
    typedef gsl_vector_complex*           gsl_cvector_t;
    typedef gsl_matrix_complex*           gsl_cmatrix_t;
    typedef gsl_eigen_nonsymmv_workspace* gsl_wspace_t;
    typedef boost::numeric::ublas::vector<complex_t> complex_vector_type;

    gsl_cvector_t evals = gsl_vector_complex_alloc(6);
    gsl_cmatrix_t evecs = gsl_matrix_complex_alloc(6, 6);
    gsl_wspace_t wspace = gsl_eigen_nonsymmv_alloc(6);
    gsl_matrix_t      M = gsl_matrix_alloc(6, 6);

    // go to GSL
    for (unsigned int i = 0; i < 6; ++i){
        for (unsigned int j = 0; j < 6; ++j) {
            gsl_matrix_set(M, i, j, Mturn(i,j));
        }
    }

    /*int status = */gsl_eigen_nonsymmv(M, evals, evecs, wspace);

//     if ( !status )
//         throw OpalException("SigmaGenerator::eigsolve_m()",
//                             "Couldn't perform eigendecomposition!");

    /*status = *///gsl_eigen_nonsymmv_sort(evals, evecs, GSL_EIGEN_SORT_ABS_ASC);

//     if ( !status )
//         throw OpalException("SigmaGenerator::eigsolve_m()",
//                             "Couldn't sort eigenvalues and eigenvectors!");

    // go to UBLAS
    for( unsigned int i = 0; i < 6; i++){
        gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evecs, i);

        for(unsigned int j = 0;j < 6; j++){
            gsl_complex zgsl = gsl_vector_complex_get(&evec_i.vector, j);
            complex_t z(GSL_REAL(zgsl), GSL_IMAG(zgsl));
            R(i,j) = z;
        }
    }

    // Sorting the Eigenvectors
    // This is an arbitrary threshold that has worked for me. (We should fix this)
    double threshold = 10e-12;
    bool isZdirection = false;
    std::vector<complex_vector_type> zVectors{};
    std::vector<complex_vector_type> xyVectors{};

    for(unsigned int i = 0; i < 6; i++){
        complex_t z = R(i,0);
        if(std::abs(z) < threshold) z = 0.;
        if(z == 0.) isZdirection = true;
        complex_vector_type v(6);
        if(isZdirection){
            for(unsigned int j = 0;j < 6; j++){
                complex_t z = R(i,j);
                v(j) = z;
            }
            zVectors.push_back(v);
        }
        else{
            for(unsigned int j = 0; j < 6; j++){
                complex_t z = R(i,j);
                v(j) = z;
            }
            xyVectors.push_back(v);
        }
        isZdirection = false;
    }

    //if z-direction not found, then the system does not behave as expected
    if(zVectors.size() != 2)
        throw OpalException("SigmaGenerator::eigsolve_m()",
                            "Couldn't find the vertical Eigenvectors.");

    // Norms the Eigenvectors
     for(unsigned int i = 0; i < 4; i++){
        double norm{0};
        for(unsigned int j = 0; j < 6; j++) norm += std::norm(xyVectors[i](j));
        for(unsigned int j = 0; j < 6; j++) xyVectors[i](j) /= std::sqrt(norm);
    }
    for(unsigned int i = 0; i < 2; i++){
        double norm{0};
        for(unsigned int j = 0; j < 6; j++) norm += std::norm(zVectors[i](j));
        for(unsigned int j = 0; j < 6; j++) zVectors[i](j) /= std::sqrt(norm);
    }

    for(double i = 0; i < 6; i++){
        R(i,0) = xyVectors[0](i);
        R(i,1) = xyVectors[1](i);
        R(i,2) = zVectors[0](i);
        R(i,3) = zVectors[1](i);
        R(i,4) = xyVectors[2](i);
        R(i,5) = xyVectors[3](i);
    }

    gsl_vector_complex_free(evals);
    gsl_matrix_complex_free(evecs);
    gsl_eigen_nonsymmv_free(wspace);
    gsl_matrix_free(M);
}


bool SigmaGenerator::invertMatrix_m(const sparse_matrix_type& R,
                                    sparse_matrix_type& invR)
{
    typedef boost::numeric::ublas::permutation_matrix<unsigned int> pmatrix_t;

    //creates a working copy of R
    cmatrix_type A(R);

    //permutation matrix for the LU-factorization
    pmatrix_t pm(A.size1());

    //LU-factorization
    int res = lu_factorize(A,pm);

    if( res != 0)
        return false;

    // create identity matrix of invR
    invR.assign(boost::numeric::ublas::identity_matrix<complex_t>(A.size1()));

    // backsubstitute to get the inverse
    boost::numeric::ublas::lu_substitute(A, pm, invR);

    return true;
}


void SigmaGenerator::decouple(const matrix_type& Mturn,
                              sparse_matrix_type& R,
                              sparse_matrix_type& invR)
{
    this->eigsolve_m(Mturn, R);

    if ( !this->invertMatrix_m(R, invR) )
        throw OpalException("SigmaGenerator::decouple()",
                            "Couldn't compute inverse matrix!");
}


double SigmaGenerator::isEigenEllipse(const matrix_type& Mturn,
                                      const matrix_type& sigma)
{
    // compute sigma matrix after one turn
    matrix_type newSigma = matt_boost::gemmm<matrix_type>(Mturn,
                                                          sigma,
                                                          boost::numeric::ublas::trans(Mturn));

    // return L2 error
    return L2ErrorNorm(sigma,newSigma);
}


void SigmaGenerator::initialize(double nuz, double ravg)
{
    /*
     * The initialization is based on the analytical solution of
     * using a spherical symmetric beam in the paper:
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
    double invbg = 1.0 / (beta_m * gamma_m);
    double c2 = Physics::c * Physics::c;

    // convert mass m_m from MeV/c^2 to eV*s^{2}/m^{2}
    double m = m_m * 1.0e6 / c2;        // [m] = eV*s^{2}/m^{2}, [m_m] = MeV/c^2

    // emittance [ex] = [ey] = [ez] = m rad
    double ex = emittance_m[0] * invbg;                        // [ex] = m rad
    double ey = emittance_m[1] * invbg;                        // [ey] = m rad
    double ez = emittance_m[2] * invbg;                        // [ez] = m rad

    // initial guess of emittance, [e] = m rad
    double e = std::cbrt(ex * ey * ez);             // cbrt computes cubic root (C++11) <cmath>

    // cyclotron radius [rcyc] = m
    double rcyc = ravg / beta_m;

    // "average" inverse bending radius
    double h = 1.0 / ravg;            // [h] = 1/m

    // formula (57)
    double lam = Physics::two_pi * Physics::c / (wo_m * nh_m); // wavelength, [lam] = m
    double K3 = 3.0 *  I_m * lam
              / (20.0 * std::sqrt(5.0) * Physics::pi * Physics::epsilon_0 * m
                      * c2 * Physics::c * beta_m * beta_m * gamma2_m * gamma_m);    // [K3] = m

    double alpha =  Physics::mu_0 * I_m / (5.0 * std::sqrt(10.0) * m * Physics::c
                 * gamma_m * nh_m) * std::sqrt(rcyc * rcyc * rcyc / (e * e * e));   // [alpha] = 1/rad --> [alpha] = 1

    double sig0 = std::sqrt(2.0 * rcyc * e) / gamma_m;                              // [sig0] = m*sqrt(rad) --> [sig0] = m

    // formula (56)
    double sig;                                     // [sig] = m
    if (alpha >= 2.5) {
        sig = sig0 * std::cbrt(1.0 + alpha);            // cbrt computes cubic root (C++11) <cmath>
    } else if (alpha >= 0) {
        sig = sig0 * (1 + alpha * (0.25 - 0.03125 * alpha));
    } else {
        throw OpalException("SigmaGenerator::initialize()",
                            "Negative alpha value: " + std::to_string(alpha) + " < 0");
    }

    // K = Kx = Ky = Kz
    double K = K3 * gamma_m / (3.0 * sig * sig * sig);   // formula (46), [K] = 1/m^{2}
    double kx = h * h * gamma2_m;                        // formula (46) (assumption of an isochronous cyclotron), [kx] = 1/m^{2}

    double a = 0.5 * kx - K;    // formula (22) (with K = Kx = Kz), [a] = 1/m^{2}
    double b = K * K;           // formula (22) (with K = Kx = Kz and kx = h^2*gamma^2), [b] = 1/m^{4}


    // b must be positive, otherwise no real-valued frequency
    if (b < 0)
        throw OpalException("SigmaGenerator::initialize()",
                            "Negative value --> No real-valued frequency.");

    double tmp = a * a - b;           // [tmp] = 1/m^{4}
    if (tmp < 0)
        throw OpalException("SigmaGenerator::initialize()",
                            "Square root of negative number.");

    tmp = std::sqrt(tmp);               // [tmp] = 1/m^{2}

    if (a < tmp)
        throw OpalException("Error in SigmaGenerator::initialize()",
                            "Square root of negative number.");

    if (h * h * nuz * nuz <= K)
        throw OpalException("SigmaGenerator::initialize()",
                            "h^{2} * nu_{z}^{2} <= K (i.e. square root of negative number)");

    double Omega = std::sqrt(a + tmp);                // formula (22), [Omega] = 1/m
    double omega = std::sqrt(a - tmp);                // formula (22), [omega] = 1/m

    double A = h / (Omega * Omega + K);           // formula (26), [A] = m
    double B = h / (omega * omega + K);           // formula (26), [B] = m
    double invAB = 1.0 / (B - A);                 // [invAB] = 1/m

    // construct initial sigma-matrix (formula (29, 30, 31)
    matrix_type sigma = boost::numeric::ublas::zero_matrix<double>(6);

    // formula (30), [sigma(0,0)] = m^2 rad
    sigma(0,0) = invAB * (B * ex / Omega + A * ez / omega);

    // [sigma(0,5)] = [sigma(5,0)] = m rad
    sigma(0,5) = sigma(5,0) = invAB * (ex / Omega + ez / omega);

    // [sigma(1,1)] = rad
    sigma(1,1) = invAB * (B * ex * Omega + A * ez * omega);

    // [sigma(1,4)] = [sigma(4,1)] = m rad
    sigma(1,4) = sigma(4,1) = invAB * (ex * Omega+ez * omega) / (K * gamma2_m);

    // formula (31), [sigma(2,2)] = m rad
    sigma(2,2) = ey / (std::sqrt(h * h * nuz * nuz - K));

    sigma(3,3) = (ey * ey) / sigma(2,2);

    // [sigma(4,4)] = m^{2} rad
    sigma(4,4) = invAB * (A * ex * Omega + B * ez * omega) / (K * gamma2_m);

    // formula (30), [sigma(5,5)] = rad
    sigma(5,5) = invAB * (ex / (B * Omega) + ez / (A * omega));

    // fill in initial guess of the sigma matrix (for each angle the same guess)
    sigmas_m.resize(nSteps_m);
    for (typename std::vector<matrix_type>::iterator it = sigmas_m.begin(); it != sigmas_m.end(); ++it) {
        *it = sigma;
    }

    if (write_m) {
        std::string energy = float2string(E_m);

        std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "maps",
            "InitialSigmaPerAngleForEnergy" + energy + "MeV.dat"
        });

        std::ofstream writeInit(fname, std::ios::app);
        writeInit << sigma << std::endl;
        writeInit.close();
    }
}



typename SigmaGenerator::matrix_type
SigmaGenerator::updateInitialSigma(const matrix_type& /*M*/,
                                   sparse_matrix_type& R,
                                   sparse_matrix_type& invR)
{
    /*
     * Function input:
     * - M: one turn transfer matrix
     * - R: transformation matrix (in paper: E)
     * - invR: inverse transformation matrix (in paper: E^{-1}
     */

    /* formula (18):
     * sigma = -E*D*E^{-1}*S
     * with diagonal matrix D (stores eigenvalues of sigma*S (emittances apart from +- i),
     * skew-symmetric matrix (formula (13)), and tranformation matrices E, E^{-1}
     */

    cmatrix_type S = boost::numeric::ublas::zero_matrix<complex_t>(6,6);
    S(0,1) = S(2,3) = S(4,5) = 1;
    S(1,0) = S(3,2) = S(5,4) = -1;

    // Build new D-Matrix
    // Section 2.4 Eq. 18 in M. Frey Semester thesis
    // D = diag(i*emx,-i*emx,i*emy,-i*emy,i*emz, -i*emz)


    cmatrix_type D = boost::numeric::ublas::zero_matrix<complex_t>(6,6);
    double invbg = 1.0 / (beta_m * gamma_m);
    complex_t im(0,1);
    for(unsigned int i = 0; i < 3; ++i){
        D(2*i, 2*i) = emittance_m[i] * invbg * im;
        D(2*i+1, 2*i+1) = -emittance_m[i] * invbg * im;
    }

    // Computing of new Sigma
    // sigma = -R*D*R^{-1}*S
    cmatrix_type csigma(6, 6);
    csigma = boost::numeric::ublas::prod(invR, S);
    csigma = boost::numeric::ublas::prod(D, csigma);
    csigma = boost::numeric::ublas::prod(-R, csigma);

    matrix_type sigma(6,6);
    for (unsigned int i = 0; i < 6; ++i){
        for (unsigned int j = 0; j < 6; ++j){
            sigma(i,j) = csigma(i,j).real();
        }
    }

    for (unsigned int i = 0; i < 6; ++i) {
        if(sigma(i,i) < 0.)
            sigma(i,i) *= -1.0;
    }

    return sigma;
}


void SigmaGenerator::updateSigma(const std::vector<matrix_type>& Mscs,
                                 const std::vector<matrix_type>& Mcycs)
{
    matrix_type M = boost::numeric::ublas::matrix<double>(6,6);

    std::ofstream writeSigma;

    if (write_m) {
        std::string energy = float2string(E_m);

        std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "maps",
            "SigmaPerAngleForEnergy" + energy + "MeV_iteration_"
            + std::to_string(niterations_m + 1)
            + ".dat"
        });

        writeSigma.open(fname,std::ios::app);
    }

    // initial sigma is already computed
    writeMatrix(writeSigma, sigmas_m[0]);

    for (unsigned int i = 1; i < nSteps_m; ++i) {
        // transfer matrix for one angle
        M = boost::numeric::ublas::prod(Mscs[i - 1],Mcycs[i - 1]);
        // transfer the matrix sigma
        sigmas_m[i] = matt_boost::gemmm<matrix_type>(M,sigmas_m[i - 1],
                                                     boost::numeric::ublas::trans(M));

        writeMatrix(writeSigma, sigmas_m[i]);
    }

    if (write_m) {
        writeSigma.close();
    }
}


double SigmaGenerator::L2ErrorNorm(const matrix_type& oldS,
                                   const matrix_type& newS)
{
    // compute difference
    matrix_type diff = newS - oldS;

    // sum squared error up and take square root
    return std::sqrt(std::inner_product(diff.data().begin(),
                                        diff.data().end(),
                                        diff.data().begin(),
                                        0.0));
}


double SigmaGenerator::L1ErrorNorm(const matrix_type& oldS,
                                   const matrix_type& newS)
{
    // compute difference
    matrix_type diff = newS - oldS;

    std::transform(diff.data().begin(), diff.data().end(), diff.data().begin(),
                  [](double& val) {
                      return std::abs(val);
                  });

    return std::accumulate(diff.data().begin(), diff.data().end(), 0.0);
}


double SigmaGenerator::LInfErrorNorm(const matrix_type& oldS,
                                     const matrix_type& newS)
{
    // compute difference
    matrix_type diff = newS - oldS;

    std::transform(diff.data().begin(), diff.data().end(), diff.data().begin(),
                  [](double& val) {
                      return std::abs(val);
                  });

    return *std::max_element(diff.data().begin(), diff.data().end());
}



std::string SigmaGenerator::float2string(double val) {
    std::ostringstream out;
    out << std::setprecision(6) << val;
    return out.str();
}


void SigmaGenerator::writeOrbitOutput_m(
    const std::pair<double,double>& tunes,
    const double& ravg,
    const double& freqError,
    const container_type& r_turn,
    const container_type& peo,
    const container_type& h_turn,
    const container_type&  fidx_turn,
    const container_type&  ds_turn)
{
    std::string fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "Tunes.dat"
    });

    // write tunes
    std::ofstream writeTunes(fname, std::ios::app);

    if(writeTunes.tellp() == 0) // if nothing yet written --> write description
        writeTunes << "energy [MeV]" << std::setw(15)
                   << "nur" << std::setw(25)
                   << "nuz" << std::endl;

    writeTunes << E_m << std::setw(30) << std::setprecision(10)
               << tunes.first << std::setw(25)
               << tunes.second << std::endl;

    // write average radius
    fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "AverageValues.dat"
    });
    std::ofstream writeAvgRadius(fname, std::ios::app);

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
    fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "FrequencyError.dat"
    });
    std::ofstream writePhase(fname, std::ios::app);

    if(writePhase.tellp() == 0) // if nothing yet written --> write description
        writePhase << "energy [MeV]" << std::setw(15)
                   << "freq. error" << std::endl;

    writePhase << E_m << std::setw(30) << std::setprecision(10)
               << freqError << std::endl;

    // write other properties
    std::string energy = float2string(E_m);
    fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "PropertiesForEnergy" + energy + "MeV.dat"
    });
    std::ofstream writeProperties(fname, std::ios::out);
    writeProperties << std::left
                    << std::setw(25) << "orbit radius"
                    << std::setw(25) << "orbit momentum"
                    << std::setw(25) << "inverse bending radius"
                    << std::setw(25) << "field index"
                    << std::setw(25) << "path length" << std::endl;

    for (unsigned int i = 0; i < r_turn.size(); ++i) {
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


void SigmaGenerator::writeMatrix(std::ofstream& os, const matrix_type& m) {
    if (!write_m)
        return;

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 6; ++j)
            os << m(i, j) << " ";
    }
    os << std::endl;
}