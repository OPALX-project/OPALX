//
// Class PartBunch
//   Particle Bunch.
//   A representation of a particle bunch as a vector of particles.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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

#include "Algorithms/PartBunch.h"

#include <cfloat>
#include <memory>
#include <utility>
#include "Ippl.h"

#ifndef ADATEST
// #include "Distribution/Distribution.h"
// #include "Particle/ParticleBalancer.h"
// #include "Structure/FieldSolver.h"
#include "Utilities/GeneralClassicException.h"
#endif

#ifdef DBG_SCALARFIELD
#include "Structure/FieldWriter.h"
#endif

extern Inform* gmsg;

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::initialize(FieldLayout_t<Dim>& fl, Mesh_t<Dim>& mesh) {
    updateLayout(fl, mesh);
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::runTests() {
    Vector_t<T, Dim> ll(-0.005);
    Vector_t<T, Dim> ur(0.005);

    setBCAllPeriodic();

    ippl::NDIndex<Dim> domain = getFieldLayout().getDomain();
    for (unsigned int i = 0; i < Dim; i++)
        nr_m[i] = domain[i].length();

    for (int i = 0; i < 3; i++)
        hr_m[i] = (ur[i] - ll[i]) / nr_m[i];

    getFieldLayout().get_Mesh().set_meshSpacing(&(hr_m[0]));
    getFieldLayout().get_Mesh().set_origin(ll);

    rho_m.initialize(getFieldLayout().getMesh(), getFieldLayout(), 1);
    eg_m.initialize(getFieldLayout().getMesh(), getFieldLayout(), 1);

    this->fs_m->solver_m->test(this);
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::do_binaryRepart() {
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::computeSelfFields(int binNumber) {
    /// Set initial charge density to zero. Create image charge
    /// potential field.
    rho_m               = binNumber;  // to make the compiler happy
    rho_m               = 0.0;
    auto imagePotential = rho_m;

    /// Set initial E field to zero.
    eg_m = Vector_t<T, Dim>(0.0);

    if (this->fs_m->hasValidSolver()) {
        /// Mesh the whole domain
        // resizeMesh();

        /// Scatter charge onto space charge grid.
        this->Q *= this->dt;

        /*
        if (!interpolationCacheSet_m) {
            if (interpolationCache_m.size() < getLocalNum()) {
                interpolationCache_m.create(getLocalNum() - interpolationCache_m.size());
            } else {
                interpolationCache_m.destroy(interpolationCache_m.size() - getLocalNum(),
                                             getLocalNum(), true);
            }
            interpolationCacheSet_m = true;
            this->Q.scatter(this->rho_m, this->R, IntrplCIC_t(), interpolationCache_m);
        } else {
            this->Q.scatter(this->rho_m, IntrplCIC_t(), interpolationCache_m);
        }
        */

        this->Q /= this->dt;
        // this->rho_m /= getdT();

        /// Calculate mesh-scale factor and get gamma for this specific slice (bin).
        double scaleFactor = 1;
        // double scaleFactor = Physics::c * getdT();
        double gammaz = 1.0;  // getBinGamma(binNumber);

        /// Scale charge density to get charge density in real units. Account for
        /// Lorentz transformation in longitudinal direction.
        double tmp2 = 1 / hr_m[0] * 1 / hr_m[1] * 1 / hr_m[2]
                      / (scaleFactor * scaleFactor * scaleFactor) / gammaz;
        rho_m *= tmp2;

        /// Scale mesh spacing to real units (meters). Lorentz transform the
        /// longitudinal direction.
        Vector_t<T, Dim> hr_scaled = hr_m * Vector_t<T, Dim>(scaleFactor);
        hr_scaled[2] *= gammaz;

        /// Find potential from charge in this bin (no image yet) using Poisson's
        /// equation (without coefficient: -1/(eps)).
        imagePotential = rho_m;

        this->fs_m->solver_m->computePotential(rho_m, hr_scaled);

        /// Scale mesh back (to same units as particle locations.)
        rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// The scalar potential is given back in rho_m
        /// and must be converted to the right units.
        // rho_m *= getCouplingConstant();

        /// IPPL Grad numerical computes gradient to find the
        /// electric field (in bin rest frame).
        eg_m = -Grad(rho_m, eg_m);

        /// Scale field. Combine Lorentz transform with conversion to proper field
        /// units.
        eg_m *= Vector_t<T, Dim>(
            gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz));

        // If desired write E-field and potential to terminal
#ifdef FIELDSTDOUT
        // Immediate debug output:
        // Output potential and e-field along the x-, y-, and z-axes
        int mx  = (int)nr_m[0];
        int mx2 = (int)nr_m[0] / 2;
        int my  = (int)nr_m[1];
        int my2 = (int)nr_m[1] / 2;
        int mz  = (int)nr_m[2];
        int mz2 = (int)nr_m[2] / 2;

        for (int i = 0; i < mx; i++)
            *gmsg << "Bin " << binNumber << ", Self Field along x axis E = " << eg_m[i][my2][mz2]
                  << ", Pot = " << rho_m[i][my2][mz2] << endl;

        for (int i = 0; i < my; i++)
            *gmsg << "Bin " << binNumber << ", Self Field along y axis E = " << eg_m[mx2][i][mz2]
                  << ", Pot = " << rho_m[mx2][i][mz2] << endl;

        for (int i = 0; i < mz; i++)
            *gmsg << "Bin " << binNumber << ", Self Field along z axis E = " << eg_m[mx2][my2][i]
                  << ", Pot = " << rho_m[mx2][my2][i] << endl;
#endif

        /// Interpolate electric field at particle positions.  We reuse the
        /// cached information about where the particles are relative to the
        /// field, since the particles have not moved since this the most recent
        /// scatter operation.
        eg_m = 0.0;
        gather(this->E, eg_m, this->R);

        // Eftmp.gather(eg_m, IntrplCIC_t(), interpolationCache_m);
        // Eftmp.gather(eg_m, this->R, IntrplCIC_t());

        /** Magnetic field in x and y direction induced by the electric field.
         *
         *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} =
         * -\frac{beta}{c}E_y \f] \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma
         * \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f] \f[ B_z = B_z^{'} = 0 \f]
         *
         */
        double betaC = std::sqrt(gammaz * gammaz - 1.0) / gammaz;  // / Physics::c;

        Bf(0) = Bf(0) - betaC * Eftmp(1);
        Bf(1) = Bf(1) + betaC * Eftmp(0);

        Ef += Eftmp;

        /// Now compute field due to image charge. This is done separately as the image charge
        /// is moving to -infinity (instead of +infinity) so the Lorentz transform is different.

        /// Find z shift for shifted Green's function.
        ippl::NDIndex<3> domain = getFieldLayout().getDomain();
        Vector_t<T, Dim> origin = rho_m.get_mesh().get_origin();
        double hz               = rho_m.get_mesh().get_meshSpacing(2);
        double zshift = -(2 * origin(2) + (domain[2].first() + domain[2].last() + 1) * hz) * gammaz
                        * scaleFactor;

        /// Find potential from image charge in this bin using Poisson's
        /// equation (without coefficient: -1/(eps)).
        this->fs_m->solver_m->computePotential(imagePotential, hr_scaled, zshift);

        /// Scale mesh back (to same units as particle locations.)
        imagePotential *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];

        /// The scalar potential is given back in rho_m
        /// and must be converted to the right units.
        imagePotential *= 1.0;  // getCouplingConstant();

#ifdef DBG_SCALARFIELD
        const int dumpFreq = 100;
        VField_t tmp_eg    = eg_m;

        if ((localTrackStep_m + 1) % dumpFreq == 0) {
            FieldWriter fwriter;
            fwriter.dumpField(rho_m, "phi", "V", localTrackStep_m / dumpFreq, &imagePotential);
        }
#endif

        /// IPPL Grad numerical computes gradient to find the
        /// electric field (in rest frame of this bin's image
        /// charge).
        eg_m = -Grad(imagePotential, eg_m);

        /// Scale field. Combine Lorentz transform with conversion to proper field
        /// units.
        eg_m *= Vector_t<T, Dim>(
            gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz));

        // If desired write E-field and potential to terminal
#ifdef FIELDSTDOUT
        // Immediate debug output:
        // Output potential and e-field along the x-, y-, and z-axes
        // int mx = (int)nr_m[0];
        // int mx2 = (int)nr_m[0] / 2;
        // int my = (int)nr_m[1];
        // int my2 = (int)nr_m[1] / 2;
        // int mz = (int)nr_m[2];
        // int mz2 = (int)nr_m[2] / 2;

        for (int i = 0; i < mx; i++)
            *gmsg << "Bin " << binNumber << ", Image Field along x axis E = " << eg_m[i][my2][mz2]
                  << ", Pot = " << rho_m[i][my2][mz2] << endl;

        for (int i = 0; i < my; i++)
            *gmsg << "Bin " << binNumber << ", Image Field along y axis E = " << eg_m[mx2][i][mz2]
                  << ", Pot = " << rho_m[mx2][i][mz2] << endl;

        for (int i = 0; i < mz; i++)
            *gmsg << "Bin " << binNumber << ", Image Field along z axis E = " << eg_m[mx2][my2][i]
                  << ", Pot = " << rho_m[mx2][my2][i] << endl;
#endif

#ifdef DBG_SCALARFIELD
        tmp_eg += eg_m;
        if ((localTrackStep_m + 1) % dumpFreq == 0) {
            FieldWriter fwriter;
            fwriter.dumpField(tmp_eg, "e", "V/m", localTrackStep_m / dumpFreq);
        }
#endif

        /// Interpolate electric field at particle positions.  We reuse the
        /// cached information about where the particles are relative to the
        /// field, since the particles have not moved since this the most recent
        /// scatter operation.

        gather(this->E, eg_m, this->R);

        //        Eftmp.gather(eg_m, IntrplCIC_t(), interpolationCache_m);
        // Eftmp.gather(eg_m, this->R, IntrplCIC_t());

        /** Magnetic field in x and y direction induced by the image charge electric field. Note
         * that beta will have the opposite sign from the bunch charge field, as the image charge is
         * moving in the opposite direction.
         *
         *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} =
         * -\frac{beta}{c}E_y \f] \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma
         * \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f] \f[ B_z = B_z^{'} = 0 \f]
         *
         */
        Bf(0) = Bf(0) + betaC * Eftmp(1);
        Bf(1) = Bf(1) - betaC * Eftmp(0);

        Ef += Eftmp;
    }
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::resizeMesh() {
    // when OPAL had AMR we need to do this
    return;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::computeSelfFields() {
    rho_m = 0.0;
    eg_m  = Vector_t<T, Dim>(0.0);

    if (this->fs_m->hasValidSolver()) {
        // mesh the whole domain
        resizeMesh();

        // scatter charges onto grid
        this->Q *= this->dt;
        scatter(Q, rho_m, this->R);
        this->Q /= this->dt;
        //        this->rho_m /= getdT();

        // calculating mesh-scale factor
        double gammaz = sum(this->P)[2] / this->getTotalNum();
        gammaz *= gammaz;
        gammaz             = std::sqrt(gammaz + 1.0);
        double scaleFactor = 1;
        // double scaleFactor = Physics::c * getdT();
        // and get meshspacings in real units [m]
        Vector_t<T, Dim> hr_scaled = hr_m * Vector_t<T, Dim>(scaleFactor);
        hr_scaled[2] *= gammaz;

        // double tmp2 = 1/hr_m[0] * 1/hr_m[1] * 1/hr_m[2] / (scaleFactor*scaleFactor*scaleFactor) /
        // gammaz;
        double tmp2 = 1 / hr_scaled[0] * 1 / hr_scaled[1] * 1 / hr_scaled[2];
        // divide charge by a 'grid-cube' volume to get [C/m^3]
        rho_m *= tmp2;

#ifdef DBG_SCALARFIELD
        FieldWriter fwriter;
        fwriter.dumpField(rho_m, "rho", "C/m^3", localTrackStep_m);
#endif

        // charge density is in rho_m
        this->fs_m->solver_m->computePotential(rho_m, hr_scaled);

        // do the multiplication of the grid-cube volume coming
        // from the discretization of the convolution integral.
        // this is only necessary for the FFT solver
        // FIXME: later move this scaling into FFTPoissonSolver
        // if (this->fs_m->getFieldSolverType() == FieldSolverType::FFT
        //    || this->fs_m->getFieldSolverType() == FieldSolverType::FFTBOX) {
        rho_m *= hr_scaled[0] * hr_scaled[1] * hr_scaled[2];
        //}

        // the scalar potential is given back in rho_m in units
        // [C/m] = [F*V/m] and must be divided by
        // 4*pi*\epsilon_0 [F/m] resulting in [V]
        // rho_m *= getCouplingConstant();

        // write out rho
#ifdef DBG_SCALARFIELD
        fwriter.dumpField(rho_m, "phi", "V", localTrackStep_m);
#endif

        // IPPL Grad divides by hr_m [m] resulting in
        // [V/m] for the electric field
        eg_m = -Grad(rho_m, eg_m);

        eg_m *= Vector_t<T, Dim>(
            gammaz / (scaleFactor), gammaz / (scaleFactor), 1.0 / (scaleFactor * gammaz));

        // write out e field
#ifdef FIELDSTDOUT
        // Immediate debug output:
        // Output potential and e-field along the x-, y-, and z-axes
        int mx  = (int)nr_m[0];
        int mx2 = (int)nr_m[0] / 2;
        int my  = (int)nr_m[1];
        int my2 = (int)nr_m[1] / 2;
        int mz  = (int)nr_m[2];
        int mz2 = (int)nr_m[2] / 2;

        for (int i = 0; i < mx; i++)
            *gmsg << "Field along x axis Ex = " << eg_m[i][my2][mz2]
                  << " Pot = " << rho_m[i][my2][mz2] << endl;

        for (int i = 0; i < my; i++)
            *gmsg << "Field along y axis Ey = " << eg_m[mx2][i][mz2]
                  << " Pot = " << rho_m[mx2][i][mz2] << endl;

        for (int i = 0; i < mz; i++)
            *gmsg << "Field along z axis Ez = " << eg_m[mx2][my2][i]
                  << " Pot = " << rho_m[mx2][my2][i] << endl;
#endif

#ifdef DBG_SCALARFIELD
        fwriter.dumpField(eg_m, "e", "V/m", localTrackStep_m);
#endif

        // interpolate electric field at particle positions.  We reuse the
        // cached information about where the particles are relative to the
        // field, since the particles have not moved since this the most recent
        // scatter operation.

        gather(this->E, eg_m, this->R);

        /** Magnetic field in x and y direction induced by the eletric field
         *
         *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} =
         * -\frac{beta}{c}E_y \f] \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma
         * \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f] \f[ B_z = B_z^{'} = 0 \f]
         *
         */
        double betaC = std::sqrt(gammaz * gammaz - 1.0) / gammaz;  // / Physics::c;

        Bf(0) = Bf(0) - betaC * Ef(1);
        Bf(1) = Bf(1) + betaC * Ef(0);
    }
}

template <class PLayout, typename T, unsigned Dim>
FieldLayout_t<Dim>& PartBunch<PLayout, T, Dim>::getFieldLayout() {
    // PLayout_t* layout = static_cast<Layout_t*>(&getLayout());
    //   return dynamic_cast<FieldLayout_t<Dim>&>(layout->getLayout().getFieldLayout());
    return this->getFieldLayout();
}

/*void PartBunch<PLayout,T,dim>::setBCAllPeriodic() {
    for (int i = 0; i < 2 * 3; ++i) {
        if (Ippl::getNodes() > 1) {
            bc_m[i] = new ParallelInterpolationFace<double, Dim, Mesh_t<Dim>, Center_t>(i);
            // std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new ParallelPeriodicFace<Vector_t<T,Dim>, Dim, Mesh_t<Dim>,
Center_t>(i); } else { bc_m[i] = new InterpolationFace<double, Dim, Mesh_t<Dim>, Center_t>(i);
            // std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new PeriodicFace<Vector_t<T,Dim>, Dim, Mesh_t<Dim>, Center_t>(i);
        }
        getBConds()[i] = ParticlePeriodicBCond;
    }
    dcBeam_m = true;
    *ippl::Info << level3 << "BC set P3M, all periodic" << endl;
}
*/
template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setBCAllOpen() {
}
template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setBCForDCBeam() {
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::updateDomainLength(Vector_t<int, 3>& grid) {
    ippl::NDIndex<3> domain = getFieldLayout().getDomain();
    for (unsigned int i = 0; i < Dim; i++)
        grid[i] = domain[i].length();
}
template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::updateFields(
    const Vector_t<T, Dim>& /*hr*/, const Vector_t<T, Dim>& origin) {
    this->getMesh().set_meshSpacing(&(hr_m[0]));
    this->getMesh().set_origin(origin);
    rho_m.initialize(this->getMesh(), getFieldLayout(), 1);
    eg_m.initialize(this->getMesh(), getFieldLayout(), 1);
}

template <class PLayout, typename T, unsigned Dim>
inline void PartBunch<PLayout, T, Dim>::resetInterpolationCache(bool clearCache) {
    //    if (clearCache) {
    //      interpolationCache_m.destroy(interpolationCache_m.size(), 0, true);
    //  }
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::switchToUnitlessPositions(bool use_dt_per_particle) {
    bool hasToReset = false;
    if (!this->R.isDirty())
        hasToReset = true;

    for (size_t i = 0; i < this->getLocalNum(); i++) {
        double dt = getdT();
        if (use_dt_per_particle)
            dt = this->dt(i);

        // this->R(i) /= Vector_t(Physics::c * dt);
        this->R(i) /= Vector_t<T, Dim>(dt);
    }

    unit_state_m = unitless;

    if (hasToReset)
        this->R.resetDirtyFlag();
}

// ADA

// FIXME: unify methods, use convention that all particles have own dt
template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::switchOffUnitlessPositions(bool use_dt_per_particle) {
    bool hasToReset = false;
    if (!this->R.isDirty())
        hasToReset = true;

    for (size_t i = 0; i < this->getLocalNum(); i++) {
        double dt = getdT();
        if (use_dt_per_particle)
            dt = this->dt(i);

        // R(i) *= Vector_t<T,Dim>(Physics::c * dt);
        this->R(i) *= Vector_t<T, Dim>(dt);
    }

    unit_state_m = units;

    if (hasToReset)
        this->R.resetDirtyFlag();
}

template <class PLayout, typename T, unsigned Dim>
bool PartBunch<PLayout, T, Dim>::isGridFixed() const {
    return fixed_grid;
}

template <class PLayout, typename T, unsigned Dim>
bool PartBunch<PLayout, T, Dim>::hasBinning() const {
    return (pbin_m != nullptr);
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::updateNumTotal() {
    size_type total_particles = 0;
    size_type local_particles = this->getLocalNum();

    MPI_Reduce(
        &local_particles, &total_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0,
        ippl::Comm->getCommunicator());
    // need to be stored
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::rebin() {
    this->Bin = 0;
    // ADA pbin_m->resetBins();
    // delete pbin_m; we did not allocate it!
    pbin_m = nullptr;
}

template <class PLayout, typename T, unsigned Dim>
int PartBunch<PLayout, T, Dim>::getLastemittedBin() {
    if (pbin_m != nullptr)
        return 0;  // ADA pbin_m->getLastemittedBin();
    else
        return 0;
}

template <class PLayout, typename T, unsigned Dim>
int PartBunch<PLayout, T, Dim>::getNumberOfEnergyBins() {
    return 0;
}

template <class PLayout, typename T, unsigned Dim>
size_t PartBunch<PLayout, T, Dim>::emitParticles(double eZ) {
    return 0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getEmissionDeltaT() {
    return 0.;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getBinGamma(int bin) {
    return 1.0;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::Rebin() {
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::calcGammas() {
}

template <class PLayout, typename T, unsigned Dim>
bool PartBunch<PLayout, T, Dim>::weHaveEnergyBins() {
    return false;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setLocalBinCount(size_t num, int bin) {
    if (pbin_m != nullptr) {
        // ADA pbin_m->setPartNum(bin, num);
    }
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setBinCharge(int bin, double q) {
    this->Q = where(eq(this->Bin, bin), q, 0.0);
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setBinCharge(int bin) {
    this->Q = where(eq(this->Bin, bin), this->qi_m, 0.0);
}

template <class PLayout, typename T, unsigned Dim>
size_t PartBunch<PLayout, T, Dim>::calcNumPartsOutside(Vector_t<T, Dim> x) {
    std::size_t localnum         = 0;
    const Vector_t<T, Dim> meanR = get_rmean();

    for (unsigned long k = 0; k < this->getLocalNum(); ++k)
        if (std::abs(this->R(k)(0) - meanR(0)) > x(0) || std::abs(this->R(k)(1) - meanR(1)) > x(1)
            || std::abs(this->R(k)(2) - meanR(2)) > x(2)) {
            ++localnum;
        }

    // ADA gather(&localnum, &globalPartPerNode_m[0], 1);

    size_t npOutside = std::accumulate(
        globalPartPerNode_m.get(), globalPartPerNode_m.get() + ippl::Comm->size(), 0,
        std::plus<size_t>());

    return npOutside;
}

/**
 * \method calcLineDensity()
 * \brief calculates the 1d line density (not normalized) and append it to a file.
 * \see ParallelTTracker
 * \warning none yet
 *
 * DETAILED TODO
 *
 */
template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::calcLineDensity(
    unsigned int nBins, std::vector<double>& lineDensity, std::pair<double, double>& meshInfo) {
    Vector_t<T, Dim> rmin, rmax;
    get_bounds(rmin, rmax);

    if (nBins < 2) {
        Vector<int, 3> /*NDIndex<3>*/ grid;
        this->updateDomainLength(grid);
        nBins = grid[2];
    }

    double length = rmax(2) - rmin(2);
    double zmin = rmin(2) - dh_m * length, zmax = rmax(2) + dh_m * length;
    double hz       = (zmax - zmin) / (nBins - 2);
    double perMeter = 1.0 / hz;  //(zmax - zmin);
    zmin -= hz;

    lineDensity.resize(nBins, 0.0);
    std::fill(lineDensity.begin(), lineDensity.end(), 0.0);

    const unsigned int lN = this->getLocalNum();
    for (unsigned int i = 0; i < lN; ++i) {
        const double z   = this->R(i)(2) - 0.5 * hz;
        unsigned int idx = (z - zmin) / hz;
        double tau       = (z - zmin) / hz - idx;

        // ADA lineDensity[idx] += Q(i) * (1.0 - tau) * perMeter;
        // lineDensity[idx + 1] += Q(i) * tau * perMeter;
    }

    // ADA reduce(&(lineDensity[0]), &(lineDensity[0]) + nBins, &(lineDensity[0]), OpAddAssign());

    meshInfo.first  = zmin;
    meshInfo.second = hz;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::boundp() {
    /*
      Assume rmin_m < 0.0
     */

    // if (!R.isDirty() && stateOfLastBoundP_ == unit_state_) return;
    if (!(this->R.isDirty() || this->ID.isDirty()) && stateOfLastBoundP_m == unit_state_m)
        return;  //-DW

    stateOfLastBoundP_m = unit_state_m;

    if (!isGridFixed()) {
        const int dimIdx = (dcBeam_m ? 2 : 3);

        /**
            In case of dcBeam_m && hr_m < 0
            this is the first call to boundp and we
            have to set hr completely i.e. x,y and z.
         */

        this->updateDomainLength(nr_m);
        get_bounds(rmin_m, rmax_m);
        Vector_t<T, Dim> len = rmax_m - rmin_m;

        double volume = 1.0;
        for (int i = 0; i < dimIdx; i++) {
            double length = std::abs(rmax_m[i] - rmin_m[i]);
            if (length < 1e-10) {
                rmax_m[i] += 1e-10;
                rmin_m[i] -= 1e-10;
            } else {
                rmax_m[i] += dh_m * length;
                rmin_m[i] -= dh_m * length;
            }
            hr_m[i] = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
        }
        if (dcBeam_m) {
            rmax_m[2] = rmin_m[2] + periodLength_m;
            hr_m[2]   = periodLength_m / (nr_m[2] - 1);
        }
        for (int i = 0; i < dimIdx; ++i) {
            volume *= std::abs(rmax_m[i] - rmin_m[i]);
        }

        // ada if (volume < 1e-21 && this->this->getTotalNum() > 1 && std::abs(sum(Q)) > 0.0) {
        //    WARNMSG(level1 << "!!! Extremely high particle density detected !!!" << endl);
        // }
        // *ippl::Info << "It is a full boundp hz= " << hr_m << " rmax= " << rmax_m << " rmin= " <<
        // rmin_m
        // << endl;

        // ADA if (hr_m[0] * hr_m[1] * hr_m[2] <= 0) {
        // throw GeneralClassicException("boundp() ", "h<0, can not build a mesh");
        // }

        Vector_t<T, Dim> origin =
            rmin_m - Vector_t<T, Dim>(hr_m[0] / 2.0, hr_m[1] / 2.0, hr_m[2] / 2.0);
        this->updateFields(hr_m, origin);
    }
    this->update();
    this->R.resetDirtyFlag();
}

template <class PLayout, typename T, unsigned Dim>
size_t PartBunch<PLayout, T, Dim>::boundp_destroyT() {
    this->updateDomainLength(nr_m);

    std::vector<size_t> tmpbinemitted;

    boundp();

    size_t ne             = 0;
    const size_t localNum = this->getLocalNum();

    double rzmean             = 0.0;  // momentsComputer_m.getMeanPosition()(2);
    double rzrms              = 0.0;  // momentsComputer_m.getStandardDeviationPosition()(2);
    const bool haveEnergyBins = weHaveEnergyBins();
    if (haveEnergyBins) {
        tmpbinemitted.resize(getNumberOfEnergyBins(), 0.0);
    }
    for (unsigned int i = 0; i < localNum; i++) {
        if (this->Bin(i)
            < 0) {  // ADA || (Options::remotePartDel > 0
                    // && std::abs(R(i)(2) - rzmean) < Options::remotePartDel * rzrms)) {
            ne++;
            this->destroy(1, i);
        } else if (haveEnergyBins) {
            tmpbinemitted[this->Bin(i)]++;
        }
    }

    boundp();

    calcBeamParameters();
    gatherLoadBalanceStatistics();

    // ADA reduce(ne, ne, OpAddAssign());
    return ne;
}

template <class PLayout, typename T, unsigned Dim>
size_t PartBunch<PLayout, T, Dim>::destroyT() {
    std::unique_ptr<size_t[]> tmpbinemitted;

    const size_t localNum = this->getLocalNum();
    const size_t totalNum = this->getTotalNum();
    size_t ne             = 0;

    if (weHaveEnergyBins()) {
        tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[getNumberOfEnergyBins()]);
        for (int i = 0; i < getNumberOfEnergyBins(); i++) {
            tmpbinemitted[i] = 0;
        }
        for (unsigned int i = 0; i < localNum; i++) {
            if (this->Bin(i) < 0) {
                // ADA destroy(1, i);
                ++ne;
            } else
                tmpbinemitted[this->Bin(i)]++;
        }
    } else {
        Inform dmsg("destroy: ", INFORM_ALL_NODES);
        for (size_t i = 0; i < localNum; i++) {
            if ((this->Bin(i) < 0)) {
                ne++;
                // ADA destroy(1, i);
            }
        }
    }
    /* ADA
    if (ne > 0) {
        performDestroy(true);
    }
    */

    calcBeamParameters();
    gatherLoadBalanceStatistics();

    size_t newTotalNum = this->getLocalNum();
    // ADA reduce(newTotalNum, newTotalNum, OpAddAssign());

    // ADA setTotalNum(newTotalNum);

    return totalNum - newTotalNum;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getPx(int /*i*/) {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getPy(int) {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getPz(int) {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getPx0(int) {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getPy0(int) {
    return 0;
}

// ff
template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getX(int i) {
    return this->R(i)(0);
}

// ff
template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getY(int i) {
    return this->R(i)(1);
}

// ff
template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getZ(int i) {
    return this->R(i)(2);
}

// ff
template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getX0(int /*i*/) {
    return 0.0;
}

// ff
template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getY0(int /*i*/) {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setZ(int /*i*/, double /*zcoo*/){};

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::get_bounds(Vector_t<T, Dim>& rmin, Vector_t<T, Dim>& rmax) const {
    this->getLocalBounds(rmin, rmax);

    double min[2 * Dim];

    for (unsigned int i = 0; i < Dim; ++i) {
        min[2 * i]     = rmin[i];
        min[2 * i + 1] = -rmax[i];
    }

    // ADA allreduce(min, 2 * Dim, std::less<double>());

    for (unsigned int i = 0; i < Dim; ++i) {
        rmin[i] = min[2 * i];
        rmax[i] = -min[2 * i + 1];
    }
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::getLocalBounds(
    Vector_t<T, Dim>& rmin, Vector_t<T, Dim>& rmax) const {
    const size_t localNum = this->getLocalNum();
    if (localNum == 0) {
        double maxValue = 1e8;
        rmin            = Vector_t<T, Dim>({maxValue, maxValue, maxValue});
        rmax            = Vector_t<T, Dim>({-maxValue, -maxValue, -maxValue});
        return;
    }

    rmin = this->R(0);
    rmax = this->R(0);
    for (size_t i = 1; i < localNum; ++i) {
        for (unsigned short d = 0; d < 3u; ++d) {
            if (rmin(d) > this->R(i)(d))
                rmin(d) = this->R(i)(d);
            else if (rmax(d) < this->R(i)(d))
                rmax(d) = this->R(i)(d);
        }
    }
}

template <class PLayout, typename T, unsigned Dim>
std::pair<Vector_t<T, Dim>, double> PartBunch<PLayout, T, Dim>::getBoundingSphere() {
    Vector_t<T, Dim> rmin, rmax;
    get_bounds(rmin, rmax);

    std::pair<Vector_t<T, Dim>, double> sphere;
    sphere.first  = 0.5 * (rmin + rmax);
    sphere.second = std::sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}

template <class PLayout, typename T, unsigned Dim>
std::pair<Vector_t<T, Dim>, double> PartBunch<PLayout, T, Dim>::getLocalBoundingSphere() {
    Vector_t<T, Dim> rmin, rmax;
    getLocalBounds(rmin, rmax);

    std::pair<Vector_t<T, Dim>, double> sphere;
    sphere.first  = 0.5 * (rmin + rmax);
    sphere.second = std::sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setdT(double dt) {
    dt_m = dt;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getdT() const {
    return dt_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setT(double t) {
    t_m = t;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::incrementT() {
    t_m += dt_m;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getT() const {
    return t_m;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::get_sPos() const {
    return spos_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::set_sPos(double s) {
    spos_m = s;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::get_gamma() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::get_meanKineticEnergy() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_origin() const {
    return rmin_m;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_maxExtent() const {
    return rmax_m;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_centroid() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_rrms() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_rprms() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_rmean() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_prms() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_pmean() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_emit() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_norm_emit() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_halo() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_68Percentile() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_95Percentile() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_99Percentile() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_99_99Percentile() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_normalizedEps_68Percentile() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_normalizedEps_95Percentile() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_normalizedEps_99Percentile() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_normalizedEps_99_99Percentile() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::get_Dx() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::get_Dy() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::get_DDx() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::get_DDy() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::get_hr() const {
    return hr_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::set_meshEnlargement(double dh) {
    dh_m = dh;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::gatherLoadBalanceStatistics() {
    for (int i = 0; i < ippl::Comm->size(); i++)
        globalPartPerNode_m[i] = 0;

    std::size_t localnum = this->getLocalNum();
    // ADA gather(&localnum, &globalPartPerNode_m[0], 1);
}

template <class PLayout, typename T, unsigned Dim>
size_t PartBunch<PLayout, T, Dim>::getLoadBalance(int p) const {
    return globalPartPerNode_m[p];
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::get_PBounds(Vector_t<T, Dim>& min, Vector_t<T, Dim>& max) const {
    bounds(this->P, min, max);
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::calcBeamParameters() {
    get_bounds(rmin_m, rmax_m);
    // momentsComputer_m.compute(*this);
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getCouplingConstant() const {
    return couplingConstant_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setCouplingConstant(double c) {
    couplingConstant_m = c;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setCharge(double q) {
    qi_m = q;
    if (this->getTotalNum() != 0)
        Q = qi_m;
    // ADA else
    //    WARNMSG("Could not set total charge in PartBunch::setCharge based on this->getTotalNum");
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setChargeZeroPart(double q) {
    qi_m = q;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setMass(double mass) {
    massPerParticle_m = mass;
    if (this->getTotalNum() != 0)
        M = mass;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setMassZeroPart(double mass) {
    massPerParticle_m = mass;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getCharge() const {
    return 0.0;  // ADA sum(this->Q);
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getChargePerParticle() const {
    return qi_m;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getMassPerParticle() const {
    return massPerParticle_m;
}
/*

  template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setSolver(FieldSolver* fs) {
    fs_m = fs;
    fs_m->initSolver(this);


if (!OpalData::getInstance()->hasBunchAllocated()) {
    this->initialize(fs_m->getFieldLayout());
    //         this->setMesh(fs_m->getMesh());
    //         this->setFieldLayout(fs_m->getFieldLayout());
}
}
*/
template <class PLayout, typename T, unsigned Dim>
bool PartBunch<PLayout, T, Dim>::hasFieldSolver() {
    if (this->fs_m)
        return this->fs_m->hasValidSolver();
    else
        return false;
}

/// \brief Return the fieldsolver type if we have a fieldsolver
/*
template <class PLayout, typename T, unsigned Dim>
FieldSolverType PartBunch<PLayout, T, Dim>::getFieldSolverType() const {
    if (fs_m) {
        return fs_m->getFieldSolverType();
    } else {
        return FieldSolverType::NONE;
    }
}
*/

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setStepsPerTurn(int n) {
    stepsPerTurn_m = n;
}

template <class PLayout, typename T, unsigned Dim>
int PartBunch<PLayout, T, Dim>::getStepsPerTurn() const {
    return stepsPerTurn_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setGlobalTrackStep(long long n) {
    globalTrackStep_m = n;
}

template <class PLayout, typename T, unsigned Dim>
long long PartBunch<PLayout, T, Dim>::getGlobalTrackStep() const {
    return globalTrackStep_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setLocalTrackStep(long long n) {
    localTrackStep_m = n;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::incTrackSteps() {
    globalTrackStep_m++;
    localTrackStep_m++;
}

template <class PLayout, typename T, unsigned Dim>
long long PartBunch<PLayout, T, Dim>::getLocalTrackStep() const {
    return localTrackStep_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setNumBunch(short n) {
    numBunch_m = n;
    bunchTotalNum_m.resize(n);
    bunchLocalNum_m.resize(n);
}

template <class PLayout, typename T, unsigned Dim>
short PartBunch<PLayout, T, Dim>::getNumBunch() const {
    return numBunch_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setGlobalMeanR(Vector_t<T, Dim> globalMeanR) {
    globalMeanR_m = globalMeanR;
}

template <class PLayout, typename T, unsigned Dim>
Vector_t<T, Dim> PartBunch<PLayout, T, Dim>::getGlobalMeanR() {
    return globalMeanR_m;
}

/*
template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion) {
    globalToLocalQuaternion_m = globalToLocalQuaternion;
}

template <class PLayout, typename T, unsigned Dim>
Quaternion_t PartBunch<PLayout, T, Dim>::getGlobalToLocalQuaternion() {
    return globalToLocalQuaternion_m;
}
*/
template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setSteptoLastInj(int n) {
    SteptoLastInj_m = n;
}

template <class PLayout, typename T, unsigned Dim>
int PartBunch<PLayout, T, Dim>::getSteptoLastInj() const {
    return SteptoLastInj_m;
}

/*
template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getQ() const {
    return reference->getQ();
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getM() const {
    return reference->getM();
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getP() const {
    return reference->getP();
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getE() const {
    return reference->getE();
}

template <class PLayout, typename T, unsigned Dim>
ParticleOrigin PartBunch<PLayout, T, Dim>::getPOrigin() const {
    return refPOrigin_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setPOrigin(ParticleOrigin origin) {
    refPOrigin_m = origin;
}

template <class PLayout, typename T, unsigned Dim>
ParticleType PartBunch<PLayout, T, Dim>::getPType() const {
    return refPType_m;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setPType(const std::string& type) {
    refPType_m = ParticleProperties::getParticleType(type);
}
*/

/*
template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::resetQ(double q) {
    const_cast<PartData*>(reference)->setQ(q);
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::resetM(double m) {
    const_cast<PartData*>(reference)->setM(m);
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getdE() const {
    return 0.0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getInitialBeta() const {
    return reference->getBeta();
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getInitialGamma() const {
    return reference->getGamma();
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getGamma(int i) {
    return 0;
}

template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::getBeta(int i) {
    return 0;
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::actT(){
    // do nothing here
};

template <class PLayout, typename T, unsigned Dim>
const PartData* PartBunch<PLayout, T, Dim>::getReference() const {
    return reference;
}


    / template <class PLayout, typename T, unsigned Dim>
    void PartBunch<PLayout, T, Dim>::calcEMean() {
    // momentsComputer_m.computeMeanKineticEnergy(*this);
}
*/

/*
  template <class PLayout, typename T, unsigned Dim>
Inform& PartBunch<PLayout, T, Dim>::print(Inform& os) {
    return PartBunch<PLayout, T, Dim>::print(os);
}
*/

template <class PLayout, typename T, unsigned Dim>
Inform& operator<<(Inform& os, PartBunch<PLayout, T, Dim>& p) {
    return p.print(os);
}

template <class PLayout, typename T, unsigned Dim>
Inform& PartBunch<PLayout, T, Dim>::print(Inform& os) {
    // if (this->getLocalNum() != 0) {  // to suppress Nans
    Inform::FmtFlags_t ff = os.flags();

    double pathLength = get_sPos();

    os << std::scientific;
    os << level1 << "\n";
    os << "* ************** B U N C H "
          "********************************************************* \n";
    os << "* NP              = " << this->getLocalNum() << "\n";
    /*
    os << "* Qtot            = " << std::setw(17) << Util::getChargeString(std::abs(sum(Q)))
       << "         "
       << "Qi    = " << std::setw(17) << Util::getChargeString(std::abs(qi_m)) << "\n";
    os << "* Ekin            = " << std::setw(17)
       << Util::getEnergyString(get_meanKineticEnergy()) << "         "
       << "dEkin = " << std::setw(17) << Util::getEnergyString(getdE()) << "\n";
    os << "* rmax            = " << Util::getLengthString(rmax_m, 5) << "\n";
    os << "* rmin            = " << Util::getLengthString(rmin_m, 5) << "\n";
    if (this->getTotalNum() >= 2) {  // to suppress Nans
        os << "* rms beam size   = " << Util::getLengthString(get_rrms(), 5) << "\n";
        os << "* rms momenta     = " << std::setw(12) << std::setprecision(5) << get_prms()
           << " [beta gamma]\n";
        os << "* mean position   = " << Util::getLengthString(get_rmean(), 5) << "\n";
        os << "* mean momenta    = " << std::setw(12) << std::setprecision(5) << get_pmean()
           << " [beta gamma]\n";
        os << "* rms emittance   = " << std::setw(12) << std::setprecision(5) << get_emit()
           << " (not normalized)\n";
        os << "* rms correlation = " << std::setw(12) << std::setprecision(5) << get_rprms()
           << "\n";
    }
    os << "* hr              = " << Util::getLengthString(get_hr(), 5) << "\n";
    os << "* dh              = " << std::setw(13) << std::setprecision(5) << dh_m * 100
       << " [%]\n";
    os << "* t               = " << std::setw(17) << Util::getTimeString(getT()) << "         "
       << "dT    = " << std::setw(17) << Util::getTimeString(getdT()) << "\n";
    os << "* spos            = " << std::setw(17) << Util::getLengthString(pathLength) << "\n";
    */
    os << "* "
          "********************************************************************************** "
       << endl;
    os.flags(ff);
    // }
    return os;
}

// angle range [0~2PI) degree
template <class PLayout, typename T, unsigned Dim>
double PartBunch<PLayout, T, Dim>::calculateAngle(double x, double y) {
    double thetaXY = atan2(y, x);

    return thetaXY >= 0 ? thetaXY : thetaXY;  // ADA + Physics::two_pi;
}

/*
 * Virtual member functions
 */

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::swap(unsigned int i, unsigned int j) {
    if (i >= this->getLocalNum() || j >= this->getLocalNum() || i == j)
        return;

    std::swap(this->R(i), this->R(j));
    std::swap(this->P(i), this->P(j));
    std::swap(this->Q(i), this->Q(j));
    std::swap(this->M(i), this->M(j));
    std::swap(this->Phi(i), this->Phi(j));
    /*
    std::swap(Ef(i), Ef[j]);
    std::swap(Eftmp(i), Eftmp[j]);
    std::swap(Bf(i), Bf[j]);
    std::swap(Bin(i), Bin[j]);
    std::swap(dt(i), dt[j]);
    std::swap(PType(i), PType[j]);
    std::swap(POrigin(i), POrigin[j]);
    std::swap(TriID(i), TriID[j]);
    std::swap(cavityGapCrossed(i), cavityGapCrossed[j]);
    std::swap(bunchNum(i), bunchNum[j]);
    */
}

template <class PLayout, typename T, unsigned Dim>
void PartBunch<PLayout, T, Dim>::setBeamFrequency(double f) {
    periodLength_m = 1.;  // ADA Physics::c / f;
}

// template <class PLayout, typename T, unsigned Dim>
// matrix_t PartBunch<PLayout, T, Dim>::getSigmaMatrix() const {}
