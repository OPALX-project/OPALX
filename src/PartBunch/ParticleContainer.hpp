
#ifndef OPAL_PARTICLE_CONTAINER_H
#define OPAL_PARTICLE_CONTAINER_H

// #include <functional>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include "Manager/BaseManager.h"

#include "PartBunch/ParticleContainerTypes.h"

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/PartData.h"
#include "Algorithms/Quaternion.hpp"
#include "PartBunch/BunchStateHandler.h"

#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include "Physics/Physics.h"

// ParticleSpatialOverlapLayout.hpp currently includes IPPL's Alpine example ParticleContainer,
// whose global class name collides with OPALX's container. The layout itself does not depend on
// that example type, so suppress only that application header while including the IPPL layout.
#ifndef IPPL_PARTICLE_CONTAINER_H
#define OPALX_SUPPRESS_IPPL_ALPINE_PARTICLE_CONTAINER
#define IPPL_PARTICLE_CONTAINER_H
#endif
#include "Particle/ParticleSpatialOverlapLayout.h"
#ifdef OPALX_SUPPRESS_IPPL_ALPINE_PARTICLE_CONTAINER
#undef IPPL_PARTICLE_CONTAINER_H
#undef OPALX_SUPPRESS_IPPL_ALPINE_PARTICLE_CONTAINER
#endif

#include "Steppers/BorisPusher.h"
#include "Steppers/SpinTBMTPusher.h"

// #include <Kokkos_Core.hpp>

template <typename T>
using ParticleAttrib = ippl::ParticleAttrib<T>;

using size_type = ippl::detail::size_type;

#include "Processes/GlobalProcesses/GlobalProcess.h"

/**
 * @class ParticleContainer
 * @brief Container for all per-particle (and per-simulation) fields tracked during OPALX tracking.
 *
 * The values tracked in Kokkos::Views during the simulation are:
 * R  - Position (from base class)
 * P  - Momentum [beta*gamma]
 * dt - Time step
 * Phi- Scalar potential
 * Bin- Energy bin
 * E  - Electric field
 * B  - Magnetic field
 *
 * Charge (Q) and mass (M):
 * - Default: QM_MODE="SINGLE" -> QMStorageMode=SingleValue
 *   * Q and M are stored as a single shared value per container (memory-efficient).
 *   * Access via `getQView()` / `getMView()` returns the shared views.
 * - Alternative: QM_MODE="ATTRIBUTES" -> QMStorageMode=Attributes
 *   * Q and M are stored as per-particle attributes.
 *   * Access via `getQView()` / `getMView()` returns the per-particle attribute views.
 *
 * Access to Q/M should be done with `getQView()` / `getMView()`.
 * They automatically select the correct underlying storage mode.
 */
template <typename T, unsigned Dim = 3>
class ParticleContainer
    : public ippl::ParticleBase<
              ippl::ParticleSpatialLayout<T, Dim, ippl::UniformCartesian<T, Dim>>,
              Kokkos::DefaultExecutionSpace::memory_space> {
    /**
     * @brief Alias for the `ippl::ParticleBase` specialization this container inherits from.
     *
     * The second template argument is a parameter pack of Kokkos view properties forwarded to
     * the optional `ID` attribute's storage. IPPL gates the `ID` attribute on
     * `sizeof...(IDProperties) > 0`, so passing any property here turns IDs ON; an empty pack
     * would leave them disabled. We pass `Kokkos::DefaultExecutionSpace::memory_space` so the
     * ID view lives in the same space as the rest of the bunch (host or device, matching the
     * build backend). With IDs enabled, IPPL auto-assigns globally unique `std::int64_t` IDs in
     * `Base::create()`.
     */
    using Base = ippl::ParticleBase<
            ippl::ParticleSpatialLayout<T, Dim, ippl::UniformCartesian<T, Dim>>,
            Kokkos::DefaultExecutionSpace::memory_space>;

private:
    /**
     * Forbid access of the following functions outside of the ParticleContainer wrappers! This is
     * for safety, since it might lead to undefined behaviour. The idea is to handle particle count
     * changes completely through the ParticleContainer!
     */
    using Base::alloc;
    using Base::create;
    using Base::destroy;

public:
    using SpatialLayout_t = ippl::ParticleSpatialLayout<T, Dim, ippl::UniformCartesian<T, Dim>>;
    using P3MLayout_t = ippl::ParticleSpatialOverlapLayout<T, Dim, ippl::UniformCartesian<T, Dim>>;

    enum class LayoutType { Spatial, SpatialOverlap };

    enum class QMStorageMode { SingleValue, Attributes };

    /// Defines which type to use as a particle bin.
    using bin_index_type = short int;  // Needed in AdaptBins class

    /// View types of Q and M values
    using qm_view_type = typename ippl::ParticleAttrib<double>::view_type;

    /// Per-particle polarization vector type: 3D vector in double precision, |Pol| in [0, 1].
    /// Double, like R and P, because Pol is rotated with every frame change (transformBunch)
    /// and float storage would round it again at each rotation.
    using spin_vector_type = ippl::Vector<double, 3>;

public:
    /// Charge view in [Cb].
    /// In `SingleValue` mode this is a rank-1 view of length 1.
    /// In `Attributes` mode this is the per-particle attribute view.
    qm_view_type getQView() const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            return QAttr.getView();
        }
        return QView_m;
    }

    /// Charge view in [Cb] over the full attribute storage.
    /// In `Attributes` mode it also covers the ghost particles that a spatial-overlap layout
    /// stores after the local ones, where `getQView()` ends at the local count.
    qm_view_type getQStorageView() const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            return qm_view_type(QAttr.getView().data(), QAttr.size());
        }
        return QView_m;
    }

    /// Mass view in [GeV].
    /// In `SingleValue` mode this is a rank-1 view of length 1.
    /// In `Attributes` mode this is the per-particle attribute view.
    qm_view_type getMView() const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            return MAttr.getView();
        }
        return MView_m;
    }

    /// timestep in [s]
    ippl::ParticleAttrib<double> dt;

    /// the scalar potential in [Cb/s]
    ippl::ParticleAttrib<double> Phi;

    /// the energy bin the particle is in
    ippl::ParticleAttrib<bin_index_type> Bin;

    /// the particle specis
    short Sp = 0;

    /// particle momenta [\beta\gamma]
    typename Base::particle_position_type P;

    /// electric field at particle position
    typename Base::particle_position_type E;

    /// electric field for gun simulation with bins
    // typename Base::particle_position_type Etmp; // TODO: might not need this...

    /// magnetic field at particle position
    typename Base::particle_position_type B;

    /// particle deletion mask (indicates which particles are deleted every timestep)
    ippl::ParticleAttrib<bool> InvalidMask;

    /// Per-particle polarization vector (rest-frame Pauli expectation values along the same
    /// axes as P, i.e. the bunch frame; |Pol| in [0, 1]). Registered only when spinEnabled_m
    /// is true.
    ippl::ParticleAttrib<spin_vector_type> Pol;

    /// Returns true when per-particle spin storage was enabled at construction
    /// (enabled per beam when the BEAM has POLARIZATION set).
    bool hasSpin() const { return spinEnabled_m; }

public:
    ParticleContainer(
            Mesh_t<Dim>& mesh, FieldLayout_t<Dim>& FL, bool spinEnabled = false,
            LayoutType layoutType = LayoutType::Spatial, T overlapCutoff = 0.0,
            ippl::BC particleBC = ippl::BC::PERIODIC)
        : spatialLayout_m(
                  layoutType == LayoutType::Spatial ? std::make_unique<SpatialLayout_t>(FL, mesh)
                                                    : nullptr),
          overlapLayout_m(
                  layoutType == LayoutType::SpatialOverlap
                          ? std::make_unique<P3MLayout_t>(FL, mesh, overlapCutoff)
                          : nullptr),
          qmStorageMode_m(
                  Options::useQMAttributes ? QMStorageMode::Attributes
                                           : QMStorageMode::SingleValue),
          QView_m("ParticleContainer::QView_m", 1),
          MView_m("ParticleContainer::MView_m", 1),
          spinEnabled_m(spinEnabled) {
        this->initialize(getPL());
        registerAttributes();
        setupBCs(particleBC);
        Kokkos::deep_copy(QView_m, 0.0);
        Kokkos::deep_copy(MView_m, 0.0);
    }

    ~ParticleContainer() {}

    void registerAttributes() {
        // register the particle attributes
        this->addAttribute(dt);
        this->addAttribute(Phi);
        this->addAttribute(Bin);
        this->addAttribute(P);
        this->addAttribute(E);
        this->addAttribute(B);
        this->addAttribute(InvalidMask);
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            this->addAttribute(QAttr);
            this->addAttribute(MAttr);
        }
        if (spinEnabled_m) {
            this->addAttribute(Pol);
        }
    }

    void setupBCs(ippl::BC particleBC) { this->setParticleBC(particleBC); }

    /// Frame the local particle data is in now.
    const ParticleFrame& frame() const { return containerState_m->frame; }

    /**
     * @brief Move the local particle data into @p target with one transform, or none if it is
     * already there.
     *
     * This is the only way particle data changes frame. It translates R and rotates P, E, B and
     * Pol; it is rank-local. Build @p target from the current reference frame (see
     * ParticleFrame).
     *
     * @param clearFields E and B are zeroed instead of rotated, also when no transform is needed.
     */
    void enterFrame(const ParticleFrame& target, bool clearFields = false) {
        if (clearFields) {
            this->E = Vector_t<T, Dim>(0.0);
            this->B = Vector_t<T, Dim>(0.0);
        }
        if (frame() == target) return;
        transformBunch(frame().to(target), clearFields);
        containerState_m->frame = target;
    }

    void enterReferenceFrame() { enterFrame(ParticleFrame()); }

    /// Holds the container in @p frame and returns it to its reference frame when the scope
    /// ends, also when an exception unwinds it.
    class FrameScope {
    public:
        FrameScope(
                ParticleContainer& container, const ParticleFrame& frame, bool clearFields = false)
            : container_m(container) {
            container_m.enterFrame(frame, clearFields);
        }
        FrameScope(const FrameScope&) = delete;
        ~FrameScope() { container_m.enterReferenceFrame(); }

    private:
        ParticleContainer& container_m;
    };

    /// @throws OpalException from @p where if the particle data is not in the reference frame.
    void requireReferenceFrame(const std::string& where) const {
        if (!frame().isReference())
            throw OpalException(
                    where, std::string("Particle data is in the ") + frame().name()
                                   + " frame, but the reference frame is required.");
    }

    /**
     * @brief Bring the local particles from index @p first on, written in @p createdIn, into
     * the current frame. The particles before @p first do not move.
     *
     * For emission: samplers append new particles in their own frame.
     */
    void adoptParticles(const ParticleFrame& createdIn, size_t first) {
        if (!(createdIn == frame())) transformBunch(createdIn.to(frame()), false, first);
    }
    SpatialLayout_t& getPL() {
        return overlapLayout_m ? static_cast<SpatialLayout_t&>(*overlapLayout_m) : *spatialLayout_m;
    }

    const SpatialLayout_t& getPL() const {
        return overlapLayout_m ? static_cast<const SpatialLayout_t&>(*overlapLayout_m)
                               : *spatialLayout_m;
    }

    bool hasP3MLayout() const { return overlapLayout_m != nullptr; }

    P3MLayout_t& getP3MLayout() {
        if (!overlapLayout_m) {
            throw OpalException(
                    "ParticleContainer::getP3MLayout",
                    "The particle container does not use ParticleSpatialOverlapLayout.");
        }
        return *overlapLayout_m;
    }

    const P3MLayout_t& getP3MLayout() const {
        if (!overlapLayout_m) {
            throw OpalException(
                    "ParticleContainer::getP3MLayout",
                    "The particle container does not use ParticleSpatialOverlapLayout.");
        }
        return *overlapLayout_m;
    }

    void updateLayout(FieldLayout_t<Dim>& FL, Mesh_t<Dim>& mesh) {
        if (overlapLayout_m) {
            overlapLayout_m->updateLayout(FL, mesh);
        } else {
            spatialLayout_m->updateLayout(FL, mesh);
        }
    }

    /**
     * @brief Update the overlap layout and its P3M cutoff together.
     * @param overlapCutoff Cutoff in metres, checked against the new rank regions; it takes
     * effect at the next update().
     */
    void updateLayout(FieldLayout_t<Dim>& FL, Mesh_t<Dim>& mesh, T overlapCutoff) {
        getP3MLayout().updateLayout(FL, mesh, overlapCutoff);
    }

    void update() {
        if (overlapLayout_m) {
            overlapLayout_m->update(*this);
        } else {
            spatialLayout_m->update(*this);
        }
    }

    void setBunchStateHandler(std::shared_ptr<BunchStateHandler> handler) {
        // We only keep the slot: per-container flags own their own sync, so
        // no back-reference to the handler is needed. The handler itself
        // manages bunch-wide state, which ParticleContainer never touches.
        // The new slot takes over the state held so far.
        const auto held   = *containerState_m;
        containerState_m  = handler->registerContainer();
        *containerState_m = held;
    }

    bool isUnitlessPositions() const { return containerState_m->unitlessPositions; }

    /**
     * @brief Set particle charge for the active Q storage mode.
     * @param q Charge value in [Cb].
     *
     * In `QMStorageMode::Attributes`, this assigns `q` to every local particle.
     * In `QMStorageMode::SingleValue`, this updates the shared scalar charge view.
     */
    void setQ(double q) {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view         = QAttr.getView();
            const size_type n = this->getLocalNum();
            if (n == 0) {
                return;
            }
            Kokkos::parallel_for(
                    "ParticleContainer::setQ", n,
                    KOKKOS_LAMBDA(const size_type i) { view(i) = q; });
            Kokkos::fence();
        } else {
            Kokkos::deep_copy(QView_m, q);
        }
    }

    /// @brief Get charge per particle [Cb].
    double getChargePerParticle() const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view = QAttr.getView();
            if (view.extent(0) == 0) {
                return 0.0;
            }
            double q = 0.0;
            Kokkos::deep_copy(q, Kokkos::subview(view, 0));
            return q;
        }
        double q = 0.0;
        Kokkos::deep_copy(q, Kokkos::subview(QView_m, 0));
        return q;
    }

    /// @brief Get total charge [Cb] in this container.
    double getTotalCharge() const { return getChargePerParticle() * this->getTotalNum(); }

    /**
     * @brief Set particle mass for the active M storage mode.
     * @param m Mass value in [GeV].
     *
     * In `QMStorageMode::Attributes`, this assigns `m` to every local particle.
     * In `QMStorageMode::SingleValue`, this updates the shared scalar mass view.
     */
    void setM(double m) {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view         = MAttr.getView();
            const size_type n = view.extent(0);

            Kokkos::parallel_for(
                    "ParticleContainer::setM", n,
                    KOKKOS_LAMBDA(const size_type i) { view(i) = m; });
            Kokkos::fence();
        } else {
            Kokkos::deep_copy(MView_m, m);
        }
    }

    /// @brief Get mass per particle [GeV].
    double getMassPerParticle() const {
        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto view = MAttr.getView();
            if (view.extent(0) == 0) {
                return 0.0;
            }
            double m = 0.0;
            Kokkos::deep_copy(m, Kokkos::subview(view, 0));
            return m;
        }
        double m = 0.0;
        Kokkos::deep_copy(m, Kokkos::subview(MView_m, 0));
        return m;
    }

    /// @brief Get total mass [GeV] in this container.
    double getTotalMass() const { return getMassPerParticle() * this->getTotalNum(); }

    /// @brief Get the reference particle position (const).
    const Vector_t<double, Dim>& getRefPartR() const { return refPartR_m; }

    /// @brief Get the reference particle position.
    Vector_t<double, Dim>& getRefPartR() { return refPartR_m; }

    /// @brief Set the reference particle position.
    void setRefPartR(const Vector_t<double, Dim>& refPartR) { refPartR_m = refPartR; }

    /// @brief Get the reference particle momentum (const).
    const Vector_t<double, Dim>& getRefPartP() const { return refPartP_m; }

    /// @brief Get the reference particle momentum.
    Vector_t<double, Dim>& getRefPartP() { return refPartP_m; }

    /// @brief Set the reference particle momentum.
    void setRefPartP(const Vector_t<double, Dim>& refPartP) { refPartP_m = refPartP; }

    /// @brief Set reference particle data; the container keeps its own copy.
    void setReference(const PartData& ref) { reference_m = ref; }

    /// @brief Get reference particle data; nullptr until setReference() is called.
    const PartData* getReference() const { return reference_m ? &*reference_m : nullptr; }

    /// @brief Set longitudinal position along design trajectory.
    void set_sPos(double sPos) { sPos_m = sPos; }

    /// @brief Get longitudinal position along design trajectory.
    double get_sPos() const { return sPos_m; }

    /// @brief Get local-to-lab coordinate transformation (const).
    const CoordinateSystemTrafo& getToLabTrafo() const { return toLabTrafo_m; }

    /// @brief Redefine the reference frame; the particle data is declared to be in it.
    void setToLabTrafo(const CoordinateSystemTrafo& toLabTrafo) {
        requireReferenceFrame("ParticleContainer::setToLabTrafo");
        toLabTrafo_m = toLabTrafo;
    }

    /**
     * @brief Advance reference/lab transform state and map bunch accordingly.
     * @return The frame change applied to the bunch (old bunch frame to new bunch frame), so
     *   containers that use this container's reference can apply the same change.
     */
    CoordinateSystemTrafo updateRefToLabCSTrafo(double bunchDT) {
        requireReferenceFrame("ParticleContainer::updateRefToLabCSTrafo");
        Vector_t<double, 3> R = toLabTrafo_m.transformFrom(refPartR_m);
        Vector_t<double, 3> P = toLabTrafo_m.rotateFrom(refPartP_m);

        const double ds =
                std::copysign(1.0, bunchDT) * std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
        sPos_m += ds;

        CoordinateSystemTrafo update(R, getQuaternion(P, Vector_t<double, 3>(0, 0, 1)));
        transformBunch(update);
        toLabTrafo_m = toLabTrafo_m * update.inverted();
        return update;
    }

    /**
     * @brief Take over another container's reference particle, frame and path length.
     *
     * Used for a decay daughter container, which has no reference particle of its own. The
     * particles are not moved; use this only when both containers already share a frame or
     * this container is empty.
     *
     * @param root Container whose reference particle this container uses.
     */
    void copyReferenceOf(const ParticleContainer& root) {
        requireReferenceFrame("ParticleContainer::copyReferenceOf");
        refPartR_m   = root.refPartR_m;
        refPartP_m   = root.refPartP_m;
        toLabTrafo_m = root.toLabTrafo_m;
        sPos_m       = root.sPos_m;
    }

    /**
     * @brief Follow a frame change of the container whose reference particle this one uses.
     *
     * Applies @p update (the value returned by @p root's updateRefToLabCSTrafo, or the
     * identity) to this container's particles, then copies @p root's reference particle,
     * frame and path length. The particles' lab positions do not change.
     *
     * @param root Container whose reference particle this container uses.
     * @param update Frame change @p root applied to its own particles in this step.
     */
    void followFrameOf(const ParticleContainer& root, const CoordinateSystemTrafo& update) {
        requireReferenceFrame("ParticleContainer::followFrameOf");
        transformBunch(update);
        copyReferenceOf(root);
    }

    /// @brief Apply a fractional Boris step and update reference/lab transform state.
    void applyFractionalStep(double tau, double pathLengthTarget) {
        requireReferenceFrame("ParticleContainer::applyFractionalStep");
        const BorisPusher pusher;
        refPartR_m /= (Physics::c * 2 * tau);
        pusher.push(refPartR_m, refPartP_m, tau);
        refPartR_m *= (Physics::c * 2 * tau);

        sPos_m = pathLengthTarget;
        CoordinateSystemTrafo update(
                refPartR_m, getQuaternion(refPartP_m, Vector_t<double, 3>(0, 0, 1)));
        toLabTrafo_m = toLabTrafo_m * update.inverted();
    }

    /// @brief Zero E and B of the local particles.
    void zeroFields() {
        E = 0;
        B = 0;
    }

    /**
     * @brief Drift every local particle over half its own dt, R += (dt / 2) v, on unitless
     *        positions R / (c dt).
     */
    void pushParticles() {
        // Per-particle dt is used so that newly emitted particles with fractional dt
        // (sampled during emission) are pushed proportionally to their sub-timestep fraction.
        switchToUnitlessPositions();

        const BorisPusher pusher;
        auto Rview             = this->R.getView();
        auto Pview             = P.getView();
        const size_type nLocal = this->getLocalNum();
        Kokkos::parallel_for(
                "ParticleContainer::pushParticles", nLocal, KOKKOS_LAMBDA(const size_type i) {
                    // Half drift: x_{n+1/2} = x_n + (dt/2) * v; unitless form via
                    // pusher.push(..., 0).
                    // TODO: verify sign convention for half-step push.
                    Vector_t<double, 3> x = Rview(i);
                    pusher.push(
                            x, Pview(i),
                            0);  // this 0 is "dt" that is not used with unitless positions!
                    Rview(i) = x;
                });

        switchOffUnitlessPositions();
        Kokkos::fence();
    }

    /**
     * @brief Boris kick of every local particle over its own dt, with E and B at the particle
     *        and the mass and charge of the reference particle.
     */
    void kickParticles() {
        Inform m("ParticleContainer::kickParticles");

        auto Pview  = P.getView();
        auto dtview = dt.getView();
        auto Efview = E.getView();
        auto Bfview = B.getView();
        m << level5 << "Got particle views for kick operation." << endl;
        // Mass (eV) and charge (proton charges) from this container's reference particle,
        // passed explicitly into BorisPusher::kick for GPU-safe kernels.

        const PartData& ref    = *getReference();
        const double mass      = ref.getM();
        const double charge    = ref.getQ();
        const size_type nLocal = this->getLocalNum();
        const BorisPusher pusher;
        Kokkos::parallel_for(
                "ParticleContainer::kickParticles", nLocal, KOKKOS_LAMBDA(const size_type i) {
                    // Boris kick: Birdsall & Langdon (1985), ch. 4-4 (non-relativistic) and ch.
                    // 15-4 (relativistic); scaled P = (v/c)*gamma, mass in rest energy units.
                    Vector_t<double, 3> p = Pview(i);
                    // TODO: consider dropping unused R/dt arguments from kick when API allows.
                    pusher.kick(0, p, Efview(i), Bfview(i), dtview(i), mass, charge);
                    Pview(i) = p;
                });
        Kokkos::fence();

        m << level5 << "Completed parallel kick operation." << endl;
    }

    /**
     * @brief Rotate Pol of every local particle by half its own dt (Thomas-BMT), with the
     *        current P, E and B. P, E, B and Pol are all along the bunch-frame axes.
     */
    void evolveSpinHalfStep() {
        const PartData& ref = *getReference();
        const double mass   = ref.getM();
        const double charge = ref.getQ();
        const double anom   = ref.getAnomaly();

        auto Polview           = Pol.getView();
        auto Pview             = P.getView();
        auto Efview            = E.getView();
        auto Bfview            = B.getView();
        auto dtview            = dt.getView();
        const size_type nLocal = this->getLocalNum();

        SpinTBMTPusher spinPusher;
        Kokkos::parallel_for(
                "ParticleContainer::evolveSpinHalfStep", nLocal, KOKKOS_LAMBDA(const size_type j) {
                    spinPusher.evolve(
                            Polview(j), Pview(j), Efview(j), Bfview(j), 0.5 * dtview(j), mass,
                            charge, anom);
                });
        Kokkos::fence();
    }

    /**
     * @brief Scale particle time-step weights by charge before scatter.
     *
     * Multiplies each local `dt[i]` by the corresponding charge used for deposition.
     * In `QMStorageMode::Attributes`, the per-particle `QAttr(i)` is used.
     * In `QMStorageMode::SingleValue`, the shared scalar `QView_m(0)` is used.
     *
     */
    void scaleDtByCharge() {
        auto dtView       = dt.getView();
        const size_type n = this->getLocalNum();
        if (n == 0) {
            return;
        }

        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto qView = QAttr.getView();
            Kokkos::parallel_for(
                    "ParticleContainer::scaleDtByCharge(attrs)", n,
                    KOKKOS_LAMBDA(const size_type i) { dtView(i) *= qView(i); });
        } else {
            auto QView = QView_m;
            Kokkos::parallel_for(
                    "ParticleContainer::scaleDtByCharge(single)", n,
                    KOKKOS_LAMBDA(const size_type i) { dtView(i) *= QView(0); });
        }
        Kokkos::fence();
    }

    /**
     * @brief Restore original `dt` values after `scaleDtByCharge()`.
     *
     * Divides each local `dt[i]` by the same charge factor applied in
     * `scaleDtByCharge()`.
     *
     */
    void unscaleDtByCharge() {
        auto dtView       = dt.getView();
        const size_type n = this->getLocalNum();
        if (n == 0) {
            return;
        }

        if (qmStorageMode_m == QMStorageMode::Attributes) {
            auto qView = QAttr.getView();
            Kokkos::parallel_for(
                    "ParticleContainer::unscaleDtByCharge(attrs)", n,
                    KOKKOS_LAMBDA(const size_type i) { dtView(i) /= qView(i); });
        } else {
            auto QView = QView_m;
            Kokkos::parallel_for(
                    "ParticleContainer::unscaleDtByCharge(single)", n,
                    KOKKOS_LAMBDA(const size_type i) { dtView(i) /= QView(0); });
        }
        Kokkos::fence();
    }

    /**
     * @brief Transform positions to unitless coordinates using each particle's dt[i].
     *
     * Applies \f$ R'_i = R_i / (c \, dt_i) \f$. Requires valid non-zero dt values per particle.
     *
     * @throws OpalException if this container is already in unitless positions.
     */
    void switchToUnitlessPositions() {
        if (containerState_m->unitlessPositions) {
            throw OpalException(
                    "ParticleContainer::switchToUnitlessPositions",
                    "ParticleContainer is already in unitless positions!");
        }
        auto Rview             = this->R.getView();
        auto dtview            = this->dt.getView();
        const size_type nLocal = this->getLocalNum();
        Kokkos::parallel_for(
                "ParticleContainer::switchToUnitlessPositions", nLocal,
                KOKKOS_LAMBDA(const size_type i) { Rview(i) *= 1.0 / (Physics::c * dtview(i)); });
        Kokkos::fence();
        containerState_m->setUnitlessPositions(true);
    }

    /**
     * @brief Restore physical positions from unitless form using each particle's dt[i].
     *
     * Applies \f$ R_i = R'_i \, c \, dt_i \f$.
     *
     * @throws OpalException if this container is not currently in unitless positions.
     */
    void switchOffUnitlessPositions() {
        if (!containerState_m->unitlessPositions) {
            throw OpalException(
                    "ParticleContainer::switchOffUnitlessPositions",
                    "ParticleContainer is already in physical positions!");
        }
        auto Rview             = this->R.getView();
        auto dtview            = this->dt.getView();
        const size_type nLocal = this->getLocalNum();
        Kokkos::parallel_for(
                "ParticleContainer::switchOffUnitlessPositions", nLocal,
                KOKKOS_LAMBDA(const size_type i) { Rview(i) *= Physics::c * dtview(i); });
        Kokkos::fence();
        containerState_m->setUnitlessPositions(false);
    }
    QMStorageMode getQMStorageMode() const { return qmStorageMode_m; }

    void setGlobalProcesses(std::vector<std::unique_ptr<GlobalProcess>> processes) {
        globalProcesses_m = std::move(processes);
    }

    const std::vector<std::unique_ptr<GlobalProcess>>& getGlobalProcesses() const {
        return globalProcesses_m;
    }

    /**
     * @brief Mark particles whose position lies outside the box [@p lower, @p upper].
     *
     * ORs the outside decision into InvalidMask; deletion is deferred to deleteInvalidParticles().
     * Collective.
     *
     * @return Global number of newly marked particles (across all MPI ranks).
     */
    size_type markParticlesOutside(
            const Vector_t<double, Dim>& lower, const Vector_t<double, Dim>& upper) {
        const size_type nLocal = this->getLocalNum();
        const double lb0 = lower[0], lb1 = lower[1], lb2 = lower[2];
        const double ub0 = upper[0], ub1 = upper[1], ub2 = upper[2];

        auto invalid = InvalidMask.getView();
        auto Rview   = this->R.getView();

        size_type localMarkedNum = 0;
        Kokkos::parallel_reduce(
                "ParticleContainer::markParticlesOutside", nLocal,
                KOKKOS_LAMBDA(const size_type i, size_type& count) {
                    bool outside = (Rview(i)[0] < lb0 || Rview(i)[0] > ub0)
                                   || (Rview(i)[1] < lb1 || Rview(i)[1] > ub1)
                                   || (Rview(i)[2] < lb2 || Rview(i)[2] > ub2);
                    const bool newlyMarked = outside && !invalid(i);
                    invalid(i)             = invalid(i) || outside;
                    count += newlyMarked ? 1 : 0;
                },
                localMarkedNum);
        Kokkos::fence();  // not needed, but also doesn't hurt

        size_type globalMarkedNum = 0;
        ippl::Comm->allreduce(localMarkedNum, globalMarkedNum, 1, std::plus<size_type>());
        return globalMarkedNum;
    }

    /**
     * @brief Mark particles at or behind the plane z = @p planeZ that move backwards (P_z < 0),
     *        with z and P_z read in the frame @p toPlaneFrame maps into.
     *
     * ORs into InvalidMask; the particles stay where they are.
     *
     * @param toPlaneFrame Transform from the particles' current frame into the plane's frame.
     * @return Local number of newly marked particles (not reduced over MPI ranks).
     */
    size_type markBackwardParticles(const RigidTransform& toPlaneFrame, double planeZ) {
        auto Rview   = this->R.getView();
        auto Pview   = P.getView();
        auto invalid = InvalidMask.getView();

        size_type localMarked  = 0;
        const size_type nLocal = static_cast<size_type>(this->getLocalNum());
        Kokkos::parallel_reduce(
                "ParticleContainer::markBackwardParticles", nLocal,
                KOKKOS_LAMBDA(const size_type i, size_type& count) {
                    const Vector_t<double, 3> localR = toPlaneFrame.pointTo(Rview(i));
                    const Vector_t<double, 3> localP = toPlaneFrame.vectorTo(Pview(i));
                    const bool backwards             = localR[2] <= planeZ && localP[2] < 0.0;
                    const bool newlyMarked           = backwards && !invalid(i);
                    invalid(i)                       = invalid(i) || backwards;
                    count += newlyMarked ? 1 : 0;
                },
                localMarked);
        Kokkos::fence();
        return localMarked;
    }

    /**
     * @brief Create/allocate a specified number of particles.
     *
     * This function creates a given number of particles in the container. It's a wrapper around the
     * non destructive IPPL particle create function, but will print out a warning if the create
     * call led to unnecessary reallocation (i.e. if the new total number of particles exceeds the
     * previous capacity).
     *
     * @note The underlying `create` is a collective call, so all MPI ranks must call this function.
     * IPPL automatically handles the short-circuit if internal capacity is sufficient.
     *
     * @param numParticles The number of particles to create.
     */
    void createParticles(size_type numParticles) {
        Inform m("ParticleContainer::createParticles");

        // Total allocated capacity of the underlying view
        size_type oldCapacity = this->R.size();
        size_type oldLocalNum = this->getLocalNum();
        this->create(numParticles, true);  // non_destructive = true
        size_type newCapacity = this->R.size();

        // Initialize newly active slots.  Emitted particles are created after
        // the per-step field reset, so E/B must be explicitly cleared here.
        /// \todo: can probably be removed later, after my flattop debugging
        auto invalid = InvalidMask.getView();
        auto dtView  = dt.getView();
        auto phiView = Phi.getView();
        auto binView = Bin.getView();
        auto eView   = E.getView();
        auto bView   = B.getView();
        Kokkos::parallel_for(
                "ParticleContainer::createParticles::initializeNewSlots", numParticles,
                KOKKOS_LAMBDA(const size_type i) {
                    const size_type idx = oldLocalNum + i;
                    invalid(idx)        = false;
                    dtView(idx)         = 0.0;
                    phiView(idx)        = 0.0;
                    binView(idx)        = 0;
                    eView(idx)          = Vector_t<double, Dim>(0.0);
                    bView(idx)          = Vector_t<double, Dim>(0.0);
                });
        // A reused slot still holds the spin of a deleted particle. Pol exists only with spin.
        if (spinEnabled_m) {
            auto polView = Pol.getView();
            Kokkos::parallel_for(
                    "ParticleContainer::createParticles::resetPol", numParticles,
                    KOKKOS_LAMBDA(const size_type i) {
                        polView(oldLocalNum + i) = spin_vector_type(0.0);
                    });
        }
        Kokkos::fence();

        // Pretty print numParticles, newCapacity and new totalNum + localNum after creation
        constexpr int labelWidth = 32;
        m << level4 << std::left << std::setw(labelWidth) << "Requested creation:" << numParticles
          << " particles" << endl;
        m << level4 << std::setw(labelWidth) << "New total number:" << this->getTotalNum()
          << " (local: " << this->getLocalNum() << ")" << endl;
        m << level4 << std::setw(labelWidth) << "Underlying view capacity:" << newCapacity << endl;

        if (newCapacity != oldCapacity) {
            m << level1
              << "WARNING: createParticles triggered a reallocation of the underlying particle "
                 "views! This can be a costly operation. To avoid this, consider increasing "
                 "preallocation (BEAM::NALLOC) or the overallocation factor."
              << endl;
        }
    }

    void allocateParticles(size_type numParticles) {
        Inform m("ParticleContainer::allocateParticles");

        // Total allocated capacity of the underlying view
        size_type oldCapacity = this->R.size();
        if (oldCapacity != 0) {
            throw OpalException(
                    "ParticleContainer::allocateParticles",
                    "Underlying views already allocated. This function is meant to be called on an "
                    "empty container, since it is destructive on existing particles. If you want "
                    "to create particles without deallocating existing ones, use createParticles() "
                    "instead.");
        }

        this->alloc(numParticles);  // alloc is always destructive

        m << level4 << std::left << std::setw(32) << "Requested allocation:" << numParticles
          << " particles" << endl;
        m << level4 << std::setw(32) << "Size of underlying view:" << this->R.size() << endl;
    }

    void printRankLoadInfo(const std::string& label = "") const {
        const int nranks = ippl::Comm->size();
        const int rank   = ippl::Comm->rank();
        const int root   = 0;

        std::vector<size_type> localParticles(static_cast<size_t>(nranks), 0);
        std::vector<size_type> localCapacity(static_cast<size_t>(nranks), 0);
        std::vector<size_type> rankParticles(static_cast<size_t>(nranks), 0);
        std::vector<size_type> rankCapacity(static_cast<size_t>(nranks), 0);

        localParticles[static_cast<size_t>(rank)] = this->getLocalNum();
        localCapacity[static_cast<size_t>(rank)]  = this->R.size();

        ippl::Comm->reduce(
                localParticles.data(), rankParticles.data(), nranks, std::plus<size_type>(), root);
        ippl::Comm->reduce(
                localCapacity.data(), rankCapacity.data(), nranks, std::plus<size_type>(), root);

        if (rank != root) {
            return;
        }

        const size_type totalParticles =
                std::accumulate(rankParticles.begin(), rankParticles.end(), size_type{0});
        const size_type totalCapacity =
                std::accumulate(rankCapacity.begin(), rankCapacity.end(), size_type{0});

        auto percent = [](size_type value, size_type total) {
            return total > 0 ? 100.0 * static_cast<double>(value) / static_cast<double>(total)
                             : 0.0;
        };
        auto localRatioPercent = [](size_type numerator, size_type denominator) {
            return denominator > 0 ? 100.0 * static_cast<double>(numerator)
                                             / static_cast<double>(denominator)
                                   : 0.0;
        };

        Inform m("ParticleContainer::printRankLoadInfo", root);
        constexpr int labelWidth = 24;
        constexpr int colWidth   = 12;

        m << level2 << "Particle load by rank";
        if (!label.empty()) {
            m << " (" << label << ")";
        }
        m << endl;
        m << level2 << std::left << std::setw(labelWidth) << "Metric";
        for (int r = 0; r < nranks; ++r) {
            m << " | " << std::right << std::setw(colWidth) << ("Rank " + std::to_string(r));
        }
        m << endl;

        m << level2 << std::left << std::setw(labelWidth) << "Particles [%]" << std::fixed
          << std::setprecision(2);
        for (int r = 0; r < nranks; ++r) {
            m << " | " << std::right << std::setw(colWidth)
              << percent(rankParticles[static_cast<size_t>(r)], totalParticles);
        }
        m << endl;

        m << level2 << std::left << std::setw(labelWidth) << "Allocated memory [%]" << std::fixed
          << std::setprecision(2);
        for (int r = 0; r < nranks; ++r) {
            m << " | " << std::right << std::setw(colWidth)
              << percent(rankCapacity[static_cast<size_t>(r)], totalCapacity);
        }
        m << endl;

        m << level2 << std::left << std::setw(labelWidth) << "Overallocation [%]" << std::fixed
          << std::setprecision(2);
        for (int r = 0; r < nranks; ++r) {
            const auto i = static_cast<size_t>(r);
            m << " | " << std::right << std::setw(colWidth)
              << localRatioPercent(rankCapacity[i], rankParticles[i]);
        }
        m << endl;
    }

    /**
     * @brief Delete particles currently marked in InvalidMask.
     *
     * This is the only ParticleContainer function that is allowed to compact the
     * IPPL particle arrays. All deletion producers must only update InvalidMask.
     *
     * @return Global number of particles deleted (across all MPI ranks).
     */
    size_type deleteInvalidParticles() {
        Inform m("ParticleContainer::deleteInvalidParticles");

        const size_type nLocal = this->getLocalNum();
        auto invalid           = InvalidMask.getView();
        if (invalid.extent(0) < nLocal) {
            throw OpalException(
                    "ParticleContainer::deleteInvalidParticles",
                    "InvalidMask extent (" + std::to_string(invalid.extent(0))
                            + ") is smaller than local particle count (" + std::to_string(nLocal)
                            + ").");
        }

        size_type localDestroyNum = 0;
        Kokkos::parallel_reduce(
                "ParticleContainer::deleteInvalidParticles::count", nLocal,
                KOKKOS_LAMBDA(const size_type i, size_type& count) { count += invalid(i) ? 1 : 0; },
                localDestroyNum);
        Kokkos::fence();

        size_type globalDestroyNum = 0;
        ippl::Comm->allreduce(localDestroyNum, globalDestroyNum, 1, std::plus<size_type>());

        // This is a collective call! All ranks must execute it. Note that this is safe with the
        // current IPPL implementation: ParticleBase::internalDestroy() first consumes the invalid
        // view to build internal deleteIndex_m and keepIndex_m arrays. Meaning, the fact that
        // destroy edits the InvalidMask while using it is not a problem.
        Base::destroy(invalid, localDestroyNum);

        // Reset particle mask after deletion.
        InvalidMask = false;
        Kokkos::fence();

        constexpr int labelWidth = 32;
        m << level4 << std::left << std::setw(labelWidth)
          << "Requested destruction:" << localDestroyNum << " particles" << endl
          << std::setw(labelWidth) << "New total number:" << this->getTotalNum()
          << " (local: " << this->getLocalNum() << ")" << endl;

        return globalDestroyNum;
    }

private:
    /// Translate R and rotate P, E, B and Pol of the local particles from index @p first on;
    /// see enterFrame().
    void transformBunch(
            const CoordinateSystemTrafo& trafo, bool fieldsCleared = false, size_t first = 0) {
        const size_t nLoc = this->getLocalNum();
        if (first >= nLoc) return;
        const auto range = std::make_pair(first, nLoc);
        const size_t n   = nLoc - first;
        trafo.transformBunchTo(Kokkos::subview(this->R.getView(), range), n);
        trafo.rotateBunchTo(Kokkos::subview(this->P.getView(), range), n);
        if (!fieldsCleared) {
            trafo.rotateBunchTo(Kokkos::subview(this->E.getView(), range), n);
            trafo.rotateBunchTo(Kokkos::subview(this->B.getView(), range), n);
        }
        if (spinEnabled_m) {
            trafo.rotateBunchTo(Kokkos::subview(this->Pol.getView(), range), n);
        }
    }

    std::unique_ptr<SpatialLayout_t> spatialLayout_m;
    std::unique_ptr<P3MLayout_t> overlapLayout_m;

    QMStorageMode qmStorageMode_m = QMStorageMode::SingleValue;

    /// Per-container state slot; standalone until `setBunchStateHandler` replaces it with a slot
    /// allocated by the handler. Owned here as the only strong reference; the handler keeps a
    /// weak_ptr, so destroying this container automatically releases the slot. The slot's own
    /// methods handle MPI consistency, so no direct handler reference is needed.
    std::shared_ptr<BunchStateHandler::ContainerState> containerState_m =
            std::make_shared<BunchStateHandler::ContainerState>();

    // Single shared scalar mode stored as a length-1 Kokkos view.
    qm_view_type QView_m;
    qm_view_type MView_m;

    // Per-particle attributes mode
    ippl::ParticleAttrib<double> QAttr;
    ippl::ParticleAttrib<double> MAttr;

    // Whether the spin attribute is registered on this container.
    bool spinEnabled_m = false;

    // Reference particle information
    Vector_t<double, Dim> refPartR_m;
    Vector_t<double, Dim> refPartP_m;

    // Reference particle to lab transformation
    CoordinateSystemTrafo toLabTrafo_m;

    // Particle reference data (!= reference particle)
    std::optional<PartData> reference_m;

    // Distance along the beamline
    double sPos_m = 0.0;

    /// Global physics processes attached to this container.
    std::vector<std::unique_ptr<GlobalProcess>> globalProcesses_m;
};

#endif
