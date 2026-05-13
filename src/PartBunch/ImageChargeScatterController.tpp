template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterScaledDtAll(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho) const {
    Inform m("ImageChargeScatterController::scatterScaledDtAll");
    const size_t nLocal = pc->getLocalNum();
    m << level5 << "all-local scatter entry: localP=" << nLocal << endl;
    pc->scaleDtByCharge();
    m << level5 << "all-local scatter: dt scaled by charge." << endl;

    using view_type    = typename RhoField_t::view_type;
    view_type rhoView  = rho.getView();
    auto dtView        = pc->dt.getView();
    auto rView         = positions.getView();
    const auto& mesh   = rho.get_mesh();
    const auto& dx     = mesh.getMeshSpacing();
    const auto& origin = mesh.getOrigin();
    const auto invdx   = 1.0 / dx;
    const auto& layout = rho.getLayout();
    const auto& lDom   = layout.getLocalNDIndex();
    const int nghost   = rho.getNghost();

    int outOfBounds = 0;
    Kokkos::parallel_reduce(
            "ImageChargeScatterController::validateScatterBounds", nLocal,
            KOKKOS_LAMBDA(const size_t i, int& update) {
                const auto l = (rView(i) - origin) * invdx + 0.5;
                ippl::Vector<int, Dim> index = l;
                ippl::Vector<int, Dim> args  = index - lDom.first() + nghost;
                bool inBounds                = true;
                for (unsigned d = 0; d < Dim; ++d) {
                    inBounds = inBounds && args[d] > 0
                               && args[d] < static_cast<int>(rhoView.extent(d));
                }
                if (!inBounds) {
                    update += 1;
                }
            },
            outOfBounds);
    m << level5 << "all-local scatter: out-of-bounds CIC particles=" << outOfBounds << endl;

    m << level5 << "all-local scatter: CIC kernel start." << endl;
    Kokkos::parallel_for(
            "ImageChargeScatterController::scatterScaledDtAllCIC", nLocal,
            KOKKOS_LAMBDA(const size_t i) {
                const auto l = (rView(i) - origin) * invdx + 0.5;
                ippl::Vector<int, Dim> index = l;
                ippl::Vector<T, Dim> whi     = l - index;
                ippl::Vector<T, Dim> wlo     = 1.0 - whi;

                ippl::Vector<int, Dim> args = index - lDom.first() + nghost;
                bool inBounds               = true;
                for (unsigned d = 0; d < Dim; ++d) {
                    // CIC touches args[d] and args[d] - 1, so valid args are
                    // [1, extent - 1]. Anything outside would underflow or
                    // overrun the field view on device.
                    inBounds = inBounds && args[d] > 0
                               && args[d] < static_cast<int>(rhoView.extent(d));
                }
                if (inBounds) {
                    ippl::Vector<size_t, Dim> viewArgs = args;
                    ippl::detail::scatterToField(
                            std::make_index_sequence<1 << Dim>{}, rhoView, wlo, whi, viewArgs,
                            dtView(i));
                }
            });
    Kokkos::fence();
    m << level5 << "all-local scatter: CIC kernel done." << endl;

    m << level5 << "all-local scatter: accumulateHalo start." << endl;
    rho.accumulateHalo();
    m << level5 << "all-local scatter: accumulateHalo done." << endl;

    pc->unscaleDtByCharge();
    m << level5 << "all-local scatter: dt restored." << endl;
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterScaledDtSubset(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho,
        const BinPolicy_t& policy, const Hash_t& hash) const {
    Inform m("ImageChargeScatterController::scatterScaledDtSubset");
    m << level5 << "subset scatter entry: localP=" << pc->getLocalNum() << ", policy=["
      << policy.begin() << "," << policy.end() << "), hashExtent=" << hash.extent(0) << endl;
    pc->scaleDtByCharge();
    m << level5 << "subset scatter: dt scaled by charge." << endl;
    ippl::ParticleAttrib<T>* dtAttrib = &pc->dt;
    scatter(*dtAttrib, rho, positions, policy, hash);
    m << level5 << "subset scatter: scatter+halo done." << endl;
    pc->unscaleDtByCharge();
    m << level5 << "subset scatter: dt restored." << endl;
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterPrimaryOnly(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho) const {
    scatterScaledDtAll(pc, positions, rho);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterPrimaryOnly(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho,
        const BinPolicy_t& policy, const Hash_t& hash) const {
    scatterScaledDtSubset(pc, positions, rho, policy, hash);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterImageOnly(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho) const {
    if (!enabled_m) {
        return;
    }
    applyMirrorTransformAll(pc, positions);
    scatterScaledDtAll(pc, positions, rho);
    restoreMirrorTransformAll(pc, positions);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterImageOnly(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho,
        const BinPolicy_t& policy, const Hash_t& hash) const {
    if (!enabled_m) {
        return;
    }
    applyMirrorTransformSubset(pc, positions, policy, hash);
    scatterScaledDtSubset(pc, positions, rho, policy, hash);
    restoreMirrorTransformSubset(pc, positions, policy, hash);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterPrimaryAndImage(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho) const {
    scatterPrimaryOnly(pc, positions, rho);
    scatterImageOnly(pc, positions, rho);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterPrimaryAndImage(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho,
        const BinPolicy_t& policy, const Hash_t& hash) const {
    scatterPrimaryOnly(pc, positions, rho, policy, hash);
    scatterImageOnly(pc, positions, rho, policy, hash);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::applyMirrorTransformAll(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions) const {
    auto rView          = positions.getView();
    const size_t nLoc   = pc->getLocalNum();
    const double planeZ = zPlane_m;
    Kokkos::parallel_for(
            "ImageChargeScatterController::applyMirrorTransformAll", nLoc,
            KOKKOS_LAMBDA(const size_t i) { rView(i)[2] = 2.0 * planeZ - rView(i)[2]; });
    flipChargeSignAll(pc);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::restoreMirrorTransformAll(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions) const {
    applyMirrorTransformAll(pc, positions);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::applyMirrorTransformSubset(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, const BinPolicy_t& policy,
        const Hash_t& hash) const {
    auto rView          = positions.getView();
    const double planeZ = zPlane_m;
    Kokkos::parallel_for(
            "ImageChargeScatterController::applyMirrorTransformSubset", policy,
            KOKKOS_LAMBDA(const size_t i) {
                const size_t idx = hash(i);
                rView(idx)[2]    = 2.0 * planeZ - rView(idx)[2];
            });
    flipChargeSignSubset(pc, policy, hash);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::restoreMirrorTransformSubset(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, const BinPolicy_t& policy,
        const Hash_t& hash) const {
    applyMirrorTransformSubset(pc, positions, policy, hash);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::flipChargeSignAll(
        std::shared_ptr<ParticleCtr_t> pc) const {
    auto qView        = pc->getQView();
    const size_t nLoc = pc->getLocalNum();
    if (pc->getQMStorageMode() == ParticleCtr_t::QMStorageMode::Attributes) {
        Kokkos::parallel_for(
                "ImageChargeScatterController::flipChargeSignAllAttributes", nLoc,
                KOKKOS_LAMBDA(const size_t i) { qView(i) = -qView(i); });
    } else {
        Kokkos::parallel_for(
                "ImageChargeScatterController::flipChargeSignAllSingle", size_t(1),
                KOKKOS_LAMBDA(const size_t) { qView(0) = -qView(0); });
    }
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::flipChargeSignSubset(
        std::shared_ptr<ParticleCtr_t> pc, const BinPolicy_t& policy, const Hash_t& hash) const {
    auto qView = pc->getQView();
    if (pc->getQMStorageMode() == ParticleCtr_t::QMStorageMode::Attributes) {
        Kokkos::parallel_for(
                "ImageChargeScatterController::flipChargeSignSubsetAttributes", policy,
                KOKKOS_LAMBDA(const size_t i) {
                    const size_t idx = hash(i);
                    qView(idx)       = -qView(idx);
                });
    } else {
        Kokkos::parallel_for(
                "ImageChargeScatterController::flipChargeSignSubsetSingle", size_t(1),
                KOKKOS_LAMBDA(const size_t) { qView(0) = -qView(0); });
    }
}
