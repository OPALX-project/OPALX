template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterScaledDtAll(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho) const {
    Inform m("ImageChargeScatterController::scatterScaledDtAll");
    m << level5 << "all-local scatter entry: localP=" << pc->getLocalNum() << endl;
    pc->scaleDtByCharge();
    m << level5 << "all-local scatter: dt scaled by charge." << endl;
    ippl::ParticleAttrib<T>* dtAttrib = &pc->dt;
    scatter(*dtAttrib, rho, positions);
    m << level5 << "all-local scatter: scatter+halo done." << endl;
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
