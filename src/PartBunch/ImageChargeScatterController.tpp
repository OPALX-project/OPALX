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

    m << level5 << "all-local scatter: host-staged accumulateHalo start." << endl;
    accumulateScalarHaloHostStaged(rho);
    m << level5 << "all-local scatter: host-staged accumulateHalo done." << endl;

    pc->unscaleDtByCharge();
    m << level5 << "all-local scatter: dt restored." << endl;
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::accumulateScalarHaloHostStaged(RhoField_t& rho) const {
    static_assert(Dim == 3, "Host-staged scalar halo accumulation currently supports Dim == 3.");

    Inform m("ImageChargeScatterController::accumulateScalarHaloHostStaged");

    auto rhoView = rho.getView();
    auto host    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rhoView);

    auto& layout          = rho.getLayout();
    const auto& neighbors = layout.getNeighbors();
    const auto& sendRange = layout.getNeighborsSendRange();
    const auto& recvRange = layout.getNeighborsRecvRange();
    MPI_Comm comm         = ippl::Comm->getCommunicator();

    const auto packRange = [&](const typename RhoField_t::Layout_t::bound_type& range) {
        std::vector<T> buffer(static_cast<size_t>(range.size()));
        size_t n = 0;
        for (long k = range.lo[2]; k < range.hi[2]; ++k) {
            for (long j = range.lo[1]; j < range.hi[1]; ++j) {
                for (long i = range.lo[0]; i < range.hi[0]; ++i) {
                    buffer[n++] = host(i, j, k);
                }
            }
        }
        return buffer;
    };

    const auto addRange =
            [&](const typename RhoField_t::Layout_t::bound_type& range,
                const std::vector<T>& buffer) {
                size_t n = 0;
                for (long k = range.lo[2]; k < range.hi[2]; ++k) {
                    for (long j = range.lo[1]; j < range.hi[1]; ++j) {
                        for (long i = range.lo[0]; i < range.hi[0]; ++i) {
                            host(i, j, k) += buffer[n++];
                        }
                    }
                }
            };

    size_t totalRequests = 0;
    for (const auto& componentNeighbors : neighbors) {
        totalRequests += componentNeighbors.size();
    }
    m << level5 << "host-staged halo: neighbor exchanges=" << totalRequests << endl;

    std::vector<std::vector<T>> sendBuffers;
    std::vector<MPI_Request> requests;
    sendBuffers.reserve(totalRequests);
    requests.reserve(totalRequests);

    constexpr size_t cubeCount = ippl::detail::countHypercubes(Dim) - 1;
    for (size_t index = 0; index < cubeCount; ++index) {
        const int tag = ippl::mpi::tag::HALO + static_cast<int>(index);
        const auto& componentNeighbors = neighbors[index];
        for (size_t i = 0; i < componentNeighbors.size(); ++i) {
            const int targetRank = componentNeighbors[i];
            sendBuffers.push_back(packRange(recvRange[index][i]));
            MPI_Request request;
            const auto& buffer = sendBuffers.back();
            MPI_Isend(
                    const_cast<T*>(buffer.data()),
                    static_cast<int>(buffer.size() * sizeof(T)),
                    MPI_BYTE, targetRank, tag, comm, &request);
            requests.push_back(request);
        }
    }

    for (size_t index = 0; index < cubeCount; ++index) {
        const int tag = ippl::mpi::tag::HALO
                        + static_cast<int>(RhoField_t::Layout_t::getMatchingIndex(index));
        const auto& componentNeighbors = neighbors[index];
        for (size_t i = 0; i < componentNeighbors.size(); ++i) {
            const int sourceRank = componentNeighbors[i];
            const auto& range    = sendRange[index][i];
            std::vector<T> recvBuffer(static_cast<size_t>(range.size()));
            MPI_Status status;
            MPI_Recv(
                    recvBuffer.data(),
                    static_cast<int>(recvBuffer.size() * sizeof(T)),
                    MPI_BYTE, sourceRank, tag, comm, &status);
            addRange(range, recvBuffer);
        }
    }

    if (!requests.empty()) {
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    }

    Kokkos::deep_copy(rhoView, host);
}

template <typename T, unsigned Dim>
void ImageChargeScatterController<T, Dim>::scatterScaledDtSubset(
        std::shared_ptr<ParticleCtr_t> pc, PositionAttr_t& positions, RhoField_t& rho,
        const BinPolicy_t& policy, const Hash_t& hash) const {
    Inform m("ImageChargeScatterController::scatterScaledDtSubset");
    const size_t nLocal = pc->getLocalNum();
    m << level5 << "subset scatter entry: localP=" << pc->getLocalNum() << ", policy=["
      << policy.begin() << "," << policy.end() << "), hashExtent=" << hash.extent(0) << endl;
    pc->scaleDtByCharge();
    m << level5 << "subset scatter: dt scaled by charge." << endl;

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

    int invalidHash = 0;
    Kokkos::parallel_reduce(
            "ImageChargeScatterController::validateSubsetHash", policy,
            KOKKOS_LAMBDA(const size_t i, int& update) {
                const size_t idx = hash(i);
                if (idx >= nLocal) {
                    update += 1;
                }
            },
            invalidHash);
    m << level5 << "subset scatter: invalid hash entries=" << invalidHash << endl;

    int outOfBounds = 0;
    Kokkos::parallel_reduce(
            "ImageChargeScatterController::validateSubsetScatterBounds", policy,
            KOKKOS_LAMBDA(const size_t i, int& update) {
                const size_t idx = hash(i);
                if (idx >= nLocal) {
                    return;
                }
                const auto l = (rView(idx) - origin) * invdx + 0.5;
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
    m << level5 << "subset scatter: out-of-bounds CIC particles=" << outOfBounds << endl;

    m << level5 << "subset scatter: CIC kernel start." << endl;
    Kokkos::parallel_for(
            "ImageChargeScatterController::scatterScaledDtSubsetCIC", policy,
            KOKKOS_LAMBDA(const size_t i) {
                const size_t idx = hash(i);
                if (idx >= nLocal) {
                    return;
                }

                const auto l = (rView(idx) - origin) * invdx + 0.5;
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
                            dtView(idx));
                }
            });
    Kokkos::fence();
    m << level5 << "subset scatter: CIC kernel done." << endl;

    m << level5 << "subset scatter: host-staged accumulateHalo start." << endl;
    accumulateScalarHaloHostStaged(rho);
    m << level5 << "subset scatter: host-staged accumulateHalo done." << endl;

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
