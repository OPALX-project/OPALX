# Add `AdaptBins` to your program:
- TODO


# What I changed in `scatter` to allow custom `Kokkos::parallel_for` iteration:

If necessary, iterate over `hash_array`. This allows to iterate over a sorted array without actually sorting them in memory, which is less efficient.
```c++
template <typename T, class... Properties>
template <typename Field, class PT, typename policy_type>
void ParticleAttrib<T, Properties...>::scatter(
    Field& f, const ParticleAttrib<Vector<PT, Field::dim>, Properties...>& pp,
    policy_type iteration_policy, hash_type hash_array) const {
    constexpr unsigned Dim = Field::dim;
    using PositionType     = typename Field::Mesh_t::value_type;

    static IpplTimings::TimerRef scatterTimer = IpplTimings::getTimer("scatter");
    IpplTimings::startTimer(scatterTimer);
    using view_type = typename Field::view_type;
    view_type view  = f.getView();

    using mesh_type       = typename Field::Mesh_t;
    const mesh_type& mesh = f.get_mesh();

    using vector_type = typename mesh_type::vector_type;
    using value_type  = typename ParticleAttrib<T, Properties...>::value_type;

    const vector_type& dx     = mesh.getMeshSpacing();
    const vector_type& origin = mesh.getOrigin();
    const vector_type invdx   = 1.0 / dx;

    const FieldLayout<Dim>& layout = f.getLayout();
    const NDIndex<Dim>& lDom       = layout.getLocalNDIndex();
    const int nghost               = f.getNghost();

    //using policy_type = Kokkos::RangePolicy<execution_space>;
    const bool useHashView = hash_array.extent(0) > 0;
    if (useHashView && (iteration_policy.end() > hash_array.extent(0))) {
        Inform m("scatter");
        m << "Hash array was passed to scatter, but size does not match iteration policy." << endl;
        ippl::Comm->abort();
    }
    Kokkos::parallel_for(
        "ParticleAttrib::scatter", iteration_policy,
        KOKKOS_CLASS_LAMBDA(const size_t idx) {
            // map index to possible hash_map
            size_t mapped_idx = useHashView ? hash_array(idx) : idx;

            // find nearest grid point
            vector_type l                        = (pp(mapped_idx) - origin) * invdx + 0.5;
            Vector<int, Field::dim> index        = l;
            Vector<PositionType, Field::dim> whi = l - index;
            Vector<PositionType, Field::dim> wlo = 1.0 - whi;

            Vector<size_t, Field::dim> args = index - lDom.first() + nghost;

            // scatter
            const value_type& val = dview_m(mapped_idx);
            detail::scatterToField(std::make_index_sequence<1 << Field::dim>{}, view, wlo, whi,
                                    args, val);
        });
    IpplTimings::stopTimer(scatterTimer);

    static IpplTimings::TimerRef accumulateHaloTimer = IpplTimings::getTimer("accumulateHalo");
    IpplTimings::startTimer(accumulateHaloTimer);
    f.accumulateHalo();
    IpplTimings::stopTimer(accumulateHaloTimer);
}

template <typename Attrib1, typename Field, typename Attrib2, typename policy_type = Kokkos::RangePolicy<typename Field::execution_space>>
inline void scatter(const Attrib1& attrib, Field& f, const Attrib2& pp) {
    attrib.scatter(f, pp, policy_type(0, attrib.getParticleCount())); 
}

// Second implementation for custom range policy, but without index array
template <typename Attrib1, typename Field, typename Attrib2, typename policy_type = Kokkos::RangePolicy<typename Field::execution_space>>
inline void scatter(const Attrib1& attrib, Field& f, const Attrib2& pp, policy_type iteration_policy, typename Attrib1::hash_type hash_array = {}) {
    attrib.scatter(f, pp, iteration_policy, hash_array);
}
```
These do not change old `scatter(...)` function calls, but allow the user to pass custom range policies for the scatter.


# Change `ParticleAttribBase` to allow overwrite unpack:

Next, this allows to sort all attributes of a particle bunch. First you need to adapt `ParticleAttribBase.h` as follows:
```c++
...
virtual void unpack(size_type, bool overwrite = false) = 0; // Add overwrite parameter (with default to make it compatible with the rest)
...
```
Next you change the implementation in `PartAttrib.hpp` as follows:
```c++
template <typename T, class... Properties>
void ParticleAttrib<T, Properties...>::unpack(size_type nrecvs, bool overwrite) {
    auto size          = dview_m.extent(0);
    size_type required = overwrite ? nrecvs : (*(this->localNum_mp) + nrecvs); // Change this (more memory efficient)!
    if (size < required) {
        int overalloc = Comm->getDefaultOverallocation();
        this->resize(required * overalloc);
    }

    size_type count   = overwrite ? 0 : *(this->localNum_mp); // Changed this!
    using policy_type = Kokkos::RangePolicy<execution_space>;
    Kokkos::parallel_for(
        "ParticleAttrib::unpack()", policy_type(0, nrecvs),
        KOKKOS_CLASS_LAMBDA(const size_t i) { dview_m(count + i) = buf_m(i); });
    Kokkos::fence();
}
```
Now when you want to apply an index array after "argsort", you can simply call the following:
```c++
bunch_m->template forAllAttributes([&]<typename Attribute>(Attribute*& attribute) {
    using memory_space    = typename Attribute::memory_space;

    // Ensure indices are in the correct memory space --> copies data ONLY when different memory spaces, so should be efficient
    auto indices_device = Kokkos::create_mirror_view_and_copy(memory_space{}, indices);

    attribute->pack(indices_device);
    attribute->unpack(localNumParticles, true);
});
```
Note that your index array needs to have e.g. the following type (execution space needs to be the same as the attribute, the create_mirror_view copies it if necessary):
```c++
using hash_type = ippl::detail::hash_type<Kokkos::DefaultExecutionSpace::memory_space>;
hash_type indices("indices", localNumParticles);
```
That way, all values in the corresponding view get put into the internal buffer (already exists in `ParticleAttrib.h`) according to the `hash`, which here is simply `indices`. Then, `unpack` puts all values from the buffer back into the main view, now in the new order. 

Make sure that `indices` has length `localNumParticles` and that `overwrite = true`. Otherwise, you will probably get segmentation faults.

Also: For some reason, `hash_type` is `Kokkos::View<int*>` instead of `Kokkos::View<size_type*>`. But I guess maybe that is enough...


# Finally change `gather` to allow `+=`

Change the following:
```c++
template <typename T, class... Properties>
template <typename Field, typename P2>
void ParticleAttrib<T, Properties...>::gather(
    Field& f, const ParticleAttrib<Vector<P2, Field::dim>, Properties...>& pp,
    const bool addToAttribute) {
    constexpr unsigned Dim = Field::dim;
    using PositionType     = typename Field::Mesh_t::value_type;

    static IpplTimings::TimerRef fillHaloTimer = IpplTimings::getTimer("fillHalo");
    IpplTimings::startTimer(fillHaloTimer);
    f.fillHalo();
    IpplTimings::stopTimer(fillHaloTimer);

    static IpplTimings::TimerRef gatherTimer = IpplTimings::getTimer("gather");
    IpplTimings::startTimer(gatherTimer);
    const typename Field::view_type view = f.getView();

    using mesh_type       = typename Field::Mesh_t;
    const mesh_type& mesh = f.get_mesh();

    using vector_type = typename mesh_type::vector_type;

    const vector_type& dx     = mesh.getMeshSpacing();
    const vector_type& origin = mesh.getOrigin();
    const vector_type invdx   = 1.0 / dx;

    const FieldLayout<Dim>& layout = f.getLayout();
    const NDIndex<Dim>& lDom       = layout.getLocalNDIndex();
    const int nghost               = f.getNghost();

    using policy_type = Kokkos::RangePolicy<execution_space>;
    Kokkos::parallel_for(
        "ParticleAttrib::gather", policy_type(0, *(this->localNum_mp)),
        KOKKOS_CLASS_LAMBDA(const size_t idx) {
            // find nearest grid point
            vector_type l                        = (pp(idx) - origin) * invdx + 0.5;
            Vector<int, Field::dim> index        = l;
            Vector<PositionType, Field::dim> whi = l - index;
            Vector<PositionType, Field::dim> wlo = 1.0 - whi;

            Vector<size_t, Field::dim> args = index - lDom.first() + nghost;

            // gather
            value_type gathered = detail::gatherFromField(std::make_index_sequence<1 << Field::dim>{},
                                                            view, wlo, whi, args);
            if (addToAttribute) dview_m(idx) += gathered;
            else                dview_m(idx)  = gathered;
        });
    IpplTimings::stopTimer(gatherTimer);
}

template <typename Attrib1, typename Field, typename Attrib2>
inline void gather(Attrib1& attrib, Field& f, const Attrib2& pp, 
                    const bool addToAttribute = false) {
    attrib.gather(f, pp, addToAttribute);
}
```
That way, the attribut can be either overwritten or added to the current one. Might be used to gather a partial field from the Lorentz transformation.
