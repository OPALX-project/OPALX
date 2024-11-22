# Add `AdaptBins` to your program:
- TODO


# What I changed in `scatter` to allow custom `Kokkos::parallel_for` iteration:

```c++
template <typename T, class... Properties>
template <typename Field, class PT, typename policy_type>
void ParticleAttrib<T, Properties...>::scatter(
    Field& f, const ParticleAttrib<Vector<PT, Field::dim>, Properties...>& pp,
    policy_type iteration_policy) const {
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
    Kokkos::parallel_for(
        "ParticleAttrib::scatter", iteration_policy,
        KOKKOS_CLASS_LAMBDA(const size_t idx) {
            // find nearest grid point
            vector_type l                        = (pp(idx) - origin) * invdx + 0.5;
            Vector<int, Field::dim> index        = l;
            Vector<PositionType, Field::dim> whi = l - index;
            Vector<PositionType, Field::dim> wlo = 1.0 - whi;

            Vector<size_t, Field::dim> args = index - lDom.first() + nghost;

            // scatter
            const value_type& val = dview_m(idx);
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
    attrib.scatter(f, pp, policy_type(0, attrib.getParticleCount())); // *(attrib.localNum_mp)
}

// Second implementation for custom range policy
template <typename Attrib1, typename Field, typename Attrib2, typename policy_type = Kokkos::RangePolicy<typename Field::execution_space>>
inline void scatter(const Attrib1& attrib, Field& f, const Attrib2& pp, policy_type iteration_policy) {
    attrib.scatter(f, pp, iteration_policy);
}
```
These do not change `scatter(...)` function calls, but allow the user to pass custom range policies for the scatter.

