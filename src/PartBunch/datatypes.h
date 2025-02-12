// some typedefs
// /todo this needs to go 

template <unsigned Dim>
using Mesh_t = ippl::UniformCartesian<double, Dim>;

template <typename T, unsigned Dim>
using PLayout_t = typename ippl::ParticleSpatialLayout<T, Dim, Mesh_t<Dim>>;

template <unsigned Dim>
using Centering_t = typename Mesh_t<Dim>::DefaultCentering;

template <unsigned Dim>
using FieldLayout_t = ippl::FieldLayout<Dim>;

template <typename T, unsigned Dim>
using Vector_t = ippl::Vector<T, Dim>;

const double pi = Kokkos::numbers::pi_v<double>;

constexpr unsigned Dim = 3;
using T                = double;
