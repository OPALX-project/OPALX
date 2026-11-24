// Compile the kernel interface independently of host TPS/context/preparation headers.
#include <tps/DeviceMap.hpp>
#include <complex>
static_assert(tps::DeviceMapScalar<double> && tps::DeviceMapScalar<float>);
static_assert(!tps::DeviceMapScalar<std::complex<double>>);
static_assert(!tps::DeviceMapScalar<int>);
static_assert(std::is_trivially_copyable_v<tps::DeviceSpan<double>>);
