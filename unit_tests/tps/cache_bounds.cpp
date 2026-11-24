#include <Kokkos_Core.hpp>
#include "tps/Algebra/Tps.hpp"
#include <climits>
#include <iostream>
#include <stdexcept>

using namespace tps;
using TPS = Tps<double>;
void require(bool ok, const char* message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void range(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const CLRangeError&) {
        rejected = true;
    }
    require(rejected, "expected CLRangeError");
}
template <class F> void size(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const SizeError&) {
        rejected = true;
    }
    require(rejected, "expected SizeError");
}
int main() {
    Kokkos::ScopeGuard kokkos;

    try {
        // Invalid dimensions must be rejected before cache indexing or allocation.
        for (int n : {-1, 100, INT_MAX}) {
            range([&] { TpsData::getTpsData(0, n); });
            range([&] { TPS x(0, n); });
            range([&] { TpsMonomial m(n); });
        }
        require(TpsData::getTpsData(0, 0) == nullptr, "scalar cache sentinel");
        for (int order : {-1, INT_MAX, INT_MAX - 1})
            range([&] { TpsData::getTpsData(order, 2); });
        range([] { TPS x(-1, 2); });
        // Highest supported variable count, including cache growth.
        TpsMonomial zero(99);
        require(zero.getIndex() == 0, "upper variable boundary");
        auto last = TPS::makeVariable(99, 98);
        require(last.getCoefficient(99) == 1, "last valid variable");
        range([] { TPS::makeVariable(99, 99); });
        range([] { TPS::makeVariable(2, -1); });
        range([] { TPS::makeVariable(0, 0); });
        auto* data = TpsData::getTpsData(2, 2);
        const int terms = data->getSize(2);
        range([&] { data->getExponents(-1); });
        range([&] { data->getExponents(terms); });
        range([&] { data->getOrder(terms); });
        range([&] { data->getProductArray(terms); });
        range([&] { data->getProductArray(-1); });
        range([&] { data->getSize(-2); });
        range([&] { data->getSize(3); });
        // Preflight rejects combinatorial and dimension overflow without poisoning cache.
        range([] { TpsData::getTpsData(50, 99); });
        range([] { TpsData::getTpsData(INT_MAX - 2, 99); });
        require(TpsData::getTpsData(1, 99)->getSize(1) == 100, "cache survived rejected growth");
        require(TpsData::getTpsData(2, 2) == data && data->getSize(2) == terms, "cache unchanged");
        TpsMonomial m(2);
        range([&] { m[-1] = 1; });
        range([&] { m[2] = 1; });
        const TpsMonomial& cm = m;
        range([&] { (void)cm[2]; });
        range([] { TpsMonomial m(2, 2); });
        m[0] = -1;
        range([&] { m.getIndex(); });
        m[0] = INT_MAX;
        m[1] = 1;
        range([&] { m.getOrder(); });
        TpsMonomial a(2), b(2);
        a[0] = INT_MAX;
        b[0] = 1;
        range([&] { (void)(a * b); });
        size([] {
            TpsMonomial a(2), b(3);
            (void)(a * b);
        });
        TpsMonomial excessive(2);
        excessive[0] = 3;
        range([&] { data->indexMonomial(excessive); });
        size([&] {
            TpsMonomial mismatch(3);
            data->indexMonomial(mismatch);
        });
        auto x = TPS::makeVariable(2, 0);
        range([&] { x.derivative(2); });
        range([&] { x.integral(-1); });
        range([&] { x.multiplyVariable(2); });
        range([&] { x.getCoefficient(-1); });
        range([&] { x.getCoefficient(x.getSize()); });
        TPS scalar(2.0);
        range([&] { scalar.setCoefficient(1, 3.0); });
        range([&] { x.setCoefficient(-1, 3.0); });
        range([&] { x.setCoefficient(data->getSize(2), 3.0); });
        TpsMonomial negative(2);
        negative[1] = -1;
        range([&] { data->indexMonomial(negative); });
        // Successful arithmetic remains available after rejected input.
        auto square = x * x;
        TpsMonomial second(2);
        second[0] = 2;
        require(square[second] == 1, "arithmetic after rejected operations");
        std::cout << "Cache bounds checks passed\n";
    } catch (const Error& e) {
        std::cerr << e.where() << ": " << e.what() << '\n';
        return 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
