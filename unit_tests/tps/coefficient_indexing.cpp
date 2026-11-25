#include <Kokkos_Core.hpp>
#include <climits>
#include <iostream>
#include <stdexcept>
#include "Utilities/CLRangeError.h"
#include "Utilities/LogicalError.h"
#include "Utilities/OpalException.h"
#include "tps/Algebra/Tps.hpp"

using namespace tps;
using TPS = Tps<double>;
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
template <class Exception, class F>
void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const Exception&) {
        rejected = true;
    }
    require(rejected, "expected indexing exception");
}
int main() {
    Kokkos::ScopeGuard kokkos;

    try {
        for (int variables : {0, 1, 2, 6})
            for (int order : {0, 1, 3}) {
                TPS a(order, variables);
                a[0]           = 4;
                const int last = a.getSize() - 1;
                a[last]        = 7;
                const TPS& ca  = a;
                require(ca[0] == a[0] && ca[last] == 7, "integer boundary access");
                for (int invalid : {-1, a.getSize(), INT_MAX}) {
                    rejects<::CLRangeError>([&] {
                        (void)ca[invalid];
                    });
                    rejects<::CLRangeError>([&] {
                        (void)a[invalid];
                    });
                    rejects<::CLRangeError>([&] {
                        a[invalid] = 99;
                    });
                }
                TpsMonomial zero(variables);
                require(ca[zero] == a[0], "constant monomial access");
                a[zero] = 5;
                require(ca[zero] == 5, "mutable monomial access");
                TpsMonomial mismatch(variables + 1);
                rejects<::LogicalError>([&] {
                    (void)ca[mismatch];
                });
                rejects<::LogicalError>([&] {
                    a[mismatch] = 99;
                });
                if (variables) {
                    TpsMonomial excessive(variables);
                    excessive[0] = order + 1;
                    rejects<::CLRangeError>([&] {
                        (void)ca[excessive];
                    });
                    rejects<::CLRangeError>([&] {
                        a[excessive] = 99;
                    });
                    rejects<::CLRangeError>([&] {
                        ca.getCoefficient(excessive);
                    });
                    TpsMonomial negative(variables);
                    negative[0] = -1;
                    rejects<::CLRangeError>([&] {
                        (void)ca[negative];
                    });
                    rejects<::CLRangeError>([&] {
                        a[negative] = 99;
                    });
                }
                require(a[0] == 5, "rejected write changed constant");
            }
        // Rejected high-order monomials must not grow the shared cache either.
        TPS fresh(1, 3);
        TpsMonomial tooHigh(3);
        tooHigh[0] = 6;
        rejects<::CLRangeError>([&] {
            (void)fresh[tooHigh];
        });
        rejects<::CLRangeError>([&] {
            static_cast<const TPS&>(fresh).getCoefficient(tooHigh);
        });
        auto* unchangedCache = TpsData::getTpsData(1, 3);
        rejects<::CLRangeError>([&] {
            unchangedCache->getSize(2);
        });
        // An expanded global cache must not authorize indexing past this object's size.
        TPS small(1, 2);
        small[1]      = 2;
        auto* cache   = TpsData::getTpsData(5, 2);
        const TPS& cs = small;
        TpsMonomial fifth(2);
        fifth[0] = 5;
        rejects<::CLRangeError>([&] {
            (void)cs[fifth];
        });
        rejects<::CLRangeError>([&] {
            small[fifth] = 9;
        });
        rejects<::CLRangeError>([&] {
            (void)small[cache->getSize(5) - 1];
        });
        require(small.getMaxOrder() == 1 && small[1] == 2, "indexing expanded or changed object");
        // Invalid mutable requests must leave both copies unchanged; valid writes detach.
        TPS copy = small;
        rejects<::CLRangeError>([&] {
            (void)small[-1];
        });
        rejects<::CLRangeError>([&] {
            small[fifth] = 9;
        });
        rejects<::LogicalError>([&] {
            TpsMonomial wrong(3);
            small[wrong] = 9;
        });
        require(cs[1] == 2 && static_cast<const TPS&>(copy)[1] == 2,
                "rejected indexing changed shared values");
        small[1] = 3;
        require(static_cast<const TPS&>(copy)[1] == 2, "integer copy independence");
        TpsMonomial x(2);
        x[0]     = 1;
        copy     = small;
        small[x] = 4;
        require(static_cast<const TPS&>(copy)[x] == 3, "monomial copy independence");
        // Explicit setter retains its documented ability to grow storage.
        small.setCoefficient(fifth, 6);
        require(small[fifth] == 6, "explicit setter growth");
        std::cout << "Coefficient indexing checks passed\n";
    } catch (const ::OpalException& e) {
        std::cerr << e.where() << ": " << e.what() << '\n';
        return 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
