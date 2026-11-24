#include <Kokkos_Core.hpp>
#include "tps/Algebra/Array2D.h"
#include <climits>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

using namespace tps;
void require(bool ok, const char* message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const CLRangeError&) {
        rejected = true;
    }
    require(rejected, "expected CLRangeError");
}
int main() {
    Kokkos::ScopeGuard kokkos;

    try {
        Array1D<std::string> a(3, "value");
        const auto& ca = a;
        a[0] = "first";
        a[2] = "last";
        require(ca[0] == "first" && ca(2) == "last", "valid 1D boundaries");
        for (int i : {-1, 3, INT_MAX}) {
            rejects([&] { (void)ca[i]; });
            rejects([&] { (void)a[i]; });
            rejects([&] { a[i] = "bad"; });
            rejects([&] { (void)ca(i); });
            rejects([&] { a(i) = "bad"; });
        }
        auto copy = a;
        copy[0] = "copy";
        require(a[0] == "first", "1D copy independence");
        rejects([&] { a.resize(-1); });
        require(a.size() == 3 && a[2] == "last", "invalid resize changed array");
        a.resize(1);
        rejects([&] { (void)a[1]; });
        a.resize(4);
        a[3] = "grown";
        require(a[0] == "first" && ca[3] == "grown", "resize and checked access");
        rejects([] { Array1D<int> bad(-1); });
        rejects([] { Array1D<int> bad(-1, 0); });
        Array1D<int> empty;
        require(empty.begin() == empty.end(), "empty 1D range");
        rejects([&] { (void)empty[0]; });
        rejects([&] { (void)static_cast<const Array1D<int>&>(empty)[0]; });
        empty.resize(0);
        empty.resize(2);
        empty[1] = 9;
        require(empty[1] == 9, "empty resize");
        Array2D<int> m(2, 3, 0);
        const auto& cm = m;
        static_assert(std::is_same<decltype(cm[0][0]), const int&>::value,
                      "const row view must be immutable");
        m[0][0] = 1;
        m[1][2] = 6;
        require(cm[0][0] == 1 && cm(1, 2) == 6, "valid 2D boundaries");
        for (int r : {-1, 2, INT_MAX}) {
            rejects([&] { (void)m[r]; });
            rejects([&] { (void)cm[r]; });
            rejects([&] { m(r, 0) = 99; });
            rejects([&] { (void)cm(r, 0); });
            rejects([&] { m.row_begin(r); });
            rejects([&] { m.row_end(r); });
            rejects([&] { cm.row_begin(r); });
            rejects([&] { cm.row_end(r); });
            rejects([&] { m.swapRows(0, r); });
        }
        for (int c : {-1, 3, INT_MAX}) {
            rejects([&] { m[0][c] = 99; });
            rejects([&] { (void)cm[0][c]; });
            rejects([&] { m(0, c) = 99; });
            rejects([&] { (void)cm(0, c); });
            rejects([&] { m.col_begin(c); });
            rejects([&] { m.col_end(c); });
            rejects([&] { cm.col_begin(c); });
            rejects([&] { cm.col_end(c); });
            rejects([&] { m.swapColumns(0, c); });
        }
        require(m[0][0] == 1 && m[1][2] == 6, "rejected writes changed 2D array");
        Array1D<int> target(1, 42);
        rejects([&] { m.getRow(target, -1); });
        rejects([&] { m.getColumn(target, 3); });
        require(target.size() == 1 && target[0] == 42, "invalid extraction changed target");
        Array1D<int> row(3, 7), column(2, 8);
        rejects([&] { m.putRow(row, -1); });
        rejects([&] { m.putColumn(column, 3); });
        m.putRow(row, 0);
        m.putColumn(column, 2);
        m.getRow(target, 0);
        require(target.size() == 3 && target[2] == 8, "row transfer");
        m.getColumn(target, 2);
        require(target.size() == 2 && target[1] == 8, "column transfer");
        m.swapRows(0, 1);
        m.swapColumns(0, 2);
        require(m[0][0] == 8 && m[1][2] == 7, "valid swaps");
        Array2D<int> matrixCopy = m;
        matrixCopy[0][0] = 20;
        require(cm[0][0] == 8, "2D copy independence");
        for (auto shape : {std::pair<int, int>{-1, 2}, {2, -1}, {INT_MAX, 2}, {2, INT_MAX}}) {
            rejects([&] { Array2D<int> bad(shape.first, shape.second); });
            rejects([&] { Array2D<int> bad(shape.first, shape.second, 0); });
        }
        Array2D<int> defaultEmpty;
        const auto& cdefaultEmpty = defaultEmpty;
        require(defaultEmpty.begin() == defaultEmpty.end(), "default empty 2D range");
        rejects([&] { (void)defaultEmpty[0]; });
        rejects([&] { (void)cdefaultEmpty[0]; });
        Array2D<int> noColumns(2, 0);
        require(noColumns.row_begin(1) == noColumns.row_end(1), "zero-column row range");
        rejects([&] { noColumns[0][0] = 1; });
        rejects([&] { (void)static_cast<const Array2D<int>&>(noColumns)[0][0]; });
        Array2D<int> noRows(0, 3);
        require(noRows.col_begin(2) == noRows.col_end(2), "zero-row column range");
        require(static_cast<const Array2D<int>&>(noRows).col_begin(2) ==
                    static_cast<const Array2D<int>&>(noRows).col_end(2),
                "const empty column");
        rejects([&] { (void)noRows[0]; });
        std::cout << "Array indexing regressions passed\n";
    } catch (const Error& e) {
        std::cerr << e.where() << ": " << e.what() << '\n';
        return 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
