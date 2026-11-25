#include <Kokkos_Core.hpp>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "Utilities/CLRangeError.h"
#include "Utilities/OpalException.h"
#include "tps/Algebra/Array2D.h"

using namespace tps;
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
template <class Iterator>
void check(Iterator begin, Iterator end, int rows, int column) {
    require(end - begin == rows && begin - end == -rows, "iterator difference");
    require(std::distance(begin, end) == rows, "std::distance");
    const Iterator fixed = begin;
    require(fixed + rows == end && end - rows == begin, "logical end arithmetic");
    auto assigned = end;
    assigned      = begin;
    assigned += rows;
    require(assigned == end, "copy assignment/add");
    assigned -= rows;
    require(assigned == begin, "subtract assign");
    auto advanced = begin;
    std::advance(advanced, rows);
    require(advanced == end, "std::advance");
    int count = 0;
    for (auto it = begin; it != end; ++it, ++count) {
        require(*it == count * 10 + column, "forward traversal");
        require(begin[count] == *it, "iterator subscript");
    }
    require(count == rows, "traversal count");
    if (rows) {
        auto it  = begin;
        auto old = it++;
        require(old == begin, "postincrement old value");
        require(it == begin + 1, "postincrement advance");
        old = it--;
        require(old == begin + 1 && it == begin, "postdecrement");
        it = end;
        --it;
        require(*it == (rows - 1) * 10 + column, "decrement end");
        require((end - 1)[0] == *it, "end-relative access");
    }
    std::vector<int> copied(rows);
    std::copy(begin, end, copied.begin());
    require(std::accumulate(begin, end, 0) == 5 * rows * (rows - 1) + rows * column, "accumulate");
    auto reverse          = std::make_reverse_iterator(end);
    const auto reverseEnd = std::make_reverse_iterator(begin);
    count                 = rows;
    for (; reverse != reverseEnd; ++reverse)
        require(*reverse == (--count) * 10 + column, "reverse traversal");
    require(count == 0, "reverse count");
}
template <class F>
void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const ::CLRangeError&) {
        rejected = true;
    }
    require(rejected, "expected iterator exception");
}
int main() {
    Kokkos::ScopeGuard kokkos;

    try {
        for (int rows : {0, 1, 2, 5})
            for (int cols : {1, 2, 4}) {
                Array2D<int> matrix(rows, cols);
                for (int r = 0; r < rows; ++r)
                    for (int c = 0; c < cols; ++c)
                        matrix(r, c) = 10 * r + c;
                const auto& cm = matrix;
                for (int c = 0; c < cols; ++c) {
                    check(matrix.col_begin(c), matrix.col_end(c), rows, c);
                    check(cm.col_begin(c), cm.col_end(c), rows, c);
                    if (rows) {
                        auto tail = matrix.col_end(c);
                        --tail;
                        *tail = 123;
                        require(matrix(rows - 1, c) == 123, "mutable tail write");
                        matrix(rows - 1, c) = 10 * (rows - 1) + c;
                    }
                    Array1D<int> copied;
                    matrix.getColumn(copied, c);
                    require(copied.size() == rows, "column copy size");
                    for (int r = 0; r < rows; ++r)
                        require(copied[r] == 10 * r + c, "column copy value");
                }
            }
        int values[] = {0, 1, 2, 3, 4, 5};
        // A reverse-stride end is represented logically, never as values-1.
        SliceIterator<int> first(values + 5, -2), last(values + 5, -2, 3);
        require(last - first == 3 && first + 3 == last, "negative stride range");
        require(*first == 5 && first[2] == 1 && *(last - 1) == 1, "negative stride values");
        ConstSliceIterator<int> cf(values + 5, -2), cl(values + 5, -2, 3);
        require(cl - cf == 3 && cf + 3 == cl && *(cl - 1) == 1, "const negative stride range");
        rejects([&] {
            SliceIterator<int> invalid(values, 0);
        });
        rejects([&] {
            ConstSliceIterator<int> invalid(values, 0);
        });
        rejects([&] {
            (void)(first - SliceIterator<int>(values, 1));
        });
        rejects([&] {
            (void)(cf - ConstSliceIterator<int>(values, 1));
        });
        SliceIterator<int> empty(nullptr, 1), emptyEnd(nullptr, 1, 0);
        require(empty == emptyEnd && emptyEnd - empty == 0, "null-origin empty range");
        ConstSliceIterator<int> ce(nullptr, 1), cend(nullptr, 1, 0);
        require(ce == cend && std::distance(ce, cend) == 0, "const null empty range");
        std::cout << "Column iterator regressions passed\n";
    } catch (const ::OpalException& e) {
        std::cerr << e.where() << ": " << e.what() << '\n';
        return 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
