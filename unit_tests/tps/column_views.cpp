#include <Kokkos_Core.hpp>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include "tps/Algebra/Array2D.h"

namespace {
    void require(bool condition, const char* message) {
        if (!condition) throw std::runtime_error(message);
    }
    void checkColumns(int rows, int columns) {
        tps::Array2D<double> matrix(rows, columns);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < columns; ++j)
                matrix(i, j) = 10 * i + j;
        const auto& constant = matrix;
        for (int j = 0; j < columns; ++j) {
            //! [column access]
            // j is a valid column; both views borrow matrix storage.
            auto column         = Kokkos::subview(matrix.view(), Kokkos::ALL(), j);
            const auto readOnly = Kokkos::subview(constant.view(), Kokkos::ALL(), j);
            static_assert(!std::is_assignable_v<decltype(readOnly(0)), double>);
            require(column.extent(0) == static_cast<std::size_t>(rows), "column extent");
            if (rows > 1)
                require(column.stride(0) == static_cast<std::size_t>(columns), "column stride");

            // deep_copy handles the strided column; raw pointer copying would not.
            Kokkos::View<double*, Kokkos::HostSpace> copy("column copy", rows);
            Kokkos::deep_copy(copy, readOnly);
            for (int i = 0; i < rows; ++i) {
                require(copy(i) == 10 * i + j, "column copy value");
                copy(i) += 100;
            }
            Kokkos::deep_copy(column, copy);
            //! [column access]
            for (int i = 0; i < rows; ++i)
                for (int c = 0; c < columns; ++c)
                    require(matrix(i, c) == 10 * i + c + (c == j ? 100 : 0),
                            "write changed another column");
            if (rows) {
                column(rows - 1) = -3;
                require(constant(rows - 1, j) == -3, "column view did not alias storage");
                require(copy(rows - 1) == 10 * (rows - 1) + j + 100,
                        "independent column copy changed");
            }
            // Restore this column before checking the next one.
            for (int i = 0; i < rows; ++i)
                column(i) = 10 * i + j;
        }
    }
}  // namespace

int main() {
    Kokkos::ScopeGuard runtime;
    try {
        for (int rows : {0, 1, 2, 4, 17})
            for (int columns : {1, 3, 5})
                checkColumns(rows, columns);
        tps::Array2D<double> empty(3, 0);
        require(empty.view().extent(0) == 3 && empty.view().extent(1) == 0, "zero-column shape");
        std::cout << "Kokkos column subview checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
