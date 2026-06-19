#ifndef OPAL_SOLVE_REQUEST_H
#define OPAL_SOLVE_REQUEST_H

#include <optional>

#include "Manager/BaseManager.h"

template <typename T, unsigned Dim>
struct SolveRequest {
    std::optional<ippl::Vector<T, Dim>> greensShift;

    bool hasShiftedGreens() const { return greensShift.has_value(); }
};

#endif  // OPAL_SOLVE_REQUEST_H
