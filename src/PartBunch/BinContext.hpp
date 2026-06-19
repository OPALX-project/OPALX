#ifndef OPALX_BIN_CONTEXT_HPP
#define OPALX_BIN_CONTEXT_HPP

#include <cstdint>

#include "PartBunch/FieldContainer.hpp"

template <typename T, unsigned Dim>
struct BinContext {
    long long binIndex                  = -1;
    std::uint64_t particleCountGlobal   = 0;
    double gamma                        = 1.0;
    Vector_t<double, Dim> pmean         = Vector_t<double, Dim>(0.0);
    Vector_t<double, Dim> meshSpacing   = Vector_t<double, Dim>(0.0);
    Vector_t<double, Dim> solveSpacing  = Vector_t<double, Dim>(0.0);
    bool imageCorrectionActive          = false;
    bool shiftedGreensCorrectionActive  = false;

    bool hasCorrection() const {
        return imageCorrectionActive || shiftedGreensCorrectionActive;
    }
};

#endif  // OPALX_BIN_CONTEXT_HPP
