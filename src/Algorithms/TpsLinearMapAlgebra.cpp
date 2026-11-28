#include "Utilities/SingularMatrixError.h"
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#include <algorithm>
#include <cmath>
#include "Algorithms/TpsLinearMapAlgebra.h"
#include "Utilities/OpalException.h"
#include "tps/LinearSolve.hpp"

namespace opalx::maps {
    namespace {
        tps::Matrix<double> dense(const matrix6x6_t& matrix) {
            tps::Matrix<double> result(6, 6, 0.0);
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) {
                    if (!std::isfinite(matrix(i, j)))
                        throw OpalException("TPS linear map algebra", "Nonfinite matrix entry.");
                    result(i, j) = matrix(i, j);
                }
            return result;
        }
    }  // namespace

    matrix6x6_t solveInputDifferences(
            const matrix6x6_t& input, const matrix6x6_t& output, double& condition) {
        const auto a = dense(input);
        auto rhs     = dense(output).transpose();
        try {
            const tps::LUFactorization<double> factor(a, 1e-18);
            condition = tps::infinityNorm(a) * tps::infinityNorm(factor.inverse());
            if (!std::isfinite(condition) || condition > 1e12)
                throw OpalException(
                        "LinearTransferMapBuilder::makeLinearTransferMap",
                        "The finite-difference input matrix is ill-conditioned (condition number "
                                + std::to_string(condition) + ").");
            factor.solveTranspose(rhs);
            matrix6x6_t result(0.0);
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) {
                    if (!std::isfinite(rhs(j, i)))
                        throw OpalException("TPS linear map algebra", "Nonfinite solve result.");
                    result(i, j) = rhs(j, i);
                }
            return result;
        } catch (const ::SingularMatrixError&) {
            throw OpalException(
                    "LinearTransferMapBuilder::makeLinearTransferMap",
                    "The finite-difference input matrix is singular.");
        } catch (const std::runtime_error& error) {
            throw OpalException("TPS linear map solve", error.what());
        }
    }

    LinearMapDiagnostics diagnoseLinearMap(const matrix6x6_t& map) {
        const auto matrix = dense(map);
        try {
            tps::Matrix<double> j(6, 6, 0.0);
            for (int pair = 0; pair < 3; ++pair) {
                j(2 * pair, 2 * pair + 1) = 1;
                j(2 * pair + 1, 2 * pair) = -1;
            }
            const auto residual = matrix.transpose() * (j * matrix) - j;
            double maximum      = 0;
            for (int i = 0; i < 6; ++i)
                for (int k = 0; k < 6; ++k) {
                    if (!std::isfinite(residual(i, k)))
                        throw OpalException(
                                "TPS linear map diagnostics", "Nonfinite canonical-form residual.");
                    maximum = std::max(maximum, std::abs(residual(i, k)));
                }
            const double determinantError = std::abs(tps::determinant(matrix) - 1);
            if (!std::isfinite(determinantError))
                throw OpalException(
                        "TPS linear map diagnostics", "Nonfinite determinant residual.");
            return {determinantError, maximum};
        } catch (const std::runtime_error& error) {
            throw OpalException("TPS linear map diagnostics", error.what());
        }
    }
}  // namespace opalx::maps
