#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "AbsBeamline/EndFieldModel/Enge.h"

namespace enge_test {
    constexpr int maxOrder = endfieldmodel::EngeConfig::max_derivative;
    using Derivatives      = std::array<long double, maxOrder + 1>;

    constexpr std::int64_t factorial(int n) {
        std::int64_t result = 1;
        for (int i = 2; i <= n; ++i)
            result *= i;
        return result;
    }

    constexpr std::int64_t binomial(int n, int k) {
        return factorial(n) / (factorial(k) * factorial(n - k));
    }

    struct EndReference {
        Derivatives g{}, reciprocal{};
    };

    /** Independent derivatives from exp(h)' = h' exp(h) and g*(1/g) = 1.
     * Uses the product rule in long double, without production index tables.
     */
    inline EndReference referenceEnd(
            const std::vector<double>& coefficients, long double x, long double lambda) {
        Derivatives h{}, exponential{};
        for (int n = 0; n <= maxOrder; ++n) {
            for (int i = n; i < static_cast<int>(coefficients.size()); ++i) {
                h[n] += static_cast<long double>(coefficients[i]) * factorial(i) / factorial(i - n)
                        * std::pow(x, i - n) / std::pow(lambda, i);
            }
        }
        exponential[0] = std::exp(h[0]);
        for (int n = 1; n <= maxOrder; ++n)
            for (int k = 1; k <= n; ++k)
                exponential[n] += binomial(n - 1, k - 1) * h[k] * exponential[n - k];
        EndReference result;
        result.g = exponential;
        result.g[0] += 1;
        result.reciprocal[0] = 1 / result.g[0];
        for (int n = 1; n <= maxOrder; ++n) {
            for (int k = 1; k <= n; ++k)
                result.reciprocal[n] -= binomial(n, k) * result.g[k] * result.reciprocal[n - k];
            result.reciprocal[n] /= result.g[0];
        }
        return result;
    }

    /** Frozen pre-compaction product-rule expansion, kept on the host only.
     * Its order-12 table has 140152 rows, including duplicate monomials.
     */
    inline std::vector<endfieldmodel::EngeDerivativeTable> expandedTables() {
        std::vector<endfieldmodel::EngeDerivativeTable> tables{{{1, -1, 0}}};
        for (int n = 1; n <= maxOrder; ++n) {
            endfieldmodel::EngeDerivativeTable next;
            for (const auto& previous : tables.back()) {
                auto term = previous;
                term[0] *= term[1];
                --term[1];
                ++term[2];
                next.push_back(std::move(term));
                for (size_t k = 2; k < previous.size(); ++k) {
                    if (previous[k] == 0) continue;
                    term = previous;
                    if (k == previous.size() - 1) term.push_back(0);
                    term[0] *= term[k];
                    --term[k];
                    ++term[k + 1];
                    next.push_back(std::move(term));
                }
            }
            tables.push_back(std::move(next));
        }
        return tables;
    }

    /// Return the expanded value and sum of absolute terms (a cancellation scale).
    inline std::pair<long double, long double> evaluateTable(
            const endfieldmodel::EngeDerivativeTable& table, const Derivatives& g) {
        // Precompute powers to keep the legacy oracle cheap even in Debug builds.
        std::array<std::array<long double, maxOrder + 1>, maxOrder + 1> powers{};
        std::array<long double, maxOrder + 2> inversePowers{};
        for (int j = 0; j <= maxOrder; ++j)
            for (int k = 0; k <= maxOrder; ++k)
                powers[j][k] = std::pow(g[j], k);
        for (int k = 1; k <= maxOrder + 1; ++k)
            inversePowers[k] = std::pow(g[0], -k);
        long double value = 0, correction = 0, absoluteTerms = 0;
        for (const auto& row : table) {
            long double term = row[0] * inversePowers[-row[1]];
            for (size_t j = 2; j < row.size(); ++j)
                term *= powers[j - 1][row[j]];
            // long double is only 64 bits on some hosts. Compensate the sum so
            // accumulation of 140152 legacy terms does not dominate the oracle.
            const long double next = value + term;
            correction += std::abs(value) >= std::abs(term) ? (value - next) + term
                                                            : (term - next) + value;
            value = next;
            absoluteTerms += std::abs(term);
        }
        return {value + correction, absoluteTerms};
    }
}  // namespace enge_test
