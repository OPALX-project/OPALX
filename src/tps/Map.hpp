#ifndef TPS_CONTEXT_MAP_HPP
#define TPS_CONTEXT_MAP_HPP
#include <span>
#include <vector>
#include "tps/Algebra/Tps.hpp"
#include "tps/Algebra/Vps.hpp"
#include "tps/MapContext.hpp"
namespace tps {
    namespace detail {
        // Bridge the legacy scalar precision setting; restore it even on exceptions.
        // Neither this bridge nor the shared polynomial cache is thread-safe.
        template <class T>
        class ScopedMapOrder {
        public:
            explicit ScopedMapOrder(int order) : previous_(Tps<T>::getGlobalTruncOrder()) {
                Tps<T>::setGlobalTruncOrder(order);
            }
            ~ScopedMapOrder() { Tps<T>::setGlobalTruncOrder(previous_); }
            ScopedMapOrder(const ScopedMapOrder&)            = delete;
            ScopedMapOrder& operator=(const ScopedMapOrder&) = delete;

        private:
            int previous_;
        };
    }  // namespace detail
    /** @brief Policy for composing finite jets versus merely stored polynomials. */
    enum class CompositionPolicy {
        /** Preserve the expansion origin; reject a nonzero inner constant. */
        preserveOrigin,
        /** Evaluate stored polynomial composition; no inherited jet-accuracy claim. */
        storedPolynomial
    };
    /** @brief Immutable context-bearing polynomial map in local chart coordinates.
     * @tparam T Coefficient/evaluation scalar type, usually double.
     * @details Physical references, scaling and frames are metadata: evaluate() takes
     * local coordinates and returns local coordinates. No conversion is implicit.
     * Construction checks component dimensions, degree and available precision.
     * Global scalar truncation is temporarily overridden/restored by operations;
     * external synchronization is still required for concurrent map algebra.
     * @test tps.map_context checks composition, chart joins and precision restoration.
     * @snippet map_context.cpp context map composition
     */
    template <class T = double>
    class Map {
    public:
        Map(MapContext context, VpsMap<T> polynomial)
            : context_(std::move(context)), polynomial_(std::move(polynomial)) {
            if (polynomial_.getDimension()
                != static_cast<int>(context_.output().coordinates.size()))
                throw std::invalid_argument("Map: output dimension mismatch");
            detail::ScopedMapOrder<T> precision(Tps<T>::EXACT);
            for (int i = 0; i < polynomial_.getDimension(); ++i) {
                const auto& component = polynomial_[i];
                if ((component.getVariables() != 0
                     && component.getVariables()
                                != static_cast<int>(context_.input().coordinates.size()))
                    || component.getMaxOrder() > context_.order()
                    || component.getTruncOrder() < context_.order())
                    throw std::invalid_argument(
                            "Map: component dimension, stored degree or precision mismatch");
                // Snapshot coefficients: even a pre-existing source reference cannot
                // bypass copy-on-write and mutate this immutable map later.
                (void)polynomial_[i][0];
            }
        }
        Map(MapContext context, std::span<const Tps<T>> components)
            : Map(context, fromComponents(context, components)) {}
        [[nodiscard]] const MapContext& context() const noexcept { return context_; }
        [[nodiscard]] const VpsMap<T>& polynomial() const noexcept { return polynomial_; }
        /** @brief Evaluate local input coordinates; output is in the output chart. */
        [[nodiscard]] std::vector<T> evaluate(std::span<const T> point) const {
            auto input = asVector(point);
            std::vector<T> output(static_cast<std::size_t>(polynomial_.getDimension()));
            for (int i = 0; i < polynomial_.getDimension(); ++i) {
                const auto& component = polynomial_[i];
                output[i] =
                        component.getVariables() == 0 ? component[0] : component.evaluate(input);
            }
            return output;
        }
        /** @brief Jacobian in local input/output coordinates, evaluated at point. */
        [[nodiscard]] Matrix<T> jacobian(std::span<const T> point) const {
            auto input = asVector(point);
            Matrix<T> result(polynomial_.getDimension(), input.size());
            detail::ScopedMapOrder<T> precision(context_.order());
            for (int i = 0; i < result.nrows(); ++i)
                for (int j = 0; j < result.ncols(); ++j) {
                    const auto& component = polynomial_[i];
                    result(i, j) = (component.getVariables() == 0 || component.getMaxOrder() == 0)
                                           ? T(0)
                                           : component.derivative(j).evaluate(input);
                }
            return result;
        }
        /** @brief Project to a lower/equal order, retaining charts and provenance. */
        [[nodiscard]] Map truncated(int order) const {
            if (order < 0 || order > context_.order())
                throw std::invalid_argument("Map: invalid projection order");
            detail::ScopedMapOrder<T> precision(order);
            VpsMap<T> result(
                    polynomial_.getDimension(),
                    static_cast<int>(context_.input().coordinates.size()));
            for (int i = 0; i < result.getDimension(); ++i)
                result[i] = polynomial_[i].getVariables() == 0 ? polynomial_[i]
                                                               : polynomial_[i].filter(0, order);
            return Map(
                    MapContext(context_.input(), context_.output(), order, context_.segments()),
                    result);
        }

    private:
        static VpsMap<T> fromComponents(
                const MapContext& context, std::span<const Tps<T>> components) {
            if (components.size() != context.output().coordinates.size())
                throw std::invalid_argument("Map: component count mismatch");
            VpsMap<T> result(
                    static_cast<int>(components.size()),
                    static_cast<int>(context.input().coordinates.size()));
            for (std::size_t i = 0; i < components.size(); ++i)
                result[static_cast<int>(i)] = components[i];
            return result;
        }
        Vector<T> asVector(std::span<const T> point) const {
            if (point.size() != context_.input().coordinates.size())
                throw std::invalid_argument("Map: point dimension mismatch");
            Vector<T> result(static_cast<int>(point.size()));
            for (std::size_t i = 0; i < point.size(); ++i)
                result[static_cast<int>(i)] = point[i];
            return result;
        }
        MapContext context_;
        VpsMap<T> polynomial_;
    };
    /** @brief Identity in a single chart, at a positive finite order. */
    template <class T = double>
    [[nodiscard]] Map<T> identityMap(const PhaseSpaceChart& chart, int order) {
        MapContext context(chart, chart, order);
        if (order < 1) throw std::invalid_argument("identityMap: order must be positive");
        VpsMap<T> identity(
                static_cast<int>(chart.coordinates.size()),
                static_cast<int>(chart.coordinates.size()));
        for (int i = 0; i < identity.getDimension(); ++i)
            identity[i] = Tps<T>::makeVariable(identity.getDimension(), i);
        return Map<T>(context, identity);
    }
    /** @brief Compute the order-p projection of outer composed with inner; inner acts first.
     * @details Adjacent charts must match exactly; order cannot exceed either map's
     * available order. Repeated segment IDs are rejected to avoid overlap double
     * counting. By default a nonzero inner constant is rejected for jet validity.
     */
    template <class T>
    [[nodiscard]] Map<T> compose(
            const Map<T>& outer, const Map<T>& inner, int order,
            CompositionPolicy policy = CompositionPolicy::preserveOrigin) {
        if (inner.context().output() != outer.context().input())
            throw std::invalid_argument("compose: boundary charts differ");
        if (order < 0 || order > inner.context().order() || order > outer.context().order())
            throw std::invalid_argument("compose: invalid requested order");
        auto segments = inner.context().segments();
        segments.insert(
                segments.end(), outer.context().segments().begin(),
                outer.context().segments().end());
        MapContext context(
                inner.context().input(), outer.context().output(), order, std::move(segments));
        if (policy == CompositionPolicy::preserveOrigin)
            for (int i = 0; i < inner.polynomial().getDimension(); ++i)
                if (inner.polynomial()[i][0] != T(0))
                    throw std::invalid_argument("compose: inner map does not preserve origin");
        detail::ScopedMapOrder<T> precision(order);
        // Build the input-variable cache before scalar substitution takes references
        // into its traversal table. Products must not rebuild that table mid-walk.
        TpsData::getTpsData(order, static_cast<int>(context.input().coordinates.size()));
        VpsMap<T> result(
                outer.polynomial().getDimension(),
                static_cast<int>(context.input().coordinates.size()));
        for (int i = 0; i < result.getDimension(); ++i) {
            const auto& component = outer.polynomial()[i];
            result[i]             = component.getVariables() == 0 ? component
                                                                  : component.substitute(inner.polynomial());
            if (result[i].getVariables() != 0) result[i] = result[i].filter(0, order);
        }
        return Map<T>(context, result);
    }
    /** @brief Compose at the minimum of the two available context orders. */
    template <class T>
    [[nodiscard]] Map<T> compose(
            const Map<T>& outer, const Map<T>& inner,
            CompositionPolicy policy = CompositionPolicy::preserveOrigin) {
        return compose(
                outer, inner, std::min(outer.context().order(), inner.context().order()), policy);
    }
    /** @brief Compose nonempty maps in tracking order: maps[0] acts first.
     * @details Empty sequences are rejected; use identityMap for an explicit identity.
     */
    template <class T>
    [[nodiscard]] Map<T> composeSequence(
            std::span<const Map<T>> maps, int order,
            CompositionPolicy policy = CompositionPolicy::preserveOrigin) {
        if (maps.empty()) throw std::invalid_argument("composeSequence: empty sequence");
        auto result = maps.front().truncated(order);
        for (std::size_t i = 1; i < maps.size(); ++i)
            result = compose(maps[i], result, order, policy);
        return result;
    }
}  // namespace tps
#endif
