/*
 *  Copyright (c) 2017, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ENDFIELDMODEL_ENGE_H_
#define ENDFIELDMODEL_ENGE_H_

#include <cstdint>
#include <iostream>
#include <vector>

#include <mutex>
#include "AbsBeamline/EndFieldModel/CompactVector.h"

#include "AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "Utilities/OpalException.h"

namespace endfieldmodel {

    /// Rows contain an integer coefficient followed by derivative powers.
    using EngeDerivativeTable = std::vector<std::vector<std::int64_t>>;

    /** Enge class is a symmetric Enge function for analytical field models
     *
     *  Enge function is
     *
     *  \f$f(x) = 1/(1+exp(h(x-x0)))+1/(1+exp(h(-x-x0)))-1\f$.
     *
     *  where h is a polynomial in x/lambda with polynomial coefficients a
     */

    /** Workspace for orders 0..max_derivative, including both endpoints.
     * Only the requested prefix is written/read during evaluation.
     */
    struct EngeWorkspace {
        static constexpr int max_derivative = 12;
        double exp_h0;
        double gNVec[max_derivative + 1];
        double hNVec[max_derivative + 1];
        double enge[max_derivative + 1];
    };
    /** GPU-ready configuration data */
    struct EngeConfig {
        static constexpr int max_derivative = EngeWorkspace::max_derivative;
        Kokkos::View<double*> a_m;
        double lambda_m = 0.0;
        double x0_m     = 0.0;
        /** Compacted derivatives of 1/g in terms of g and its derivatives.
         * Identical monomials are combined exactly before transfer to the device.
         * Order 12 needs 64-bit coefficients (the largest is 5388768000).
         * All supported coefficients are exactly representable as doubles.
         * Combining terms changes floating-point summation order, not the formula.
         */
        Kokkos::View<std::int64_t**> gIndices[max_derivative + 1];
        Kokkos::View<int**> hIndices[max_derivative + 1];
    };

    class Enge : public EndFieldModel {
    public:
        /** Default constructor */
        Enge() {}
        /** Builds Enge function with parameters a_0, a_1, ..., lambda and x0.
         *
         *  Note that this class is in the inner loop of tracking, so many function
         *  calls are _not_ checked for correct indexing. Call setMaximumDerivative
         *  before use.
         */
        Enge(const std::vector<double> a, double x0, double lambda);
        Enge(Kokkos::View<double*> a, double x0, double lambda);

        /** Destructor - no mallocs, so does nothing */
        ~Enge();

        /** Inheritable copy constructor - no mallocs, so does nothing */
        Enge* clone() const override;

        /** Rescale so Enge(x) -> Enge(scaleFactor*x)
         *
         *  Sets x0 to scaleFactor*x0 and lambda to scaleFactor*lambda
         */
        void rescale(double scaleFactor) override;

        /** Return the value of enge(x+x0) + enge(-x-x0) at some point x */
        double function(double x, int n) const override;

        /** Evaluate a host point without allocating device work buffers.
         * Configuration mirrors are cached until coefficients/indices change.
         */
        void function(double x, int n, double* values) const override;

        /** GPU-aware version of function
         *
         *  @param xView: list of x values
         *  @param n: order of the derivative
         *  @param values: fill values with the derivatives. Should have same
         *         number of rows as xView and n+1 columns. Orders 0..12 are supported;
         *         invalid orders or extents throw OpalException before launching work.
         */
        inline void function(
                const Kokkos::View<double*>& xView, const int n,
                Kokkos::View<double**>& values) const override;

        /** Host side static wrapper */
        static inline void functionHost(
                const EngeConfig& config, const Kokkos::View<double*>& xView, const int n,
                Kokkos::View<double**>& values);

        /// Validate a host request before computing n+1 or constructing output views.
        static inline void validateDerivativeOrder(int n) {
            if (n < 0 || n > EngeConfig::max_derivative)
                throw OpalException("Enge::function", "Derivative order must be between 0 and 12.");
        }

        /** Potentially device-side function call
         *
         *  @param config: configuration
         *  @param x: returns d^n f(x)/dx^n
         *  @param n: order of derivative
         *  @param values: Kokkos::View<double*> of length n+1
         *  @pre Order n is in 0..12 and coefficient/index tables have been initialized.
         *  Invalid device orders abort; host entry points throw before launching work.
         *
         *  I couldn't figure out the type of values so I just make it a template
         */
        template <class C>
        static KOKKOS_INLINE_FUNCTION void functionDevice(
                const EngeConfig& config, double x, int n, C values);

        /** Nominal end length is lambda */
        inline double getEndLength() const override;

        /** Nominal centre length is x0/2 */
        inline double getCentreLength() const override;

        /** Print human-readable version of enge */
        std::ostream& print(std::ostream& out) const override;

        /** Returns the enge polynomial coefficients (a_i) */
        std::vector<double> getCoefficients() const;

        /** Sets the enge polynomial coefficients (a_i) */
        void setCoefficients(std::vector<double> a) {
            config_m.a_m = makeView(a, "efm::Enge::a_m");
            invalidateHostConfig();
        }

        /** Returns the value of lambda */
        double getLambda() const { return config_m.lambda_m; }

        /** Sets the value of lambda */
        inline void setLambda(double lambda) { config_m.lambda_m = lambda; }

        /** Returns the value of x0 */
        double getX0() const { return config_m.x0_m; }

        /** Sets the value of x0 */
        inline void setX0(double x0) { config_m.x0_m = x0; }

        /** Calls setEngeDiffIndices to set the maximum derivative */
        inline void setMaximumDerivative(size_t n) override;

        /** Get a copy of the Config object */
        inline EngeConfig getConfig() const { return config_m; }

        /** Get a copy of the Config object */
        inline void setConfig(EngeConfig& config) { config_m = config; }

        /** Returns the value of the Enge function or its \f$n^{th}\f$ derivative.
         *
         *  Please call setEngeDiffIndices(n) before calling if n > max_index
         */
        static KOKKOS_INLINE_FUNCTION double getEnge(
                const EngeConfig& config, double x, int n, EngeWorkspace& ws);

        /** Setup Enge workspace for an Enge & derivatives calculation */
        static KOKKOS_INLINE_FUNCTION void setupWorkspace(
                const EngeConfig& config, double x, int n, EngeWorkspace& ws);

        /** Returns \f$h(x)\f$ or its \f$n^{th}\f$ derivative.
         *
         *  Here \f$h(x) = a_0 + a_1 x/\lambda + a_2 x^2/lambda^2 + \ldots \f$
         *  Please call setEngeDiffIndices(n) before calling if n > max_index
         */
        static KOKKOS_INLINE_FUNCTION double hN(const EngeConfig& config, double x, int n);

        /** Returns \f$g(x)\f$ or its \f$n^{th}\f$ derivative.
         *
         *  Here \f$g(x) = 1+exp(h(x))\f$.
         *  Please call setEngeDiffIndices(n) before calling if n > max_index
         */
        static KOKKOS_INLINE_FUNCTION double gN(const EngeConfig& config, int n, EngeWorkspace& ws);

        /** Recursively calculate the indices for Enge and H
         *
         *  This will calculate the indices for Enge and H that are required to
         *  calculate the differential up to order n.
         */
        void setEngeDiffIndices(size_t n);
        static void setEngeDiffIndices(size_t n, EngeConfig& config);

        /** Return compact terms for the nth derivative of 1/g, orders 0..12.
         * A row is [coefficient, power of g, power of g', ...]. For multiplicities
         * m_j of derivatives g^(j), k=sum(m_j), sum(j*m_j)=n, the coefficient is
         * \f$(-1)^k n! k! / \prod_j(m_j! (j!)^{m_j})\f$ and the power of g is -k-1.
         * Coefficients use 64-bit integers to include order 12 without overflow.
         */
        inline static EngeDerivativeTable getQIndex(int n);

        /** Return the indices for calculating the nth derivative of g(x) ito h(x) */
        inline static std::vector<std::vector<int>> getHIndex(int n);

        /** Helper function for making 2D view from stl vectors */
        static void copyVectorToView(
                const std::vector<std::vector<int>>& src, Kokkos::View<int**>& dest);
        static void copyVectorToView(
                const EngeDerivativeTable& src, Kokkos::View<std::int64_t**>& dest);

        /** Helper function for making from stl vector */
        static Kokkos::View<double*> makeView(
                const std::vector<double>& src, const std::string& label);
        /** Helper function for making stl vector from view */
        static std::vector<double> makeVector(const Kokkos::View<double*>& src);

        /** Initialise checks that one time setup is done for e.g. initialising
         *  coefficients for derivative calculation
         */
        void initialise() override;

        /** Release Kokkos allocations */
        void finalise() override;
        static void releaseEngeConfig(EngeConfig& config);

    private:
        Enge(const Enge& enge);
        Enge& operator=(const Enge& enge);
        EngeConfig config_m;

        struct HostConfig {
            using HostCoefficients = decltype(Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace{}, Kokkos::View<double*>{}));
            using HostGIndices = decltype(Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace{}, Kokkos::View<std::int64_t**>{}));
            using HostHIndices = decltype(Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace{}, Kokkos::View<int**>{}));
            HostCoefficients a_m;
            double lambda_m = 0.0, x0_m = 0.0;
            HostGIndices gIndices[EngeConfig::max_derivative + 1];
            HostHIndices hIndices[EngeConfig::max_derivative + 1];
        };
        mutable HostConfig host_config_m;
        mutable int host_derivative_m    = -1;
        mutable bool host_coefficients_m = false;
        mutable std::mutex host_mutex_m;
        void invalidateHostConfig() const;
        void prepareHostConfig(int n) const;  // Caller holds host_mutex_m.

        template <class Config, class Values>
        static KOKKOS_INLINE_FUNCTION void functionPoint(
                const Config& config, double x, int n, Values values);
        template <class Config>
        static KOKKOS_INLINE_FUNCTION void setupWorkspaceImpl(
                const Config& config, double x, int n, EngeWorkspace& ws);
        template <class Config>
        static KOKKOS_INLINE_FUNCTION double getEngeImpl(
                const Config& config, double x, int n, EngeWorkspace& ws);
        template <class Config>
        static KOKKOS_INLINE_FUNCTION double hNImpl(const Config& config, double x, int n);
        template <class Config>
        static KOKKOS_INLINE_FUNCTION double gNImpl(const Config& config, int n, EngeWorkspace& ws);

        static void setEngeDiffIndicesQM(size_t n);
        static void setEngeDiffIndicesHM(size_t n);
        static void setEngeDiffIndicesKokkos(size_t n, EngeConfig& config);

        /** Indexes the derivatives of enge in terms of g */
        static std::vector<EngeDerivativeTable> q_m;
        /** Indexes the derivatives of g in terms of h */
        static std::vector<std::vector<std::vector<int>>> h_m;
    };

}  // namespace endfieldmodel

#include "AbsBeamline/EndFieldModel/Enge-inl.icc"

#endif
