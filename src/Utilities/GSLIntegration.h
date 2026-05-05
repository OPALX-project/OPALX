//
// GSL Integration compatibility to replace gsl_integration
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPAL_GSL_INTEGRATION_HH
#define OPAL_GSL_INTEGRATION_HH

#include <cmath>
#include <functional>
#include <vector>

/// \see Documentation on https://www.gnu.org/software/gsl/doc/html/integration.html
/// \see Implementation on https://www.gnu.org/software/gsl/
/// \brief Workspace for adaptive integration routines.
struct gsl_integration_workspace {
    size_t limit;
    std::vector<double> alist;
    std::vector<double> blist;
    std::vector<double> rlist;
    std::vector<double> elist;
    std::vector<size_t> order;
    size_t size;
    size_t nrmax;
    size_t i;
    size_t maximum_level;

    gsl_integration_workspace() : limit(0), size(0), nrmax(0), i(0), maximum_level(0) {}
};

/// \see https://www.gnu.org/software/gsl/doc/html/integration.html
/// \brief Allocate integration workspace with limit \f$n\f$.
/// \param n Input: maximum number of subintervals.
/// \return Output: pointer to newly allocated workspace.
inline gsl_integration_workspace* gsl_integration_workspace_alloc(size_t n) {
    gsl_integration_workspace* w = new gsl_integration_workspace();
    w->limit                     = n;
    w->alist.resize(n);
    w->blist.resize(n);
    w->rlist.resize(n);
    w->elist.resize(n);
    w->order.resize(n);
    return w;
}

/// \see https://www.gnu.org/software/gsl/doc/html/integration.html
/// \brief Free a workspace allocated by \c gsl_integration_workspace_alloc.
/// \param w Input: workspace to release (can be null).
inline void gsl_integration_workspace_free(gsl_integration_workspace* w) { delete w; }

/// \see https://www.gnu.org/software/gsl/doc/html/roots.html
/// \brief Function wrapper used by integration routines.
/// \param function Input: callback \f$f(x,\mathrm{params})\f$.
/// \param params Input: opaque parameters passed to \p function.
struct gsl_function {
    double (*function)(double x, void* params);
    void* params;
};

// Helper function for adaptive integration.
namespace {
    std::pair<double, double> adaptive_integrate(
            const std::function<double(double)>& func, double a, double b, double tol,
            size_t max_depth, size_t depth) {
        if (depth > max_depth) {
            return {0.0, tol};
        }

        // Simpson's rule
        auto simpson = [&func](double a, double b) -> double {
            double h  = (b - a) / 2.0;
            double fa = func(a);
            double fm = func(a + h);
            double fb = func(b);
            return (h / 3.0) * (fa + 4.0 * fm + fb);
        };

        double c     = (a + b) / 2.0;
        double sab   = simpson(a, b);
        double sac   = simpson(a, c);
        double scb   = simpson(c, b);
        double error = std::abs(sab - sac - scb) / 15.0;

        if (error < tol) {
            return {sac + scb, error};
        } else {
            auto left  = adaptive_integrate(func, a, c, tol / 2.0, max_depth, depth + 1);
            auto right = adaptive_integrate(func, c, b, tol / 2.0, max_depth, depth + 1);
            return {left.first + right.first, left.second + right.second};
        }
    }
}  // namespace

/// \see https://www.gnu.org/software/gsl/doc/html/integration.html
/// \brief Adaptive integration on \f$[a,b]\f$ using Simpson refinement.
/// \details Computes an approximation of \f$\int_a^b f(x)\,dx\f$ with absolute
/// tolerance \p epsabs and relative tolerance \p epsrel.
/// \param f Input: integrand.
/// \param a Input: lower limit.
/// \param b Input: upper limit.
/// \param epsabs Input: absolute error tolerance.
/// \param epsrel Input: relative error tolerance.
/// \param limit Input: maximum subintervals (unused).
/// \param key Input: quadrature rule selector (unused).
/// \param workspace Input: workspace (unused).
/// \param result Output: estimated integral value.
/// \param abserr Output: estimated absolute error.
/// \return Output: 0 on success.
inline int gsl_integration_qag(
        const gsl_function* f, double a, double b, double epsabs, double epsrel, size_t /* limit */,
        int /* key */, gsl_integration_workspace* /* workspace */, double* result, double* abserr) {
    // Simple adaptive Simpson's rule
    const size_t max_depth = 20;
    double tolerance       = std::max(epsabs, epsrel * std::abs(b - a));

    std::function<double(double)> func = [f](double x) {
        return f->function(x, f->params);
    };

    auto res = adaptive_integrate(func, a, b, tolerance, max_depth, 0);
    *result  = res.first;
    *abserr  = res.second;

    return 0;  // Success
}

/// \brief Integration rule identifiers (accepted but not used in this implementation).
constexpr int GSL_INTEG_GAUSS15 = 1;
constexpr int GSL_INTEG_GAUSS21 = 2;
constexpr int GSL_INTEG_GAUSS31 = 3;
constexpr int GSL_INTEG_GAUSS41 = 4;
constexpr int GSL_INTEG_GAUSS51 = 5;
constexpr int GSL_INTEG_GAUSS61 = 6;

#endif  // OPAL_GSL_INTEGRATION_HH
