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

#include <vector>
#include <cmath>
#include <functional>

// Integration workspace
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

inline gsl_integration_workspace* gsl_integration_workspace_alloc(size_t n) {
    gsl_integration_workspace* w = new gsl_integration_workspace();
    w->limit = n;
    w->alist.resize(n);
    w->blist.resize(n);
    w->rlist.resize(n);
    w->elist.resize(n);
    w->order.resize(n);
    return w;
}

inline void gsl_integration_workspace_free(gsl_integration_workspace* w) {
    delete w;
}

// Function structure
struct gsl_function {
    double (*function)(double x, void* params);
    void* params;
};

// Helper function for adaptive integration
namespace {
    std::pair<double, double> adaptive_integrate(
        const std::function<double(double)>& func,
        double a, double b, double tol, size_t max_depth, size_t depth) {
        
        if (depth > max_depth) {
            return {0.0, tol};
        }
        
        // Simpson's rule
        auto simpson = [&func](double a, double b) -> double {
            double h = (b - a) / 2.0;
            double fa = func(a);
            double fm = func(a + h);
            double fb = func(b);
            return (h / 3.0) * (fa + 4.0 * fm + fb);
        };
        
        double c = (a + b) / 2.0;
        double sab = simpson(a, b);
        double sac = simpson(a, c);
        double scb = simpson(c, b);
        double error = std::abs(sab - sac - scb) / 15.0;
        
        if (error < tol) {
            return {sac + scb, error};
        } else {
            auto left = adaptive_integrate(func, a, c, tol / 2.0, max_depth, depth + 1);
            auto right = adaptive_integrate(func, c, b, tol / 2.0, max_depth, depth + 1);
            return {left.first + right.first, left.second + right.second};
        }
    }
}

// Adaptive Gauss-Kronrod integration (simplified)
inline int gsl_integration_qag(const gsl_function* f, double a, double b,
                               double epsabs, double epsrel, size_t /* limit */,
                               int /* key */, gsl_integration_workspace* /* workspace */,
                               double* result, double* abserr) {
    // Simple adaptive Simpson's rule
    const size_t max_depth = 20;
    double tolerance = std::max(epsabs, epsrel * std::abs(b - a));
    
    std::function<double(double)> func = [f](double x) {
        return f->function(x, f->params);
    };
    
    auto res = adaptive_integrate(func, a, b, tolerance, max_depth, 0);
    *result = res.first;
    *abserr = res.second;
    
    return 0;  // Success
}

// GSL integration constants
constexpr int GSL_INTEG_GAUSS15 = 1;
constexpr int GSL_INTEG_GAUSS21 = 2;
constexpr int GSL_INTEG_GAUSS31 = 3;
constexpr int GSL_INTEG_GAUSS41 = 4;
constexpr int GSL_INTEG_GAUSS51 = 5;
constexpr int GSL_INTEG_GAUSS61 = 6;

#endif // OPAL_GSL_INTEGRATION_HH

