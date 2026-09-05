// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_COMPENSATED_SUM_H
#define OPAL_COMPENSATED_SUM_H

namespace compensated {
    /**
     * @brief Kahan addition retaining increments below the current sum's spacing.
     * The represented value is sum - correction. Do not reassociate these operations
     * (e.g. with fast-math). Both values must be retained when copying a tracked state.
     */
    inline void add(double increment, double& sum, double& correction) {
        const double adjusted = increment - correction;
        const double next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }

    /// Difference of two compensated values, preserving their small residuals.
    inline double difference(double left, double leftCorrection, double right,
                             double rightCorrection) {
        return (left - right) - (leftCorrection - rightCorrection);
    }
}
#endif
