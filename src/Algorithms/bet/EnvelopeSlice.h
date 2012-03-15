#ifndef _ENVELOPE_SLICE_H
#define _ENVELOPE_SLICE_H

#include <cmath>

enum EnvelopeSliceParameter {
    /// z slice position [m]
    SLI_z    = 0,
    /// beta normalized velocity (total) [-]
    SLI_beta = 1,
    /// x beam size x (rms) [m]
    SLI_x    = 2,
    /// px beam divergence x [rad]
    SLI_px   = 3,
    /// y beam size y (rms) [m]
    SLI_y    = 4,
    /// py beam divergence y [rad]
    SLI_py   = 5,
    /// X0 position centroid x [m]
    SLI_x0   = 6,
    /// pX0 angular deflection centriod x
    SLI_px0  = 7,
    /// Y0 position centroid y [m]
    SLI_y0   = 8,
    /// pY0 angular deflection centriod y
    SLI_py0  = 9,
    /// number of slice parameters
    SLNPAR
};

/**
 * @brief class encpasulating an envelope slice
 *
 * Important note: The index must match definition in RK routine!
 */
class EnvelopeSlice {

public:

    /**
     * @brief parameters
     */
    double p[SLNPAR];
    /**
     * @brief parameters before last integration step (backup)
     */
    double p_old[SLNPAR];

    /**
     * @brief flag holding validity of slice
     */
    int valid;
    /**
     * @brief backup copy of valid
     */
    int was_valid;

    EnvelopeSlice();

    ~EnvelopeSlice() {}

    /**
     * @brief set the memory buffer to a new location
     *
     * @param b pointer to new buffer location
     */
    void setBuffer(double *b);

    /**
     * @brief calculates \gamma from \beta
     *
     * @return \gamma of the slice
     */
    inline double gamma() { return sqrt(1.0 / (1.0 - p[SLI_beta] * p[SLI_beta])); }

    /**
     * @brief backup the slice
     */
    void backup();
    /**
     * @brief restore a slice
     */
    void restore();

    /**
     * @brief check if a slice has changed
     *
     * @return type of change
     */
    int check();
    /**
     * @brief check if a slice is still valid
     *
     * @return 1 if valid, 0 otherwise
     */
    inline int is_valid() { return valid; }
};

typedef EnvelopeSlice *EnvelopeSliceP;

#endif

