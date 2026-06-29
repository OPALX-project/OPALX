#ifndef OPAL_OpalBend_HH
#define OPAL_OpalBend_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalBend
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

// Class OpalBend
// ------------------------------------------------------------------------
/// Base class for all bending magnets.
//  This class factors out the special behaviour for the DOOM interface
//  and the printing in OPAL format, as well as the bend attributes.

class OpalBend : public OpalElement {
public:
    /// The attribute numbers for OPAL bend magnets.
    enum {
        ANGLE = COMMON,  // The bend angle.
        K0,
        K0S,  // The multipole coefficients; must be in this order.
        K1,
        K1S,
        K2,
        K2S,
        K3,
        K3S,
        E1,
        E2,  // The edge angles.
        H1,
        H2,  // The pole face curvatures.
        HGAP,
        FINT,  // The fringing field parameters.
        SLICES,
        STEPSIZE,       // Parameters used to determine slicing.
        FMAPFN,         // File name containing on-axis field.
        GAP,            // Full gap of magnet.
        HAPERT,         // Horizontal aperture of magnet.
        ROTATION,       // Magnet rotation about z axis.
        DESIGNENERGY,   // the design energy of the particles
        GREATERTHANPI,  // Boolean flag set to true if bend angle is greater
                        // than 180 degrees.
        NSLICES,        // The number of slices / steps per element for map tracking
        SIZE            // Total number of attributes.
    };

    /// Exemplar constructor.
    OpalBend(const char* name, const char* help);

    virtual ~OpalBend();

    /// Print the bend magnet.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream&) const;

protected:
    /// Clone constructor.
    OpalBend(const std::string& name, OpalBend* parent);

    /**
     * @brief Return the dipole coefficient implied by an analytic bend angle.
     *
     * The geometry is defined by the total bend angle \f$\theta\f$ and an
     * element-specific field normalization length \f$L_\mathrm{norm}\f$. The
     * corresponding ideal dipole coefficient is therefore
     * \f[
     * k_0 = \frac{\theta}{L_\mathrm{norm}}.
     * \f]
     * `SBEND` uses the historical OPAL chord-to-arc normalization, while
     * `RBEND` uses the integrated Enge field-envelope length.
     *
     * When no positive body length is available, the angle itself is returned
     * as the legacy zero-length fallback used by the historical OPAL parser.
     *
     * @param fieldNormalizationLength Element-specific field normalization length.
     * @param angle Total bend angle in radians.
     * @return Ideal analytic dipole coefficient.
     */
    static double deriveAnalyticDipoleCoefficient(double fieldNormalizationLength, double angle);

    /**
     * @brief Validate the analytic bend-definition inputs.
     *
     * The OPALX analytic `SBEND`/`RBEND` model uses `ANGLE` as the primary
     * geometric quantity. `K0` is treated as a secondary quantity that may be
     * supplied only for consistency checking against
     * \f$k_0 = \theta / L_\mathrm{norm}\f$.
     *
     * @param elementName Name used in exception messages.
     * @param hasAngle True if the input explicitly defines `ANGLE`.
     * @param hasK0 True if the input explicitly defines `K0`.
     * @param fieldNormalizationLength Element-specific field normalization length.
     * @param angle Total bend angle in radians.
     * @param k0Input User-provided dipole coefficient.
     */
    static void validateAnalyticBendDefinition(
            const std::string& elementName, bool hasAngle, bool hasK0,
            double fieldNormalizationLength, double angle, double k0Input);

private:
    // Not implemented.
    OpalBend(const OpalBend&);
    void operator=(const OpalBend&);
};

#endif  // OPAL_OpalBend_HH
