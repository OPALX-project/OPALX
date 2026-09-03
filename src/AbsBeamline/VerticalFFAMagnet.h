//
// Header file for VerticalFFAMagnet Component
//
// Copyright (c) 2019 Chris Rogers
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#include "AbsBeamline/ElementBase.h"
#include "BeamlineGeometry/Geometry.h"
#include "PartBunch/PartBunch.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"

#ifndef ABSBEAMLINE_VerticalFFAMagnet_H
#define ABSBEAMLINE_VerticalFFAMagnet_H

/** VerticalFFA field calculation data
 */
template <class EFM>
struct VerticalFFAMagnetConfig {
    size_t maxOrder_m   = 0;
    double k_m          = 0.;
    double Bz_m         = 0.;
    double zNegExtent_m = 0.;  // extent downwards from the midplane
    double zPosExtent_m = 0.;  // extent upwards from the midplane
    double halfWidth_m  = 0.;  // extent in either +x or -x
    double bbLength_m   = 0.;
    EFM endField_m;
    std::vector<std::vector<double> > dfCoefficients_m;
};

/** Bending magnet with an exponential dependence on field in the vertical plane
 *
 *  VerticalFFAMagnet makes a rectangular bending magnet with a dipole field
 *  that has a dependence like B0 exp(mz)
 */
template <class EFM> // EndFieldModel
class VerticalFFAMagnet : public ElementBase {
public:
    /** Construct a new VerticalFFAMagnet
     *
     *  \param name User-defined name of the VerticalFFAMagnet
     */
    explicit VerticalFFAMagnet(const std::string& name);

    /** Destructor - deletes the field */
    ~VerticalFFAMagnet();

    /** Inheritable copy constructor */
    ElementBase* clone() const override;

    /** Calculate the field at the position of the ith particle
     *
     *  \param i index of the particle event; field is calculated at this
     *         position
     *  \param t time at which the field is to be calculated
     *  \param E calculated electric field - always 0 (no E-field)
     *  \param B calculated magnetic field
     */
    inline void apply(const std::shared_ptr<ParticleContainer_t>& pc) override;

    inline void apply(
            const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B);

    /** Calculate the field at some arbitrary position
     *
     *  \param R position in the local coordinate system of the magnet
     *  \param P not used
     *  \param t not used
     *  \param E not used
     *  \param B calculated magnetic field
     */
    inline void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /** Calculate the field at some arbitrary position in cartesian coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cartesian coordinates defined like (x, y, z)
     *  \param B calculated magnetic field defined like (Bx, By, Bz)
     *  \returns true if particle is outside the field map, else false
     */
    bool getFieldValue(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const;

    /** Calculate the field for particles in the container
     *
     *  \param pc the set of particles for which the field is calculated
     *
     *  This is a static function so that it can call the GPU
     */
    static void getFieldValue(const VerticalFFAMagnetConfig<EFM>& config, const std::shared_ptr<ParticleContainer_t>& pc);

    /** Calculate the field at some arbitrary position in cartesian coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cartesian coordinates defined like (x, y, z)
     *  \param B calculated magnetic field defined like (Bx, By, Bz)
     *  \returns true if particle is outside the field map, else false
     */
    KOKKOS_INLINE_FUNCTION static bool getFieldValue(const VerticalFFAMagnetConfig<EFM>& config, 
                                                     const Vector_t<double, 3>& R,
                                                     Vector_t<double, 3>& B);

    /** Initialise the VerticalFFAMagnet
     *
     *  \param bunch the global bunch object (but not used)
     */
    void initialise(PartBunch_t* bunch) override;

    /** Initialise the VerticalFFAMagnet
     *
     *  Sets up the field expansion and the geometry; call after changing any
     *  field parameters
     */
    void initialise();

    /** Finalise the VerticalFFAMagnet - sets bunch to nullptr */
    void finalise() override;

    /** Return false - VerticalFFAMagnet is a straight magnet
     *
     *  Nb: the VerticalFFAMagnet geometry is straight even though trajectories
     *      are not
     */

    /** Not implemented */
    void getFieldExtent(double& /*zBegin*/, double& /*zEnd*/) const override {}

    /** Return the cell geometry */
    Geometry& getGeometry() override;

    /** Return the cell geometry */
    const Geometry& getGeometry() const override;

    /** Accept a beamline visitor */
    void accept(BeamlineVisitor& visitor) const override;

    /** Get the fringe field
     *
     *  Returns the fringe field model; VerticalFFAMagnet retains ownership of
     *  the returned memory.
     */
    EFM getEndField() const { return config_m.endField_m; }

    /** Set the fringe field
     *
     * - endField: the new fringe field; VerticalFFAMagnet takes ownership of
     *   the memory associated with endField.
     */
    void setEndField(EFM endField);

    /** Get the maximum power of x used in the off-midplane expansion;
     */
    size_t getMaxOrder() const { return config_m.maxOrder_m; }

    /** Set the maximum power of x used in the off-midplane expansion;
     */
    void setMaxOrder(size_t maxOrder);

    /** Get the centre field at z=0 */
    double getB0() const { return config_m.Bz_m; }

    /** Set the centre field at z=0 */
    void setB0(double Bz) { config_m.Bz_m = Bz; }

    /** Get the field index */
    double getFieldIndex() const { return config_m.k_m; }  // units are [m^{-1}]

    /** Set the field index */
    void setFieldIndex(double index) { config_m.k_m = index; }

    /** Get the maximum extent below z = 0 */
    double getNegativeVerticalExtent() const { return config_m.zNegExtent_m; }

    /** Set the maximum extent below z = 0 */
    inline void setNegativeVerticalExtent(double negativeExtent);

    /** Get the maximum extent above z = 0 */
    double getPositiveVerticalExtent() const { return config_m.zPosExtent_m; }

    /** set the maximum extent above z = 0 */
    inline void setPositiveVerticalExtent(double positiveExtent);

    /** Get the length of the bounding box (centred on magnet centre) */
    double getBBLength() const { return config_m.bbLength_m; }

    /** Set the length of the bounding box (centred on magnet centre) */
    void setBBLength(double bbLength) { config_m.bbLength_m = bbLength; }

    /** Get the full width of the bounding box (centred on magnet centre) */
    double getWidth() const { return config_m.halfWidth_m * 2.; }

    /** Set the full width of the bounding box (centred on magnet centre) */
    void setWidth(double width) { config_m.halfWidth_m = width / 2; }

    /** Get the coefficients used for the field expansion
     *
     *  B_y is given by
     *     sum_n B_0 exp(ky) f_n x^n
     *  where
     *     f_n = sum_k c_{nk} partial_k f_0
     *
     *  Returns a vector of vectors, like c[n][k]. The expansion for the other
     *  field elements can be related back to c[n][k] (see elsewhere for details).
     */
    inline std::vector<std::vector<double> > getDfCoefficients() const;

private:
    void calculateDfCoefficients();

    /** Copy constructor */
    VerticalFFAMagnet(const VerticalFFAMagnet& right);

    VerticalFFAMagnet& operator=(const VerticalFFAMagnet& rhs);
    Geometry straightGeometry_m{Geometry::makeStraight(1.)};
    VerticalFFAMagnetConfig<EFM> config_m;
};

template class VerticalFFAMagnet<endfieldmodel::Tanh>;

template <class EFM>
void VerticalFFAMagnet<EFM>::setNegativeVerticalExtent(double negativeExtent) {
    config_m.zNegExtent_m = negativeExtent;
}

template <class EFM>
void VerticalFFAMagnet<EFM>::setPositiveVerticalExtent(double positiveExtent) {
    config_m.zPosExtent_m = positiveExtent;
}

template <class EFM>
void VerticalFFAMagnet<EFM>::apply(const std::shared_ptr<ParticleContainer_t>& /*pc*/) {}

template <class EFM>
void VerticalFFAMagnet<EFM>::apply(
        const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    auto Rview                              = pc->R.getView();
    auto Pview                              = pc->P.getView();
    const Vector_t<double, 3> R             = Rview(i);
    const Vector_t<double, 3> P             = Pview(i);
    apply(R, P, t, E, B);
}

template <class EFM>
void VerticalFFAMagnet<EFM>::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double&,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B) {
    getFieldValue(R, B);
}

template <class EFM>
std::vector<std::vector<double> > VerticalFFAMagnet<EFM>::getDfCoefficients() const {
    return config_m.dfCoefficients_m;
}

template <class EFM>
void VerticalFFAMagnet<EFM>::getFieldValue(const VerticalFFAMagnetConfig<EFM>& config, 
                                           const std::shared_ptr<ParticleContainer_t>& pc) {
    const Kokkos::View<Vector_t<double, 3>*> R = pc->R.getView();
    const Kokkos::View<Vector_t<double, 3>*> B = pc->B.getView();
    const size_t count = pc->getLocalNum();
        Kokkos::parallel_for(
            "VerticalFFAMagnet<>::getFieldValue()", count, KOKKOS_LAMBDA(const size_t i) {
                getFieldValue(config, R(i), B(i));
            });
}


template <class EFM>
bool VerticalFFAMagnet<EFM>::getFieldValue(const VerticalFFAMagnetConfig<EFM>& config_m, 
                   const Vector_t<double, 3>& R,
                   Vector_t<double, 3>& B) {
    if (std::abs(R[0]) > config_m.halfWidth_m || R[2] < 0. || R[2] > config_m.bbLength_m || R[1] < -config_m.zNegExtent_m
        || R[1] > config_m.zPosExtent_m) {
        return true;
    }
    std::vector<double> fringeDerivatives(config_m.maxOrder_m + 2, 0.);
    double zRel = R[2] - config_m.bbLength_m / 2.;  // z relative to centre of magnet
    for (size_t i = 0; i < fringeDerivatives.size(); ++i) {
        fringeDerivatives[i] = config_m.endField_m.function(zRel, i);  // d^i_phi f
    }

    std::vector<double> x_n(config_m.maxOrder_m + 1);  // x^n
    x_n[0] = 1.;                              // x^0
    for (size_t i = 1; i < x_n.size(); ++i) {
        x_n[i] = x_n[i - 1] * R[0];
    }

    // note that the last element is always 0, because dfCoefficients_m is
    // of size maxOrder_m+1. This leads to better Maxwellianness in testing.
    std::vector<double> f_n(config_m.maxOrder_m + 2, 0.);
    std::vector<double> dz_f_n(config_m.maxOrder_m + 1, 0.);
    for (size_t n = 0; n < config_m.dfCoefficients_m.size(); ++n) {
        const std::vector<double>& coefficients = config_m.dfCoefficients_m[n];
        for (size_t i = 0; i < coefficients.size(); ++i) {
            f_n[n] += coefficients[i] * fringeDerivatives[i];
            dz_f_n[n] += coefficients[i] * fringeDerivatives[i + 1];
        }
    }
    double bref = config_m.Bz_m * exp(config_m.k_m * R[1]);
    B[0]        = 0.;
    B[1]        = 0.;
    B[2]        = 0.;
    for (size_t n = 0; n < x_n.size(); ++n) {
        B[0] += bref * f_n[n + 1] * (n + 1) / config_m.k_m * x_n[n];
        B[1] += bref * f_n[n] * x_n[n];
        B[2] += bref * dz_f_n[n] / config_m.k_m * x_n[n];
    }
    return false;
}

#endif
