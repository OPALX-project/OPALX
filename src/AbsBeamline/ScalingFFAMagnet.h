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

#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "BeamlineGeometry/Geometry.h"

#ifndef ABSBEAMLINE_ScalingFFAMagnet_H
#define ABSBEAMLINE_ScalingFFAMagnet_H

/** Sector bending magnet with an FFA-style field index and spiral end shape
 *
 *  Note about placement and end field; in order to get a user-defined end field
 *  shape, we do a lookup from user defined end fields, which is pulled from
 *  the EndFieldModel object and done by setupEndField(). Because end fields tie
 *  into the geometry, this has to be done before placements; but it also needs
 *  to be done after parsing has finished. setupEndField() is called by e.g.
 *  ParallelCyclotronTracker just before the object is handed to OpalRing for
 *  placement.
 *
 *  Note about ScalingFFAMagnetConfig; GPU can't handle member data well so I 
 *  make a struct to hold member data which we hand off for the actual field
 *  lookup.
 */

template <class EFM> // EndFieldModel
struct ScalingFFAMagnetConfig {
    size_t maxOrder_m        = 0;
    double tanDelta_m        = 0.;
    double k_m               = 0.;
    double Bz_m              = 0.;
    double r0_m              = 0.;
    double rMin_m            = 0.;  // minimum radius
    double rMax_m            = 0.;  // maximum radius
    double phiStart_m        = 0.;  // offsets this element
    double phiEnd_m          = 0.;  // used for placement of next element
    double azimuthalExtent_m = 0.;  // maximum distance used for field calculation
    double verticalExtent_m  = 0.;  // maximum allowed distance from the midplane
    Vector_t<double, 3> centre_m;
    EFM endField_m;
    std::string endFieldName_m               = "";
    const double fp_tolerance                = 1e-18;
    std::vector<std::vector<double> > dfCoefficients_m;
};

template <class EFM> // EndFieldModel
class ScalingFFAMagnet : public ElementBase {
public:
    /** Construct a new ScalingFFAMagnet
     *
     *  \param name User-defined name of the ScalingFFAMagnet
     */
    explicit ScalingFFAMagnet(const std::string& name);

    /** Destructor - deletes map */
    ~ScalingFFAMagnet();

    /** Inheritable copy constructor */
    ScalingFFAMagnet* clone() const override;

    /** Calculate the field for particles in the container
     *
     *  \param pc the set of particles for which the field is calculated
     */
    void apply(const std::shared_ptr<ParticleContainer_t>& pc) override;

    /** Calculate the field at some arbitrary position
     *
     *  \param R position in the local coordinate system of the bend
     *  \param P not used
     *  \param t not used
     *  \param E not used
     *  \param B calculated magnetic field
     */
    void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /** Calculate the field for particles in the container
     *
     *  \param pc the set of particles for which the field is calculated
     *
     *  This is a static function so that it can call the GPU
     */
    static void getFieldValue(const ScalingFFAMagnetConfig<EFM>& config, const std::shared_ptr<ParticleContainer_t>& pc);

    /** Calculate the field at some arbitrary position in cartesian coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cartesian coordinates defined like (x, y, z)
     *  \param B calculated magnetic field defined like (Bx, By, Bz)
     *  \returns true if particle is outside the field map, else false
     */
    KOKKOS_INLINE_FUNCTION static bool getFieldValue(const ScalingFFAMagnetConfig<EFM>& config, 
                                                     const Vector_t<double, 3>& R,
                                                     Vector_t<double, 3>& B);
    
    /** Calculate the field at some arbitrary position in cylindrical coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cylindrical polar coordinates defined like (r, y, phi)
     *  \param B calculated magnetic field defined like (Br, By, Bphi)
     *  \returns true if particle is outside the field map, else false
     */
    KOKKOS_INLINE_FUNCTION static bool getFieldValueCylindrical(const ScalingFFAMagnetConfig<EFM>& config,
                                                                const Vector_t<double, 3>& R,
                                                                Vector_t<double, 3>& B);

    /** Calculate the field at some arbitrary position in cylindrical coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cylindrical polar coordinates defined like (r, y, phi)
     *  \param B calculated magnetic field defined like (Br, By, Bphi)
     *  \returns true if particle is outside the field map, else false
     */
    bool getFieldValue(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const;

    /** Initialise the ScalingFFAMagnet
     *
     *  \param bunch the global bunch object
     */
    void initialise(PartBunch_t* bunch) override;

    /** Initialise the ScalingFFAMagnet
     *
     *  Sets up the field expansion and the geometry; call after changing any
     *  field parameters
     */
    void initialise();

    /** Finalise the ScalingFFAMagnet - sets bunch to nullptr */
    void finalise() override;

    /** Not implemented */
    void getFieldExtent(double& /*zBegin*/, double& /*zEnd*/) const override {}

    /** Return the cell geometry */
    Geometry& getGeometry() override;

    /** Return the cell geometry */
    const Geometry& getGeometry() const override;

    /** Accept a beamline visitor */
    void accept(BeamlineVisitor& visitor) const override;

    /** Return the field lookup configuration */
    ScalingFFAMagnetConfig<EFM> getConfig() const {return config_m;}

    /** Get tan delta - delta is the spiral angle */
    double getTanDelta() const { return config_m.tanDelta_m; }

    /** Set tan delta - delta is the spiral angle */
    void setTanDelta(double tanDelta) { config_m.tanDelta_m = tanDelta; }

    /** Get the field index k */
    double getFieldIndex() const { return config_m.k_m; }

    /** Set the field index k */
    void setFieldIndex(double k) { config_m.k_m = k; }

    /** Get the dipole constant B_0 */
    double getDipoleConstant() const { return config_m.Bz_m; }

    /** Set the dipole constant B_0 */
    void setDipoleConstant(double Bz) { config_m.Bz_m = Bz; }

    /** Get the radius constant R_0 */
    double getR0() const { return config_m.r0_m; }

    /** Set the radius constant R_0 */
    void setR0(double r0) { config_m.r0_m = r0; }

    /** Get the centre of the sector */
    Vector_t<double, 3> getCentre() const { return config_m.centre_m; }

    /** Set the centre of the sector */
    void setCentre(Vector_t<double, 3> centre) { config_m.centre_m = centre; }

    /** Get the fringe field
     *
     *  Returns the fringe field model; ScalingFFAMagnet retains ownership of the
     *  returned memory.
     */
    EFM getEndField() const { return config_m.endField_m; }

    /** Set the fringe field
     *
     * - endField: the new fringe field; ScalingFFAMagnet takes ownership of the
     *   memory associated with endField.
     */
    void setEndField(EFM endField);

    /** Get the maximum power of y modelled in the off-midplane expansion;
     */
    size_t getMaxOrder() const { return config_m.maxOrder_m; }

    /** Set the maximum power of y modelled in the off-midplane expansion;
     */
    void setMaxOrder(size_t maxOrder) { config_m.maxOrder_m = maxOrder; }

    /** Get the offset of the magnet centre from the start
     */
    double getPhiStart() const { return config_m.phiStart_m; }

    /** Set the offset of the magnet centre from the start
     */
    void setPhiStart(double phiStart) { config_m.phiStart_m = phiStart; }

    /** Get the offset of the magnet end from the start
     */
    double getPhiEnd() const { return config_m.phiEnd_m; }

    /** Set the offset of the magnet end from the start
     */
    void setPhiEnd(double phiEnd) { config_m.phiEnd_m = phiEnd; }

    /** Get the maximum radius
     */
    double getRMin() const { return config_m.rMin_m; }

    /** Set the maximum radius
     */
    void setRMin(double rMin) { config_m.rMin_m = rMin; }

    /** Get the maximum radius
     */
    double getRMax() const { return config_m.rMax_m; }

    /** Set the maximum radius
     */
    void setRMax(double rMax) { config_m.rMax_m = rMax; }

    /** Get the maximum azimuthal displacement from \psi=0
     */
    double getAzimuthalExtent() const { return config_m.azimuthalExtent_m; }

    /** Set the maximum azimuthal displacement from \psi=0
     */
    void setAzimuthalExtent(double azimuthalExtent) { config_m.azimuthalExtent_m = azimuthalExtent; }

    /** Get the maximum vertical displacement from the midplane
     */
    double getVerticalExtent() const { return config_m.verticalExtent_m; }

    /** Set the maximum vertical displacement from the midplane
     */
    void setVerticalExtent(double verticalExtent) { config_m.verticalExtent_m = verticalExtent; }

    /** Return the calculated df coefficients */
    std::vector<std::vector<double> > getDfCoefficients() { return config_m.dfCoefficients_m; }

    /** setupEndField does some end field and geometry set-up
     *
     *  This is normally called just before the magnet is placed; can only set
     *  up the end field after everything has been parsed from input (otherwise
     *  OPAL may not know about an end field model).
     *
     *  sets PhiStart, PhiEnd, AzimuthalExtent and the end field model itself.
     */
    void setupEndField();

    /** Set the end field name.
     *
     *  Called during parsing of the input file; OPAL looks for the endFieldName
     *  when setupEndField() is called.
     */
    void setEndFieldName(std::string name) { config_m.endFieldName_m = name; }

    /** Return the end field name. */
    std::string getEndFieldName() const { return config_m.endFieldName_m; }

private:
    /** Calculate the df coefficients, ready for field generation
     *
     *  Must be called following any update to the the field parameters, in
     *  order for correct field to be calculated.
     */
    void calculateDfCoefficients();

    /** Copy constructor */
    ScalingFFAMagnet(const ScalingFFAMagnet& right);

    ScalingFFAMagnet& operator=(const ScalingFFAMagnet& rhs);
    Geometry planarArcGeometry_m{Geometry::makeSBend(1., 1.)};

    ScalingFFAMagnetConfig<EFM> config_m;
};

template <class EFM>
bool ScalingFFAMagnet<EFM>::getFieldValue(
    const ScalingFFAMagnetConfig<EFM>& config, const Vector_t<double, 3>& R, Vector_t<double, 3>& B) {
    Vector_t<double, 3> pos = R - config.centre_m;
    double r                = std::sqrt(pos[0] * pos[0] + pos[2] * pos[2]);
    double phi              = std::atan2(
            pos[2], pos[0]);  // angle between y-axis and position vector in anticlockwise direction
    Vector_t<double, 3> posCyl({r, pos[1], phi});
    Vector_t<double, 3> bCyl({0., 0., 0.});  // br bz bphi
    bool outOfBounds = getFieldValueCylindrical(config, posCyl, bCyl);
    // this is cartesian coordinates
    B[1] += bCyl[1];
    B[0] += bCyl[0] * std::cos(phi) - bCyl[2] * std::sin(phi);
    B[2] += bCyl[0] * std::sin(phi) + bCyl[2] * std::cos(phi);
    return outOfBounds;
}

template <class EFM>
bool ScalingFFAMagnet<EFM>::getFieldValueCylindrical(
    const ScalingFFAMagnetConfig<EFM>& config, const Vector_t<double, 3>& pos, Vector_t<double, 3>& B) {
    double r   = pos[0];
    double z   = pos[1];
    double phi = pos[2];
    if (r < config.rMin_m || r > config.rMax_m) {
        return true;
    }

    double normRadius = r / config.r0_m;
    double g          = config.tanDelta_m * std::log(normRadius);
    double phiSpiral  = phi - g - config.phiStart_m;
    double h          = std::pow(normRadius, config.k_m) * config.Bz_m;
    if (phiSpiral < -config.azimuthalExtent_m || phiSpiral > config.azimuthalExtent_m) {
        return true;
    }
    if (z < -config.verticalExtent_m || z > config.verticalExtent_m) {
        return true;
    }
    // std::cerr << "ScalingFFAMagnet::getFieldValueCylindrical " << phiSpiral << " "
    //           << config.endField_m->function(phiSpiral, 0) << " " << config.endField_m->getEndLength()
    //           << " " << config.endField_m->getCentreLength()  << std::endl;
    std::vector<double> fringeDerivatives(config.maxOrder_m + 1, 0.);
    for (size_t i = 0; i < fringeDerivatives.size(); ++i) {
        fringeDerivatives[i] = config.endField_m.function(phiSpiral, i);  // d^i_phi f
    }
    for (size_t n = 0; n < config.dfCoefficients_m.size(); n += 2) {
        double f2n = 0;
        Vector_t<double, 3> deltaB;
        for (size_t i = 0; i < config.dfCoefficients_m[n].size(); ++i) {
            f2n += config.dfCoefficients_m[n][i] * fringeDerivatives[i];
        }
        deltaB[1] = f2n * h * std::pow(z / r, n);  // Bz = sum(f_2n * h * (z/r)^2n
        if (config.maxOrder_m >= n + 1) {
            double f2nplus1 = 0;
            for (size_t i = 0;
                 i < config.dfCoefficients_m[n + 1].size() && n + 1 < config.dfCoefficients_m.size(); ++i) {
                f2nplus1 += config.dfCoefficients_m[n + 1][i] * fringeDerivatives[i];
            }
            deltaB[0] = (f2n * (config.k_m - n) / (n + 1) - config.tanDelta_m * f2nplus1) * h
                        * std::pow(z / r, n + 1);  // Br
            deltaB[2] =
                    f2nplus1 * h * std::pow(z / r, n + 1);  // Bphi = sum(f_2n+1 * h * (z/r)^2n+1
        }
        B += deltaB;
    }
    return false;
}

template class ScalingFFAMagnet<endfieldmodel::Tanh>;

#endif
