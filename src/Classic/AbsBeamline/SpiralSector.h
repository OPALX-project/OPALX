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

#include "Fields/BMultipoleField.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "AbsBeamline/Component.h"

#ifndef ABSBEAMLINE_SPIRALSECTOR_H
#define ABSBEAMLINE_SPIRALSECTOR_H

/** Sector bending magnet with an FFAG-style field index and spiral end shape
 */

class SpiralSector : public Component {
  public:
    /** Construct a new SpiralSector
     *
     *  \param name User-defined name of the SpiralSector
     */
    explicit SpiralSector(const std::string &name);

    /** Destructor - deletes map */
    ~SpiralSector();

    /** Inheritable copy constructor */
    ElementBase* clone() const;

    /** Calculate the field at the position of the ith particle
     *
     *  \param i index of the particle event; field is calculated at this
     *         position
     *  \param t time at which the field is to be calculated
     *  \param E calculated electric field - always 0 (no E-field)
     *  \param B calculated magnetic field
     *  \returns true if particle is outside the field map
     */
    bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    /** Calculate the field at some arbitrary position
     *
     *  \param R position in the local coordinate system of the bend
     *  \param P not used
     *  \param t time at which the field is to be calculated
     *  \param E calculated electric field - always 0 (no E-field)
     *  \param B calculated magnetic field
     *  \returns true if particle is outside the field map, else false
     */
    bool apply(const Vector_t &R, const Vector_t &P, const double &t,
               Vector_t &E, Vector_t &B);

    /** Calculate the field at some arbitrary position in cylindrical coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cylindrical polar coordinates defined like (r, y, phi)
     *  \param P not used
     *  \param t not used (field is time independent)
     *  \param E not used (no E-field)
     *  \param B calculated magnetic field defined like (Br, By, Bphi)
     *  \returns true if particle is outside the field map, else false
     */
    bool applyCylindrical(const Vector_t &R, const Vector_t &P,
                    const double &t, Vector_t &E, Vector_t &B);

     /** Initialise the SpiralSector
      *
      *  \param bunch the global bunch object
      *  \param startField not used
      *  \param endField not used
      */
    void initialise(PartBunch *bunch, double &startField, double &endField);

     /** Finalise the SpiralSector - sets bunch to NULL */
    void finalise();

    /** Return true - SpiralSector always bends the reference particle */
    inline bool bends() const;

    /** Not implemented */
    void getDimensions(double &zBegin, double &zEnd) const {}

    /** Return the cell geometry */
    BGeometryBase& getGeometry();

    /** Return the cell geometry */
    const BGeometryBase& getGeometry() const;

    /** Return a dummy (0.) field value (what is this for?) */
    EMField &getField();

    /** Return a dummy (0.) field value (what is this for?) */
    const EMField &getField() const;

    /** Accept a beamline visitor */
    void accept(BeamlineVisitor& visitor) const;

    /** Get tan delta - delta is the spiral angle */
    Vector_t getTanDelta() const {return tanDelta_m;}

    /** Set tan delta - delta is the spiral angle */
    void setTanDelta(double tanDelta) {tanDelta_m = tanDelta;}

    /** Get the field index k */
    double getFieldIndex() const {return k_m;}

    /** Set the field index k */
    void setFieldIndex(double k) {k_m = k;}

    /** Get the dipole constant B_0 */
    double getDipoleConstant() const {return Bz_m;}

    /** Set the dipole constant B_0 */
    void setDipoleConstant(double Bz) {Bz_m = Bz;}

    /** Get the radius constant R_0 */
    double getR0() const {return r0_m;}

    /** Set the radius constant R_0 */
    void setR0(double r0) {r0_m = r0;}

    /** Get the centre of the sector */
    Vector_t getCentre() const {return centre_m;}

    /** Set the centre of the sector */
    void setCentre(Vector_t centre) {centre_m = centre;}

    /** Get the fringe field
     *
     *  Returns the fringe field model; SpiralSector retains ownership of the
     *  returned memory.
     */
    endfieldmodel::EndFieldModel* getEndField() const {return endField_m;}

    /** Set the fringe field
      * 
      * - endField: the new fringe field; SpiralSector takes ownership of the
      *   memory associated with endField.
      */
    void setEndField(endfieldmodel::EndFieldModel* endField);

    /** Get the maximum pole modelled in the off-midplane expansion; 0 is dipole; 2 is quadrupole; etc */
    size_t getMaxOrder() const {return maxOrder_m;}

    /** Set the maximum pole modelled in the off-midplane expansion; 0 is dipole; 2 is quadrupole; etc */
    void setMaxOrder(size_t maxOrder) {maxOrder_m = maxOrder;}

    /** Calculate the df coefficients, ready for field generation
     *
     *  Must be called following any update to the the field parameters, in
     *  order for correct field to be calculated.
     */
    void calculateDfCoefficients();

    /** Return the calculated df coefficients */
    std::vector<std::vector<double> > getDfCoefficients() {return dfCoefficients_m;}

  private:
    /** Copy constructor */
    SpiralSector(const SpiralSector &right);

    SpiralSector& operator=(const SpiralSector& rhs);
    PlanarArcGeometry planarArcGeometry_m;
    BMultipoleField dummy;

    size_t maxOrder_m = 0;
    double tanDelta_m = 0.;
    double k_m;
    double Bz_m;
    double r0_m;
    Vector_t centre_m;
    endfieldmodel::EndFieldModel* endField_m = NULL;
    const double fp_tolerance = 1e-18;
    std::vector<std::vector<double> > dfCoefficients_m;
};

#endif

