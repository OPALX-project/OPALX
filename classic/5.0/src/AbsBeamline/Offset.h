/* 
 *  Copyright (c) 2014, Chris Rogers
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


#ifndef CLASSIC_ABSBEAMLINE_Offset_HH
#define CLASSIC_ABSBEAMLINE_Offset_HH

#include "Fields/EMField.h"
#include "BeamlineGeometry/Euclid3DGeometry.h"
#include "AbsBeamline/Component.h"

class Fieldmap;

/** @Class Offset
 *
 *  Enables user to define a placement, either in global coordinates or in the
 *  coordinate system of the previously placed object
 *
 *  @param _end_position final position of the offset
 *  @param _end_direction normal vector to entry face
 *  @param _is_local parameter is True if everything is in the coordinate system
 *  of the last placed object. It is expected that everything will be in a local 
 *  coordinate system before tracking begins (this is expected by, for example,
 *  OpalRing). This is expected to be set by the visit function in e.g.
 *  ParallelCyclotronTracker.
 *  @param geometry the geometry that OpalRingSection uses to do displacements
 *  This has to be a pointer because SRotatedGeometry
 *  does not have an defined assignment op. SRotatedGeometry does not have
 *  a assignment op because BGeometryBase does not have an assignment op...
 *  and there is no way (I think) to reassign a previously assigned variable
 *  in C++ without using =
 *  @param wrappedGeometry Something to do with the reference handling in
 *  SRotatedGeometry - looks like it doesn't make a copy but holds a reference, 
 *  so unless I keep that alive in Offset I get memory errors (seg fault)
 *  @float_tolerance bends() and operator==(...) use float_tolerance when
 *  evaluating equality between doubles.
 */
class Offset : public Component {
  public:
    /** Constructor sets everything to 0., makes a default geometry with
     *  everything set to 0.
     */
    explicit Offset(const std::string &name);

    /** Default constructor sets everything to 0., makes a default geometry
     *  with everything set to 0.
     */
    Offset();

    /** Copy constructor; deep copies geometry_m; all other
     *  stuff is copied as well
     */
    Offset(std::string name, const Offset &);

    /** Copy constructor; deep copies geometry_m; all other
     *  stuff is copied as well
     */
    Offset(const Offset &);

    /** Assignment operator deep copies geometry and wrappedGeometry; all other
     *  stuff is copied as well
     */
    Offset& operator=(const Offset &);

    /** Factory method to make an offset in Cylindrical coordinates local to the
     *  end of the previous element
     *   - name name of the offset
     *   - theta_in angle between the previous element and the displacement
     *     vector
     *   - theta_out angle between the displacement vector and the next element
     *   - displacement length of the displacement vector
     */
    static Offset localCylindricalOffset(std::string name,
                                         double theta_in,
                                         double theta_out,
                                         double displacement);

    /** Factory method to make an offset in global cylindrical polar coordinates
     *   - name name of the offset
     *   - radius_out radius of the end of the offset
     *   - phi_out azimuthal angle of the end of the offset
     *   - theta_out angle relative to the tangent of the end of the offset
     *  Call updateGeometry(Vector_t, Vector_t) using the end of the previous
     *  element before placement; Offset will convert to the local coordinate
     *  system.
     */
    static Offset globalCylindricalOffset(std::string name,
                                          double radius_out,
                                          double phi_out,
                                          double theta_out);

    /** Factory method to make an offset in cartesian coordinates local to the
     *  end of the previous element
     *   - name name of the offset
     *   - end_position position of the end of the offset
     *   - end_direction direction of the end of the offset
     */
    static Offset localCartesianOffset(std::string name,
                                       Vector_t end_position,
                                       Vector_t end_direction);

    /** Factory method to make an offset in global cartesian coordinates
     *   - name name of the offset
     *   - end_position position of the end of the offset
     *   - end_direction direction of the end of the offset
     *  Call updateGeometry(Vector_t, Vector_t) using the end of the previous
     *  element before placement; Offset will convert to the local coordinate
     *  system.
     */
    static Offset globalCartesianOffset(std::string name,
                                        Vector_t end_position,
                                        Vector_t end_direction);

    /** deletes geometry and wrappedGeometry */
    ~Offset();

    /** Apply visitor to Offset.
     *
     *  Sets ring radius
     */
    void accept(BeamlineVisitor &) const;

    /** Just calls the copy constructor on *this */
    ElementBase* clone() const;

    /** Return false
     *
     *  No bound checking is performed and no field value is calculated as this
     *  element does not have an associated field or beam pipe aperture.
     */
    bool apply(const size_t &i, const double &t, double E[], double B[]);

    /** Return false
     *
     *  No bound checking is performed and no field value is calculated as this
     *  element does not have an associated field or beam pipe aperture.
     */
    bool apply(const size_t &i, const double &t,
                       Vector_t &E, Vector_t &B);

    /** Return false
     *
     *  No bound checking is performed and no field value is calculated as this
     *  element does not have an associated field or beam pipe aperture.
     */
    bool apply(const Vector_t &R, const Vector_t &centroid,
                       const double &t, Vector_t &E, Vector_t &B);

    /** Returns true if either input angle or output angle are greater than
     *  float_tolerance
     */
    bool bends() const;

    void initialise(PartBunch *bunch, double &startField,
                            double &endField, const double &scaleFactor);
    void finalise();
    void getDimensions(double &zBegin, double &zEnd) const {}

    void setEndPosition(Vector_t position);
    Vector_t getEndPosition() const;

    void setEndDirection(Vector_t direction);
    Vector_t getEndDirection() const;

    /** Set to true if stored coordinates are in the local coordinate system of
     *  the last placed object
     */
    void setIsLocal(bool isLocal);

    /** Returns true if stored coordinates are in the local coordinate system of
     *  the last placed object
     */
    bool getIsLocal() const;

    Euclid3DGeometry& getGeometry();
    const Euclid3DGeometry& getGeometry() const;
    void updateGeometry(Vector_t startPosition, Vector_t startDirection);
    void updateGeometry();

    /// Not implemented - throws OpalException
    EMField &getField();
    /// Not implemented - throws OpalException
    const EMField &getField() const;
    /** Calculate the angle between vectors on the midplane
     *
     *  Returns theta in domain -pi, pi. A positive angle means a rotation
     *  anticlockwise from vec1 to vec2.
     *
     *  Throws an OpalException if vec1, vec2 are not in the midplane i.e.
     *  non-zero z.
     */
    static double getTheta(Vector_t vec1, Vector_t vec2);

    /** Rotate vec anticlockwise by angle theta about z axis; return the rotated
     *  vector
     */
    static Vector_t rotate(Vector_t vec, double theta);

    static double float_tolerance;
  private:
    Vector_t _end_position;
    Vector_t _end_direction;
    bool     _is_local;
    // The offset's geometry.
    Euclid3DGeometry* geometry_m;
};

/** Return true if off1 is equal to off2 within static floating point tolerance
 */
bool operator==(const Offset& off1, const Offset& off2);

/** Return not of operator == */
bool operator!=(const Offset& off1, const Offset& off2);

/** Print Offset off1 to the ostream */
std::ostream& operator<<(std::ostream& out, const Offset& off1);

#endif  // CLASSIC_ABSBEAMLINE_Offset_HH
