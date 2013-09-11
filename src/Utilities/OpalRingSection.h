/* 
 *  Copyright (c) 2012, Chris Rogers
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

#ifndef OPAL_RING_SECTION_H
#define OPAL_RING_SECTION_H

#include <vector>

#include "AbsBeamline/Component.h"

/** \class OpalRingSection Component placement handler in ring geometry
 *
 *  OpalRingSection handles placement of a component when it is placed in a ring
 *  geometry. Here the primary index for section placement is azimuthal angle.
 *  OpalRingSection assumes that ring objects occupy a space defined by a start
 *  plane and another end plane.
 *
 *  All vectors should be in cartesian coordinates (x,y,z); with z being the
 *  axis of the ring.
 *
 *  \param component_m component that the OpalSection wraps - this is a borrowed
 *                   reference (OpalRingSection does not own the memory)
 *  \param startPosiion_m position of the centre of the start face in cylindrical
 *                   polar coordinates
 *  \param startOrientation_m vector normal to the start face pointing towards
 *                   the field map
 *  \param endPosition_m position of the centre of the end face in cylindrical
 *                   polar coordinates
 *  \param endOrientation_m vector normal to the end face pointing away from the
 *                   field map
 *  \param componentPosition_m field map position V
 *  \param componentRotation_m field map rotation R
 *
 *  So to be clear field maps are calculated in U_local = V+R*U_global local
 *  coordinate system where R is rotation matrix and V is componentPosition_m.
 *  Return field values are returned like B_global = R^{-1}*B_local
 */

class OpalRingSection {
public:
    /** Construct a ring section - positions, orientations etc default to 0.
     */
    OpalRingSection();

    /** Destructor - does nothing*/
    ~OpalRingSection();

    /** Return true if pos is on or past start plane
     *
     *  \param   pos position to test
     *  \returns true if phi(pos) >= phi(start plane), after taking account of
     *           the position and rotation of the start face
     */
    bool isOnOrPastStartPlane(const Vector_t& pos);

    /** Return true if pos is past end plane
     *
     *  \returns true if phi > phi(end plane), after taking account of the
     *           position and rotation of the start face
     */
    bool isPastEndPlane(const Vector_t& pos);

    /** Return field value in global coordinate system
     *
     *  \param pos Position in global Cartesian coordinates
     *  \param t time in lab frame
     *  \param centroid Not sure what this is
     *  \param E Vector to be filled with electric field values; will always
     *         overwrite with 0. before filling field
     *  \param B Vector to be filled with magnetic field values; will always
     *         overwrite with 0. before filling field
     *  \returns true if pos is outside of the field bounding box, else false
     */
    bool getFieldValue(const Vector_t& pos, const Vector_t& centroid,
                       const double& t, Vector_t& E, Vector_t& B);

    /** Set the component wrapped by OpalRingSection
     *
     *  This borrows the Component* pointer (caller is responsible for cleanup)
     */
    inline void setComponent(Component* component) {component_m = component;}

    /** Get the component wrapped by OpalRingSection
     *
     *  Component* is not owned by caller or OpalRingSection
     */
    inline Component* getComponent() const {return component_m;}

    /** Set a position on the plane of the section start */
    inline void setStartPosition(Vector3D pos) {startPosition_m = convert(pos);}

    /** Get a position on the plane of the section start */
    inline Vector3D getStartPosition() const {return convert(startPosition_m);}

    /** Set the normal vector to the section start plane */
    inline void setStartNormal(Vector3D orientation) {startOrientation_m = convert(orientation);} 

    /** Get the normal vector to the section start plane */
    inline Vector3D getStartNormal() const {return convert(startOrientation_m);}

    /** Set a position on the section end plane */
    inline void setEndPosition(Vector3D pos) {endPosition_m = convert(pos);}

    /** Get a position on the section end plane */
    inline Vector3D getEndPosition() const {return convert(endPosition_m);}

    /** Set the normal vector to the section end plane  */
    inline void setEndNormal(Vector3D orientation) {endOrientation_m = convert(orientation);}

    /** Get the normal vector to the section end plane */
    inline Vector3D getEndNormal() const {return convert(endOrientation_m);}

    /** Set the displacement for the component relative to the section start */
    inline void setComponentPosition(Vector3D position) {componentPosition_m = convert(position);}

    /** Get the displacement for the component relative to the section start */
    inline Vector3D getComponentPosition() const {return convert(componentPosition_m);}

    /** Set the rotation for the component relative to the section start */
    inline void setComponentOrientation(Vector3D orientation);

    /** Get the rotation for the component relative to the section start */
    inline Vector3D getComponentOrientation() const {return convert(componentOrientation_m);}

private:
    inline Vector_t convert(const Vector3D& vec) const;
    inline Vector3D convert(const Vector_t& vec) const;
    inline void rotate(Vector_t& vector) const;
    inline void rotate_back(Vector_t& vector) const;

    Component* component_m;

    Vector_t componentPosition_m;
    Vector_t componentOrientation_m;
    Rotation3D componentOrientation3D_m;

    Vector_t startPosition_m;
    Vector_t startOrientation_m;

    Vector_t endPosition_m;
    Vector_t endOrientation_m;

    void updateComponentOrientation();
    double sin0_m;
    double cos0_m;
    double sin1_m;
    double cos1_m;
    double sin2_m;
    double cos2_m;
};

Vector_t OpalRingSection::convert(const Vector3D& vec) const {
    return Vector_t(vec(0), vec(1), vec(2));
}

Vector3D OpalRingSection::convert(const Vector_t& vec) const {
    return Vector3D(vec(0), vec(1), vec(2));
}

typedef std::vector<OpalRingSection*> RingSectionList;

inline void OpalRingSection::setComponentOrientation(Vector3D orientation) {
    componentOrientation_m = convert(orientation);
    updateComponentOrientation();
}

inline void OpalRingSection::rotate(Vector_t& vector) const {
    const Vector_t temp(vector);
    vector(0) = +cos2_m * temp(0) + sin2_m * temp(1);
    vector(1) = -sin2_m * temp(0) + cos2_m * temp(1);
/* Out-of-plane rotations are commented for now
    vector(0) = (c0 * c2) * temp(0) +
                (-s0 * s1 * c2 + c1 * s2) * temp(1) +
                (-s0 * c1 * c2 - s1 * s2) * temp(2);
    vector(1) = (-c0 * s2) * temp(0) +
                (s0 * s1 * s2 + c1 * c2) * temp(1) +
                (s0 * c1 * s2 - s1 * c2) * temp(2);
    vector(2) = (s0) * temp(0) +
                (c0 * s1) * temp(1) +
                (c0 * c1) * temp(2);
*/
}

inline void OpalRingSection::rotate_back(Vector_t& vector) const {
    const Vector_t temp(vector);
    vector(0) = +cos2_m * temp(0) - sin2_m * temp(1);
    vector(1) = +sin2_m * temp(0) + cos2_m * temp(1);
/* Out-of-plane rotations are commented for now
    vector(0) = (c0 * c2) * temp(0) +
                (-s0 * s1 * c2 + c1 * s2) * temp(1) +
                (-s0 * c1 * c2 - s1 * s2) * temp(2);
    vector(1) = (-c0 * s2) * temp(0) +
                (s0 * s1 * s2 + c1 * c2) * temp(1) +
                (s0 * c1 * s2 - s1 * c2) * temp(2);
    vector(2) = (s0) * temp(0) +
                (c0 * s1) * temp(1) +
                (c0 * c1) * temp(2);
*/
}

#endif //OPAL_RING_SECTION_H

