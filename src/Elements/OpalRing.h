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

#ifndef OPALRING_H
#define OPALRING_H

#include <string>

#include "Physics/Physics.h"
#include "Algorithms/Tracker.h"
#include "Fields/ConstBField.h"
#include "AbsBeamline/Component.h"
#include "AbsBeamline/SBend3D.h"

#include "BeamlineCore/RBendRep.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"

#include "Utilities/OpalRingSection.h"
#include "Utilities/OpalField.h"
#include "Utilities/OpalException.h"

class LossDataSink;
class PartBunch;
class FieldMap;

/** @class OpalRing describes a ring type geometry for tracking
 *
 *  OpalRing describes a ring type geometry for tracking. OpalRing provides the
 *  necessary interfaces for e.g. OPAL-CYCL to track through the ring geometry,
 *  while enabling the user to add arbitrary field elements in a closed
 *  geometry.
 *
 *  OpalRing uses similar routines to OPAL-T OpalBeamline class to set up
 *  geometry; note that as OPAL-CYCL places the beam in an x-y geometry we place
 *  in an x-y geometry for OpalRing (i.e. the axis of the ring is by default the
 *  z-axis). So far only placement of elements in the midplane is supported. It
 *  is not possible to give vertical displacements or rotations, or add elements
 *  that might create vertical displacements and rotations. 
 *
 *  Also aim to maintain backwards compatibility with Cyclotron (i.e. use
 *  ParallelCyclotronTracker)
 */
class OpalRing : public Component {
  public:
    /** Constructor
     *
     *  \param ring Name of the ring as defined in the input file
     */
    OpalRing(std::string ring);

    /** Copy constructor
     *
     *  Can't copy LossDataSink so throw exception if this is set
     */
    OpalRing(const OpalRing& ring);

    /** Destructor - deletes lossDS_m if not NULL */
    virtual ~OpalRing();

    /** Overwrite data in array E and B with electric and magnetic fields and
     *  flag particles outside of the ring aperture
     *
     *  @param id index of item in RefPartBunch_m - particle bunch
     *  @param t time
     *  @param E array where electric field vector is stored - any
     *         existing data is overwritten
     *  @param B array where magnetic field vector is stored - any
     *         existing data is overwritten
     *
     *  @returns false if particle is outside of the field map apertures, else
     *  true. If particle is off the field maps, then set flag on the particle
     *  "Bin" data to -1
     */
    virtual bool apply(const size_t &id, const double &t, double E[],
                       double B[]);

    /** Overwrite data in vector E and B with electric and magnetic field
     *
     *  @param i index of item in RefPartBunch_m - particle bunch
     *  @param t time
     *  @param E array where electric field vector is stored - any
     *         existing data is overwritten
     *  @param B array where magnetic field vector is stored - any
     *         existing data is overwritten
     *
     *  @returns false if particle is outside of the field map apertures, else
     *  true. If particle is off the field maps, then set flag on the particle
     *  "Bin" data to -1 and store in LossDataSink
     */
    virtual bool apply(const size_t &id, const double &t, Vector_t &E,
                       Vector_t &B);

    /** Overwrite data in vector E and B with electromagnetic field at point R
     *
     *  @param R 3 vector position at which the field is found in Cartesian
     *         coordinates (i.e. x, y, z with z=vertical)
     *  @param centroid unknown, but not used - bunch mean maybe?
     *  @param t time
     *  @param E vector where electric field vector will be stored - any
     *         existing data is overwritten
     *  @param B vector where magnetic field vector will be stored - any
     *         existing data is overwritten
     *
     *  @returns false if particle is still in the ring aperture, else true
     */
    virtual bool apply(const Vector_t &R, const Vector_t &centroid,
                       const double &t, Vector_t &E, Vector_t &B);

    /** Initialise the OpalRing
     *
     *  @param bunch the particle bunch. OpalRing borrows this pointer (caller
     *         owns memory)
     *  @param startField - not used
     *  @param endField - not used
     *  @param scaleFactor - not used
     */
    virtual void initialise(PartBunch *bunch, double &startField,
                            double &endField, const double &scaleFactor);

    /** Initialise the OpalRing - set the bunch and allocate a new LossDataSink
     *
     *  @param bunch the particle bunch. OpalRing borrows this pointer (caller
     *         owns memory)
     */
    virtual void initialise(PartBunch *bunch);

    /** Clean up the OpalRing
     *
     *  OpalRing relinquishes RefPartBunch pointer and deletes LossDataSink
     */
    virtual void finalise();

    /** Returns true - OpalRing is assumed to bend particles, being a ring */
    virtual bool bends() const {return true;}

    /** Accept the BeamlineVisitor
     *
     *  Just calls visitOpalRing function on the visitor. I guess the point of
     *  this function is that it enables us to store a pointer to the visitor
     *  object or something
     */
    virtual void accept(BeamlineVisitor& visitor) const;

    /** Not implemented - always throws an exception */
    virtual void getDimensions(double &zBegin, double &zEnd) const;

    /** Inherited copy constructor */
    virtual ElementBase* clone() const {return new OpalRing(*this);}

    /** Add element to the ring
     *
     *  Add element to the ring. Elements are assumed to occupy a region of
     *  space defined by a (flat) plane at the start and a plane at the end,
     *  both infinite in extent. The position and rotation of these planes are
     *  defined according to the Component geometry.
     *
     *  Borrows a reference to the element (relinquished when the ring is
     *  deleted).
     *
     *  Throws an exception if the geometry would bend the element out of the
     *  midplane (elements out of midplane are not yet supported).
     */
    void appendElement(const Component &element);

    /** Not implemented */
    virtual EMField &getField() {return constBField_m.getField();}

    /** Not implemented */
    virtual const EMField &getField() const {return constBField_m.getField();}

    /** Not implemented */
    virtual PlanarArcGeometry &getGeometry() {return planarArcGeometry_m;}

    /** Not implemented */
    virtual const PlanarArcGeometry &getGeometry() const {return planarArcGeometry_m;}

    /** Set LossDataSink to sink.
     *
     *  @param sink The LossDataSink. OpalRing takes ownership of memory
     *         allocated to sink
     */
    void setLossDataSink(LossDataSink* sink);

    /** Get pointer to lossDataSink.
     *
     *  OpalRing still owns the memory to which lossDataSink points.
     */
    PartBunch* getLossDataSink() const;

    /** Set RefPartBunch to bunch.
     *
     *  @param sink The Bunch. OpalRing borrows memory allocated to bunch.
     *
     *  Note for compliance with style guide and compatibility with parent two
     *  pointer to RefPartBunch are stored; this keeps them aligned
     */
    void setRefPartBunch(PartBunch* bunch);

    /** Get pointer to RefPartBunch from the bunch.
     *
     *  OpalRing does not own this memory (so neither does caller).
     */
    PartBunch* getRefPartBunch() const;

    /** Set the harmonic number for RF (number of bunches in the ring) */
    void setHarmonicNumber(double cyclHarm) {cyclHarm_m = cyclHarm;}

    /** Get the harmonic number for RF (number of bunches in the ring) */
    double getHarmonicNumber() {return cyclHarm_m;}
    // note this is not a const method to follow parent

    /** Set the nominal RF frequency */
    void setRFFreq(double rfFreq) {rfFreq_m = rfFreq;}

    /** Get the nominal RF frequency */
    double getRFFreq() const {return rfFreq_m;}

    /** Set the initial beam radius */
    void setBeamRInit(double rInit) {beamRInit_m = rInit;}

    /** Get the initial beam radius */
    double getBeamRInit() const {return beamRInit_m;}

    /** Set the initial beam azimuthal angle */
    void setBeamPhiInit(double phiInit) {beamPhiInit_m = phiInit;}

    /** Get the initial beam azimuthal angle */
    double getBeamPhiInit() const {return beamPhiInit_m;}

    /** Set the initial beam radial momentum */
    void setBeamPRInit(double pRInit) {beamPRInit_m = pRInit;}

    /** Get the initial beam radial momentum */
    double getBeamPRInit() const {return beamPRInit_m;}

    /** Set the initial element's radius */
    void setLatticeRInit(double rInit) {latticeRInit_m = rInit;}

    /** Get the initial element's radius */
    double getLatticeRInit() const {return latticeRInit_m;}

    /** Set the initial element's azimuthal angle */
    void setLatticePhiInit(double phiInit) {latticePhiInit_m = phiInit;}

    /** Get the initial  element's azimuthal angle */
    double getLatticePhiInit() const {return latticePhiInit_m;}
 
    /** Set the first element's horizontal angle
     *
     *  Set the angle in the ring plane with respect to the tangent vector
     */
    void setLatticeThetaInit(double thetaInit) {latticeThetaInit_m = thetaInit;}

    /** Get the first element's horizontal angle
     *
     *  Get the angle in the ring plane with respect to the tangent vector
     */
    double getLatticeThetaInit() const {return latticeThetaInit_m;}
   
    /** Set the rotational symmetry of the ring (number of cells) */
    void setSymmetry(double symmetry) {symmetry_m = symmetry;}

    /** Get the rotational symmetry of the ring (number of cells) */
    double getSymmetry() const {return symmetry_m;}
   
    /** Set flag for closure checking */
    void setIsClosed(bool isClosed) {isClosed_m = isClosed;}

    /** Get flag for closure checking */
    double getIsClosed() const {return isClosed_m;}

    /** Lock the ring
     *
     *  Lock the ring; apply closure checks and symmetry properties as required.
     *  Sets isLocked_m to true. New elements can no longer be added (as it may
     *  break the symmetry/bound checking)
     */
    void lockRing();

  private:
    /** Get the section at position pos */
    OpalRingSection* getSectionAt(const Vector_t& pos);

    /** Disabled */
    OpalRing();

    /** Disabled */
    OpalRing& operator=(const OpalRing& ring);

    void checkMidplane(Euclid3D delta) const;
    Rotation3D getRotationStartToEnd(Euclid3D delta) const;

    RBendRep constBField_m;
    PlanarArcGeometry planarArcGeometry_m;

    // points to same location as RefPartBunch_m on the child bunch, but we
    // rename to keep in line with style guide
    //
    // OpalRing borrows this memory
    PartBunch *refPartBunch_m;

    // store for particles out of the aperture 
    //
    // OpalRing owns this memory
    LossDataSink *lossDS_m;

    // initial position of the beam
    double beamRInit_m;
    double beamPRInit_m;
    double beamPhiInit_m;

    // position, orientation of the first lattice element
    double latticeRInit_m;
    double latticePhiInit_m;
    double latticeThetaInit_m;

    // Ring is locked - new elements cannot be added
    bool isLocked_m;

    // Set to false to enable a non-circular ring (e.g. some weird spiral
    // geometry)
    bool isClosed_m;

    // number of cells/rotational symmetry of the ring
    int symmetry_m;

    // rf harmonic number
    double cyclHarm_m;

    // nominal rf frequency
    double rfFreq_m;

    // element wrappers
    RingSectionList::iterator current_section_m;
    RingSectionList section_list_m; // vector of OpalRingSection (Component placement)

    // tolerance on checking for geometry consistency
    static const double lengthTolerance_m;
    static const double angleTolerance_m;
};

#endif //#ifndef OPALRING_H

