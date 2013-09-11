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

#include <sstream>
#include <limits>

#include "Fields/EMField.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Structure/LossDataSink.h"
#include "Algorithms/PartBunch.h"
#include "Utilities/OpalException.h"
#include "Utilities/OpalRingSection.h"
#include "AbsBeamline/SBend3D.h"

#include "Elements/OpalRing.h"

// fairly generous tolerance here... maybe too generous? Maybe should be
// user parameter?
const double OpalRing::lengthTolerance_m = 1e-2;
const double OpalRing::angleTolerance_m = 1e-2;

OpalRing::OpalRing(std::string ring)
        : Component(ring), constBField_m(), planarArcGeometry_m(1, 1),
          refPartBunch_m(NULL),
          lossDS_m(NULL),
          beamRInit_m(0.), beamPRInit_m(0.), beamPhiInit_m(0.),
          latticeRInit_m(0.), latticePhiInit_m(0.), latticeThetaInit_m(0.),
          isLocked_m(false), isClosed_m(true),
          symmetry_m(0), cyclHarm_m(0), rfFreq_m(0), current_section_m(0) {
    setRefPartBunch(NULL);
    constBField_m.setB(1.);
}

void OpalRing::accept(BeamlineVisitor& visitor) const {
    visitor.visitOpalRing(*this);
}

OpalRing::OpalRing(const OpalRing& ring)
               : constBField_m(ring.constBField_m),
                 planarArcGeometry_m(ring.planarArcGeometry_m),
                 lossDS_m(NULL),
                 beamRInit_m(ring.beamRInit_m),
                 beamPRInit_m(ring.beamPRInit_m),
                 beamPhiInit_m(ring.beamPhiInit_m),
                 latticeRInit_m(ring.latticeRInit_m),
                 latticePhiInit_m(ring.latticePhiInit_m),
                 latticeThetaInit_m(ring.latticeThetaInit_m),
                 isLocked_m(ring.isLocked_m),
                 isClosed_m(ring.isClosed_m),
                 symmetry_m(ring.symmetry_m),
                 cyclHarm_m(ring.cyclHarm_m),
                 current_section_m(ring.current_section_m) {
    setRefPartBunch(ring.refPartBunch_m);
    if (ring.lossDS_m != NULL)
        throw OpalException("OpalRing::OpalRing(const OpalRing&)",
                 "Can't copy construct LossDataSink so copy constructor fails");
}

OpalRing::~OpalRing() {
    if (lossDS_m != NULL)
        delete lossDS_m;
}

bool OpalRing::apply(const size_t &id, const double &t, Vector_t &E,
                                                                  Vector_t &B) {
  bool flagNeedUpdate =
          apply(refPartBunch_m->R[id], refPartBunch_m->get_centroid(), t, E, B);
  if(flagNeedUpdate) {
      Inform gmsgALL("OPAL ", INFORM_ALL_NODES);
      gmsgALL << getName() << ": particle " << id
              << " at " << refPartBunch_m->R[id]
              << " out of the field map boundary" << endl;
      lossDS_m->addParticle(refPartBunch_m->R[id], refPartBunch_m->P[id], id);
      lossDS_m->save();
      refPartBunch_m->Bin[id] = -1;
  }
  return flagNeedUpdate;
}

bool OpalRing::apply(const size_t &id, const double &t,
                                                       double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);

    if(apply(id, t, Ev, Bv))
        return true;

    E[0] = Ev(0);
    E[1] = Ev(1);
    E[2] = Ev(2);
    B[0] = Bv(0);
    B[1] = Bv(1);
    B[2] = Bv(2);

    return false;
}

bool OpalRing::apply(const Vector_t &R, const Vector_t &centroid,
                                    const double &t, Vector_t &E, Vector_t &B) {
    B = Vector_t(0.0);
    E = Vector_t(0.0);

    OpalRingSection* section = getSectionAt(R);
    return section->getFieldValue(R, centroid, t, E, B);
}

void OpalRing::setLossDataSink(LossDataSink* sink) {
    if (lossDS_m != NULL) {
        delete lossDS_m;
    }
    lossDS_m = sink;
}

void OpalRing::getDimensions(double &zBegin, double &zEnd) const {
    throw OpalException("OpalRing::getDimensions",
                        "Cannot get s-dimension of a ring");
}

void OpalRing::initialise(PartBunch *bunch) {
    online_m = true;
    setRefPartBunch(bunch);
    setLossDataSink(new LossDataSink(getName(), false));
}

void OpalRing::initialise(PartBunch * bunch, double &startField,
                                  double &endField, const double &scaleFactor) {
    initialise(bunch);
}

void OpalRing::finalise() {
    online_m = false;
    setRefPartBunch(NULL);
    setLossDataSink(NULL);
}

void OpalRing::setRefPartBunch(PartBunch* bunch) {
    RefPartBunch_m = bunch; // inherited from Component
    refPartBunch_m = bunch; // private data (obeys style guide)
}

// potential BUG - what if start plane/end plane intersect with non-zero radius
// needs an inner radial and outer radial aperture also to be defined
OpalRingSection* OpalRing::getSectionAt(const Vector_t& r) {
    size_t n_iterations = 0;
    size_t sec_list_size = section_list_m.size();
    while (n_iterations < 2*sec_list_size) {
        // std::cerr << "getSectionAt "
        //          << current_section_m-section_list_m.begin() << std::endl;
        if ((*current_section_m)->isOnOrPastStartPlane(r)) {
            if ((*current_section_m)->isPastEndPlane(r)) {
                ++current_section_m;
                if (current_section_m == section_list_m.end())
                    current_section_m = section_list_m.begin();
            } else {
                return *current_section_m;
            }
        } else {
            if (current_section_m == section_list_m.begin())
                current_section_m = section_list_m.end()-1;
            else
                --current_section_m;
        }
        ++n_iterations;
    }
    std::stringstream r_stream;
    r_stream << r;
    throw OpalException("OpalRing::getSectionAt(const Vector_t& R)",
                        "Failed to locate section at "+r_stream.str());
}

Rotation3D OpalRing::getRotationStartToEnd(Euclid3D delta) const {
    // rotation from start to end in local coordinates
    // Euclid3D/Rotation3D doesnt have a "getAngle" method so we use fairly
    // obscure technique to extract it.
    Vector3D rotationTest(1., 0., 0.);
    rotationTest = delta.getRotation()*rotationTest;
    double deltaAngle = atan2(rotationTest(2), rotationTest(0));
    Rotation3D elementRotation = Rotation3D::ZRotation(deltaAngle);
    return elementRotation;
}

void OpalRing::checkMidplane(Euclid3D delta) const {
    if (fabs(delta.getVector()(1)) > lengthTolerance_m ||
        fabs(delta.getRotation().getAxis()(0)) > angleTolerance_m ||
        fabs(delta.getRotation().getAxis()(2)) > angleTolerance_m) {
        throw OpalException("OpalRing::checkMidplane",
                            std::string("Placement of elements out of the ")+
                            "midplane is not supported by OpalRing");
    }
}

void OpalRing::appendElement(const Component &element) {
    if (isLocked_m) {
        throw OpalException("OpalRing::appendElement",
                            "Attempt to append element "+element.getName()+
                            " when ring is locked");
    }
    // delta is transform from start of bend to end with x, z as horizontal
    Euclid3D delta = element.getGeometry().getTotalTransform();
    checkMidplane(delta);

    OpalRingSection* section = new OpalRingSection();

    Vector3D start_pos(latticeRInit_m*cos(latticePhiInit_m),
                       latticeRInit_m*sin(latticePhiInit_m),
                       0.);
    Vector3D start_norm(sin(latticePhiInit_m+latticeThetaInit_m),
                        cos(latticePhiInit_m+latticeThetaInit_m),
                        0.);
    if (section_list_m.size() > 0) {
        start_pos = section_list_m.back()->getEndPosition();
        start_norm = section_list_m.back()->getEndNormal();
    }

    section->setComponent(dynamic_cast<Component*>(element.clone()));
    section->setStartPosition(start_pos);
    section->setStartNormal(start_norm);
    // translation from start to end in local coordinates (note we Euclid3D
    // works in x-z as horizontal plane; but OPAL-CYCL works with x-y as
    // horizontal plane. So we need to rotate the translation into midplane
    Vector3D elementTranslation = Rotation3D::XRotation(-Physics::pi/2.)
                                                             *delta.getVector();
    /*double r = 2350., pi=3.14159265359;
    double theta = pi/4.;
    std::cerr << "DX (local, no rot) " << delta.getVector()(0) << " " << delta.getVector()(1) << " " << delta.getVector()(2) << std::endl;
    std::cerr << "DX (local) " << elementTranslation(0) << " " << elementTranslation(1) << " " << elementTranslation(2) << std::endl;
    std::cerr << "DX (test local) " << -2.*r*sin(theta/2.)*sin(theta/2.) << " " << 2.*r*sin(theta/2.)*cos(theta/2.) << " " << 0. << std::endl;*/
    // rotation to element start
    double dphi = atan2(-start_norm(0), start_norm(1));
    Rotation3D elementStart = Rotation3D::ZRotation(dphi);
    // std::cerr << "DPHI to start " << dphi << " [rad] " << dphi/Physics::pi*180. << " [deg] " << std::endl;
    // translation from start to end in global coordinates
    Vector3D globalTrans = elementStart*elementTranslation;
    // std::cerr << "DX (global) " << globalTrans(0) << " " << globalTrans(1) << " " << globalTrans(2) << std::endl;

    section->setEndPosition(start_pos + globalTrans);

    // rotation from start to end in global coordinates
    Rotation3D elementRotation = getRotationStartToEnd(delta);
    section->setEndNormal(elementRotation*start_norm);
    // dphi = atan2(-section->getEndNormal()(0), section->getEndNormal()(1));
    section->setComponentPosition(Vector3D(0., 0., 0.));
    section->setComponentOrientation(Vector3D(0., 0., -2.*Physics::pi/8.+dphi));

    section_list_m.push_back(section);

    Inform msg("OPAL ");
    msg << "Added " << element.getName() << " to OpalRing" << endl;
    msg << "  Start position ("
        << section->getStartPosition()(0) << ", "
        << section->getStartPosition()(1) << ", "
        << section->getStartPosition()(2) << ") normal ("
        << section->getStartNormal()(0) << ", "
        << section->getStartNormal()(1) << ", "
        << section->getStartNormal()(2) << "), dphi " << dphi << endl;
    msg << "  End position ("
        << section->getEndPosition()(0) << ", "
        << section->getEndPosition()(1) << ", "
        << section->getEndPosition()(2) << ") normal ("
        << section->getEndNormal()(0) << ", "
        << section->getEndNormal()(1) << ", "
        << section->getEndNormal()(2) << ")" << endl;
    current_section_m = section_list_m.begin();
}

void OpalRing::lockRing() {
    Inform msg("OPAL ");
    if (isLocked_m) {
        throw OpalException("OpalRing::lockRing",
                            "Attempt to lock ring when it is already locked");
    }
    // check for any elements at all
    size_t sectionListSize = section_list_m.size();
    if (sectionListSize == 0) {
        throw OpalException("OpalRing::lockRing",
                            "Failed to find any elements in OpalRing");
    }
    current_section_m = section_list_m.begin();
    // Apply symmetry properties; I think it is fastest to just duplicate
    // elements rather than try to do some angle manipulations when we do field
    // lookups because we do field lookups in Cartesian coordinates in general.
    msg << "Applying symmetry to OpalRing" << endl;  
    for (int i = 1; i < symmetry_m; ++i) {
        for (size_t j = 0; j < sectionListSize; ++j) {
            appendElement(*section_list_m[j]->getComponent());
        }
    }
    // Check that the ring is closed within tolerance
    if (isClosed_m) {
        Vector3D last_pos = section_list_m.back()->getEndPosition();
        Vector3D last_norm = section_list_m.back()->getEndNormal();
        Vector3D first_pos = section_list_m[0]->getStartPosition();
        Vector3D first_norm = section_list_m[0]->getStartNormal();
        for (int i = 0; i < 3; ++i) {
            if (fabs(first_pos(i) - last_pos(i)) > lengthTolerance_m || 
                fabs(first_norm(i) - last_norm(i)) > angleTolerance_m)
                throw OpalException("OpalRing::lockRing",
                                    "Ring is not closed");
        }
        // Exactly close the ring
        section_list_m.back()
                        ->setEndPosition(section_list_m[0]->getStartPosition());
        section_list_m.back()
                        ->setEndNormal(section_list_m[0]->getStartNormal());
        // std::cerr << "IsClosed" << std::endl;
    }
    isLocked_m = true;
    getSectionAt(Vector_t(2349 , -0.000232433 , 2.87833e-12));
    // test_f();  // poor man's unit tests
}

Vector_t convert(Vector3D vec_3d) {
    return Vector_t(vec_3d(0), vec_3d(1), vec_3d(2));
}

Vector3D convert(Vector_t vec_t) {
    return Vector3D(vec_t[0], vec_t[0], vec_t[0]);
}

void OpalRing::test_f() {
    std::ofstream fout("TestOpalRing.dat");
    std::string testpass = "pass";
    for (size_t i = 0; i < section_list_m.size(); ++i) {
        OpalRingSection* sec = section_list_m[i];
        Vector3D start = sec->getStartPosition();
        Vector3D end = sec->getEndPosition();
        bool my_test = getSectionAt(convert(start)) == sec;
        my_test = getSectionAt(convert(end)) == sec;
        my_test = getSectionAt(convert(start+end)/2.) == sec;
        if (!my_test)
            testpass = "fail";
        Vector_t pos, e_loc, b_loc, e_glob, b_glob;
        pos[0] = (sec->getStartPosition()(0)+sec->getEndPosition()(0))/2.;
        pos[1] = (sec->getStartPosition()(1)+sec->getEndPosition()(1))/2.;
        pos[2] = (sec->getStartPosition()(2)+sec->getEndPosition()(2))/2.;
        bool oob_loc = sec->getComponent()->apply(Vector_t(2350., 0., 0.), (0., 0., 0.), 0., e_loc, b_loc);
        bool oob_glob = apply(pos, pos, 0., e_glob, b_glob);
        fout << "Local " << oob_loc << " " << e_loc << " " << b_loc << std::endl;
        fout << "Global " << oob_glob << " " << e_glob << " " << b_glob << std::endl;
    }
    double pi = 3.14159265359, r=2350.;
    for (double phi = -pi/8.; phi < 17.*pi/8.; phi += pi/16.)
      for (double r = 2350.; r < 2351.; r += 50.)
        for (double y = 50.; y < 51.; y += 50.) {
            Vector_t pos(r*cos(phi), r*sin(phi), y);
            Vector_t b, e;
            bool oob_glob = apply(pos, pos, 0., e, b);
            fout << "Position, field " << pos << " " << b << " " << sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]) << std::endl;
        }
    fout << "TestOpalRing " << testpass << std::endl;
    std::cerr << "TestOpalRing " << testpass << std::endl;
}

