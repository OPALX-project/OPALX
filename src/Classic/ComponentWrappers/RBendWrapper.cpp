// ------------------------------------------------------------------------
// $RCSfile: RBendWrapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBendWrapper
//   Defines a representation for a modified RBend.
//
// ------------------------------------------------------------------------
// Class category: ComponentWrappers
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:54 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "ComponentWrappers/RBendWrapper.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class RBendWrapper
// ------------------------------------------------------------------------

RBendWrapper::RBendWrapper(RBend *ideal):
    RBend(*ideal),
    itsDesign(ideal),
    itsError(),
    tempField() {
    shareFlag = false;
}


RBendWrapper::RBendWrapper(const RBendWrapper &rhs):
    RBend(rhs),
    itsDesign(rhs.itsDesign),
    itsError(rhs.itsError),
    tempField() {
    shareFlag = false;
}


RBendWrapper::~RBendWrapper()
{}


void RBendWrapper::accept(BeamlineVisitor &visitor) const {
    visitor.visitRBendWrapper(*this);
}


ElementBase *RBendWrapper::clone() const {
    RBendWrapper *rbw = new RBendWrapper(*this);
    rbw->itsDesign = dynamic_cast<RBend *>(itsDesign->clone());
    return rbw;
}


ElementBase *RBendWrapper::copyStructure() {
    if(isSharable()) {
        return this;
    } else {
        RBendWrapper *rbw = new RBendWrapper(*this);
        rbw->itsDesign = dynamic_cast<RBend *>(itsDesign->copyStructure());
        return rbw;
    }
}


BMultipoleField &RBendWrapper::errorField() const {
    return itsError;
}


void RBendWrapper::makeSharable() {
    shareFlag = true;
    itsDesign->makeSharable();
}


BMultipoleField &RBendWrapper::getField() {
    tempField = itsDesign->getField();
    tempField.addField(itsError);
    return tempField;
}


const BMultipoleField &RBendWrapper::getField() const {
    tempField = itsDesign->getField();
    tempField.addField(itsError);
    return tempField;
}


RBendGeometry &RBendWrapper::getGeometry() {
    return itsDesign->getGeometry();
}


const RBendGeometry &RBendWrapper::getGeometry() const {
    return itsDesign->getGeometry();
}


double RBendWrapper::getB() const {
    getField();
    return tempField.getNormalComponent(1);
}


double RBendWrapper::getEntryFaceRotation() const {
    return itsDesign->getEntryFaceRotation();
}


double RBendWrapper::getExitFaceRotation()  const {
    return itsDesign->getExitFaceRotation();
}


double RBendWrapper::getEntryFaceCurvature() const {
    return itsDesign->getEntryFaceCurvature();
}


double RBendWrapper::getExitFaceCurvature()  const {
    return itsDesign->getExitFaceCurvature();
}


double RBendWrapper::getSlices()  const {
    return itsDesign->getSlices();
}


double RBendWrapper::getStepsize()  const {
    return itsDesign->getStepsize();
}


ElementBase::ElementType RBendWrapper::getType() const {
    return RBENDWRAPPER;
}


const RBend &RBendWrapper::getDesign() const {
    return *itsDesign;
}


RBend &RBendWrapper::getDesign() {
    return *itsDesign;
}


ElementBase *RBendWrapper::makeFieldWrapper() {
    return this;
}


ElementBase *RBendWrapper::removeFieldWrapper() {
    return &*itsDesign;
}


const ElementBase *RBendWrapper::removeFieldWrapper() const {
    return &*itsDesign;
}


ElementBase *RBendWrapper::removeWrappers() {
    return &*itsDesign;
}


const ElementBase *RBendWrapper::removeWrappers() const {
    return &*itsDesign;
}