// ------------------------------------------------------------------------
// $RCSfile: SBendWrapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SBendWrapper
//   Defines a representation for a modified SBend.
//
// ------------------------------------------------------------------------
// Class category: ComponentWrappers
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:54 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "ComponentWrappers/SBendWrapper.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class SBendWrapper
// ------------------------------------------------------------------------

SBendWrapper::SBendWrapper(SBend *ideal):
    SBend(*ideal),
    itsDesign(ideal),
    itsError(),
    tempField() {
    shareFlag = false;
}


SBendWrapper::SBendWrapper(const SBendWrapper &rhs):
    SBend(rhs),
    itsDesign(rhs.itsDesign),
    itsError(rhs.itsError),
    tempField() {
    shareFlag = false;
}


SBendWrapper::~SBendWrapper()
{}


void SBendWrapper::accept(BeamlineVisitor &visitor) const {
    visitor.visitSBendWrapper(*this);
}


ElementBase *SBendWrapper::clone() const {
    SBendWrapper *sbw = new SBendWrapper(*this);
    sbw->itsDesign = dynamic_cast<SBend *>(itsDesign->clone());
    return sbw;
}


ElementBase *SBendWrapper::copyStructure() {
    if(isSharable()) {
        return this;
    } else {
        SBendWrapper *sbw = new SBendWrapper(*this);
        sbw->itsDesign = dynamic_cast<SBend *>(itsDesign->copyStructure());
        return sbw;
    }
}


BMultipoleField &SBendWrapper::errorField() const {
    return itsError;
}


void SBendWrapper::makeSharable() {
    shareFlag = true;
    itsDesign->makeSharable();
}


BMultipoleField &SBendWrapper::getField() {
    tempField = itsDesign->getField();
    tempField.addField(itsError);
    return tempField;
}


const BMultipoleField &SBendWrapper::getField() const {
    tempField = itsDesign->getField();
    tempField.addField(itsError);
    return tempField;
}


PlanarArcGeometry &SBendWrapper::getGeometry() {
    return itsDesign->getGeometry();
}


const PlanarArcGeometry &SBendWrapper::getGeometry() const {
    return itsDesign->getGeometry();
}


double SBendWrapper::getB() const {
    getField();
    return tempField.getNormalComponent(1);
}


double SBendWrapper::getEntryFaceRotation() const {
    return itsDesign->getEntryFaceRotation();
}


double SBendWrapper::getExitFaceRotation()  const {
    return itsDesign->getExitFaceRotation();
}


double SBendWrapper::getEntryFaceCurvature() const {
    return itsDesign->getEntryFaceCurvature();
}


double SBendWrapper::getExitFaceCurvature()  const {
    return itsDesign->getExitFaceCurvature();
}


double SBendWrapper::getSlices()  const {
    return itsDesign->getSlices();
}


double SBendWrapper::getStepsize()  const {
    return itsDesign->getStepsize();
}


ElementBase::ElementType SBendWrapper::getType() const {
    return SBENDWRAPPER;
}


const SBend &SBendWrapper::getDesign() const {
    return *itsDesign;
}


SBend &SBendWrapper::getDesign() {
    return *itsDesign;
}


ElementBase *SBendWrapper::makeFieldWrapper() {
    return this;
}


ElementBase *SBendWrapper::removeFieldWrapper() {
    return &*itsDesign;
}


const ElementBase *SBendWrapper::removeFieldWrapper() const {
    return &*itsDesign;
}


ElementBase *SBendWrapper::removeWrappers() {
    return &*itsDesign;
}


const ElementBase *SBendWrapper::removeWrappers() const {
    return &*itsDesign;
}