// ------------------------------------------------------------------------
// $RCSfile: CorrectorWrapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CorrectorWrapper
//   Defines a representation for a modified corrector.
//
// ------------------------------------------------------------------------
// Class category: ComponentWrappers
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "ComponentWrappers/CorrectorWrapper.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class CorrectorWrapper
//   A CorrectorWrapper represents a unique instance of a corrector magnet
//   in the accelerator model.
// ------------------------------------------------------------------------

CorrectorWrapper::CorrectorWrapper(Corrector *corr):
    Corrector(*corr),
    itsDesign(corr),
    itsError(),
    tempField() {
    shareFlag = false;
}


CorrectorWrapper::CorrectorWrapper(const CorrectorWrapper &rhs):
    Corrector(rhs),
    itsDesign(rhs.itsDesign),
    itsError(rhs.itsError),
    tempField() {
    shareFlag = false;
}


CorrectorWrapper::~CorrectorWrapper()
{}


void CorrectorWrapper::accept(BeamlineVisitor &visitor) const {
    visitor.visitCorrectorWrapper(*this);
}


ElementBase *CorrectorWrapper::clone() const {
    CorrectorWrapper *cw = new CorrectorWrapper(*this);
    cw->itsDesign = dynamic_cast<Corrector *>(itsDesign->clone());
    return cw;
}


ElementBase *CorrectorWrapper::copyStructure() {
    if(isSharable()) {
        return this;
    } else {
        CorrectorWrapper *cw = new CorrectorWrapper(*this);
        cw->itsDesign = dynamic_cast<Corrector *>(itsDesign->copyStructure());
        return cw;
    }
}


void CorrectorWrapper::makeSharable() {
    shareFlag = true;
    itsDesign->makeSharable();
}


StraightGeometry &CorrectorWrapper::getGeometry() {
    return itsDesign->getGeometry();
}


const StraightGeometry &CorrectorWrapper::getGeometry() const {
    return itsDesign->getGeometry();
}


Corrector::Plane CorrectorWrapper::getPlane() const {
    return itsDesign->getPlane();
}


ElementBase::ElementType CorrectorWrapper::getType() const {
    return CORRECTORWRAPPER;
}


const Corrector &CorrectorWrapper::getDesign() const {
    return *itsDesign;
}


Corrector &CorrectorWrapper::getDesign() {
    return *itsDesign;
}


BDipoleField &CorrectorWrapper::errorField() const {
    return itsError;
}


ElementBase *CorrectorWrapper::makeFieldWrapper() {
    return this;
}


ElementBase *CorrectorWrapper::removeFieldWrapper() {
    return &*itsDesign;
}


const ElementBase *CorrectorWrapper::removeFieldWrapper() const {
    return &*itsDesign;
}


ElementBase *CorrectorWrapper::removeWrappers() {
    return &*itsDesign;
}


const ElementBase *CorrectorWrapper::removeWrappers() const {
    return &*itsDesign;
}


BDipoleField &CorrectorWrapper::getField() {
    tempField = itsDesign->getField();
    tempField.addField(itsError);
    return tempField;
}


const BDipoleField &CorrectorWrapper::getField() const {
    tempField = itsDesign->getField();
    tempField.addField(itsError);
    return tempField;
}