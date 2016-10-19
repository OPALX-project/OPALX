// ------------------------------------------------------------------------
// $RCSfile: AlignWrapper.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignWrapper
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class AlignWrapper
// ------------------------------------------------------------------------


void AlignWrapper::accept(BeamlineVisitor &visitor) const {
    visitor.visitAlignWrapper(*this);
}


AlignWrapper *AlignWrapper::clone() const {
    return new AlignWrapper(*this);
}


ElementBase *AlignWrapper::copyStructure() {
    if(isSharable()) {
        return this;
    } else {
        AlignWrapper *wrap = new AlignWrapper(*this);
        wrap->itsElement = itsElement->copyStructure();
        return wrap;
    }
}


void AlignWrapper::makeSharable() {
    shareFlag = true;
    itsElement->makeSharable();
}


Euclid3D AlignWrapper::getEntranceTransform() const {
    if(itsOffset.isIdentity()) {
        return itsOffset;
    } else {
        Euclid3D frame = itsElement->getEntranceFrame();
        return frame.inverse() * itsOffset * frame;
    }
}


Euclid3D AlignWrapper::getExitTransform() const {
    if(itsOffset.isIdentity()) {
        return itsOffset;
    } else {
        Euclid3D frame = itsElement->getExitFrame();
        return frame * Inverse(itsOffset) * Inverse(frame);
    }
}


ElementBase *AlignWrapper::getElement() const {
    return &*itsElement;
}


void AlignWrapper::setElement(ElementBase *elem) {
    itsElement = elem;
}


BGeometryBase &AlignWrapper::getGeometry() {
    return itsElement->getGeometry();
}


const BGeometryBase &AlignWrapper::getGeometry() const {
    return itsElement->getGeometry();
}


ElementBase::ElementBase::ElementType AlignWrapper::getType() const {
    return itsElement->getType();
}


Euclid3D &AlignWrapper::offset() const {
    return itsOffset;
}


ElementBase *AlignWrapper::makeAlignWrapper() {
    return this;
}


ElementBase *AlignWrapper::makeFieldWrapper() {
    itsElement = itsElement->makeFieldWrapper();
    itsElement->setName(itsElement->getName());
    return this;
}


ElementBase *AlignWrapper::removeAlignWrapper() {
    return &*itsElement;
}


const ElementBase *AlignWrapper::removeAlignWrapper() const {
    return &*itsElement;
}


ElementBase *AlignWrapper::removeFieldWrapper() {
    return new AlignWrapper(itsElement->removeFieldWrapper());
}


const ElementBase *AlignWrapper::removeFieldWrapper() const {
    return new AlignWrapper(itsElement->removeFieldWrapper());
}


ElementBase *AlignWrapper::removeWrappers() {
    return itsElement->removeWrappers();
}


const ElementBase *AlignWrapper::removeWrappers() const {
    return itsElement->removeWrappers();
}


AlignWrapper::AlignWrapper(const AlignWrapper &rhs):
    ElementBase(rhs), itsElement(rhs.itsElement), itsOffset(rhs.itsOffset) {
    shareFlag = false;
}


AlignWrapper::~AlignWrapper()
{}


AlignWrapper::AlignWrapper(ElementBase *elem):
    ElementBase(elem->getName()), itsElement(elem), itsOffset() {
    shareFlag = false;
}