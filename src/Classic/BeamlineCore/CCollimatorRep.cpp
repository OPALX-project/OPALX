// ------------------------------------------------------------------------
// $RCSfile: CCollimatorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CCollimatorRep
//   Defines a concrete collimator representation.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/CCollimatorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(CCollimatorRep::*get)() const;
        void (CCollimatorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &CCollimatorRep::getElementLength,
            &CCollimatorRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class CCollimatorRep
// ------------------------------------------------------------------------

CCollimatorRep::CCollimatorRep():
    CCollimator(),
    geometry(0.0)
{}


CCollimatorRep::CCollimatorRep(const CCollimatorRep &right):
    CCollimator(right),
    geometry(right.geometry)
{}


CCollimatorRep::CCollimatorRep(const std::string &name):
    CCollimator(name),
    geometry()
{}


CCollimatorRep::~CCollimatorRep()
{}


ElementBase *CCollimatorRep::clone() const {
    return new CCollimatorRep(*this);
}


Channel *CCollimatorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<CCollimatorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &CCollimatorRep::getField() {
    return field;
}

const NullField &CCollimatorRep::getField() const {
    return field;
}


StraightGeometry &CCollimatorRep::getGeometry() {
    return geometry;
}

const StraightGeometry &CCollimatorRep::getGeometry() const {
    return geometry;
}


ElementImage *CCollimatorRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


/*
double CCollimatorRep::getXsize() const
{
  return xSize;
}

double CCollimatorRep::getYsize() const
{
  return ySize;
}

void CCollimatorRep::setXsize(double size)
{
  INFOMSG("void CCollimatorRep::setXsize(double size) " << xSize << endl;);
  xSize = size;
}

void CCollimatorRep::setYsize(double size)
{
  ySize = size;
}
*/