// ------------------------------------------------------------------------
// $RCSfile: CollimatorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CollimatorRep
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

#include "BeamlineCore/CollimatorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(CollimatorRep::*get)() const;
        void (CollimatorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &CollimatorRep::getElementLength,
            &CollimatorRep::setElementLength
        },
        {
            "XSIZE",
            &CollimatorRep::getXsize,
            &CollimatorRep::setXsize
        },
        {
            "YSIZE",
            &CollimatorRep::getYsize,
            &CollimatorRep::setYsize
        },
        { 0, 0, 0 }
    };
}


// Class CollimatorRep
// ------------------------------------------------------------------------

CollimatorRep::CollimatorRep():
    Collimator(),
    geometry(0.0)
{}


CollimatorRep::CollimatorRep(const CollimatorRep &right):
    Collimator(right),
    geometry(right.geometry)
{}


CollimatorRep::CollimatorRep(const std::string &name):
    Collimator(name),
    geometry()
{}


CollimatorRep::~CollimatorRep()
{}


ElementBase *CollimatorRep::clone() const {
    return new CollimatorRep(*this);
}


Channel *CollimatorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<CollimatorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &CollimatorRep::getField() {
    return field;
}

const NullField &CollimatorRep::getField() const {
    return field;
}


StraightGeometry &CollimatorRep::getGeometry() {
    return geometry;
}

const StraightGeometry &CollimatorRep::getGeometry() const {
    return geometry;
}


ElementImage *CollimatorRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


/*
double CollimatorRep::getXsize() const
{
  return xSize;
}

double CollimatorRep::getYsize() const
{
  return ySize;
}

void CollimatorRep::setXsize(double size)
{
  INFOMSG("void CollimatorRep::setXsize(double size) " << xSize << endl;);
  xSize = size;
}

void CollimatorRep::setYsize(double size)
{
  ySize = size;
}
*/
