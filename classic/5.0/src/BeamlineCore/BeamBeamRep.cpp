// ------------------------------------------------------------------------
// $RCSfile: BeamBeamRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamBeamRep
//   Defines a concrete beam-beam interaction.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/BeamBeamRep.h"
#include "AbsBeamline/ElementImage.h"
#include "BeamlineGeometry/Matrix3D.h"
#include "BeamlineGeometry/Vector3D.h"
#include "Channels/IndexedChannel.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(BeamBeamRep::*get)() const;
        void (BeamBeamRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &BeamBeamRep::getElementLength,
            0
        },
        {
            "CHARGE",
            &BeamBeamRep::getBunchCharge,
            &BeamBeamRep::setBunchCharge
        },
        { 0, 0, 0 }
    };
}


// Class BeamBeamRep
// ------------------------------------------------------------------------

BeamBeamRep::BeamBeamRep():
    BeamBeam(), geometry()
{}


BeamBeamRep::BeamBeamRep(const BeamBeamRep &right):
    BeamBeam(right), geometry()
{}


BeamBeamRep::BeamBeamRep(const string &name):
    BeamBeam(name), geometry()
{}


BeamBeamRep::~BeamBeamRep()
{}


ElementBase *BeamBeamRep::clone() const {
    return new BeamBeamRep(*this);
}


Channel *BeamBeamRep::getChannel(const string &aKey, bool create) {
    static char plane[] = "XYS";
    if(aKey[0] == 'D') {
        for(int i = 0; i < 3; i++) {
            if(aKey[1] == plane[i]) {
                return new IndexedChannel<BeamBeamRep>
                       (*this, &BeamBeamRep::getDisplacement,
                        &BeamBeamRep::setDisplacement, i);
            }
        }
    } else if(aKey[0] == 'R') {
        for(int i = 0; i < 3; i++) {
            if(aKey[1] == plane[i]) {
                for(int j = 0; j < 3; j++) {
                    if(aKey[2] == plane[j]) {
                        int index = 3 * i + j;
                        return new IndexedChannel<BeamBeamRep>
                               (*this, &BeamBeamRep::getMoment,
                                &BeamBeamRep::setMoment, index);
                    }
                }
            }
        }
    }

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<BeamBeamRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


double BeamBeamRep::getBunchCharge() const {
    return Q;
}


const Matrix3D &BeamBeamRep::getBunchMoment() const {
    return sigma;
}


const Vector3D &BeamBeamRep::getBunchDisplacement() const {
    return delta;
}


NullField &BeamBeamRep::getField() {
    return field;
}


const NullField &BeamBeamRep::getField() const {
    return field;
}


NullGeometry &BeamBeamRep::getGeometry() {
    return geometry;
}

const NullGeometry &BeamBeamRep::getGeometry() const {
    return geometry;
}


ElementImage *BeamBeamRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}


void BeamBeamRep::setBunchCharge(double charge) {
    Q = charge;
}


void BeamBeamRep::setBunchMoment(const Matrix3D &moments) {
    sigma = moments;
}


void BeamBeamRep::setBunchDisplacement(const Vector3D &disp) {
    delta = disp;
}


double BeamBeamRep::getMoment(int index) const {
    int i1 = index / 3;
    int i2 = index % 3;
    return sigma(i1, i2);
}


double BeamBeamRep::getDisplacement(int index) const {
    return delta(index);
}


void BeamBeamRep::setMoment(int index, double moment) {
    int i1 = index / 3;
    int i2 = index % 3;
    sigma(i1, i2) = sigma(i2, i1) = moment;
}


void BeamBeamRep::setDisplacement(int index, double disp) {
    delta(index) = disp;
}
