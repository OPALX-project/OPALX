// ------------------------------------------------------------------------
// $RCSfile: CorrectionBase.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Classes: CorrectionBase
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:21 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Tables/CorrectionBase.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Monitor.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/IdealMapper.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalHKicker.h"
#include "Elements/OpalKicker.h"
#include "Elements/OpalVKicker.h"
#include "Fields/BDipoleField.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"


// Nested class CorrectionBase::Row
// ------------------------------------------------------------------------


// Class CorrectionBase
// ------------------------------------------------------------------------

CorrectionBase::Row::Row(ElementBase *elem, int occur):
    FlaggedElmPtr(elem, occur), orbit(), matrix()
{}


CorrectionBase::Row::Row(FlaggedElmPtr const &rhs):
    FlaggedElmPtr(rhs), orbit(), matrix()
{}


CorrectionBase::Row::~Row(void)
{}


// Class CorrectionBase
// ------------------------------------------------------------------------
// Abstract base class for all orbit correction commands.
// Factors out all common behaviour for these algorithms.

CorrectionBase::CorrectionBase(int size, const char *name, const char *help):
    Action(size, name, help) {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "The beam line to be Micadoed");
    itsAttr[BEAM] = Attributes::makeString
                    ("BEAM", "The beam to be used", "UNNAMED_BEAM");
    itsAttr[RANGE] = Attributes::makeRange
                     ("RANGE", "The range in the lattice");
}


CorrectionBase::CorrectionBase(const string &name, CorrectionBase *parent):
    Action(name, parent)
{}


CorrectionBase::~CorrectionBase()
{}


void CorrectionBase::addKick(int plane, Row &row, double kick) {
    // Find object corresponding to current corrector.
    const string &name = row.getElement()->getName();
    Object *elem = OPAL.find(name);

    // Decide on type of corrector.
    if(dynamic_cast<OpalKicker *>(elem)) {
        int attrNumber = (plane == 0) ? OpalKicker::HKICK : OpalKicker::VKICK;
        kick += Attributes::getReal(elem->itsAttr[attrNumber]);
        Attributes::setReal(elem->itsAttr[attrNumber], kick);
    } else if(dynamic_cast<OpalHKicker *>(elem)) {
        kick += Attributes::getReal(elem->itsAttr[OpalHKicker::KICK]);
        Attributes::setReal(elem->itsAttr[OpalHKicker::KICK], kick);
    } else if(dynamic_cast<OpalVKicker *>(elem)) {
        kick += Attributes::getReal(elem->itsAttr[OpalVKicker::KICK]);
        Attributes::setReal(elem->itsAttr[OpalVKicker::KICK], kick);
    } else {
        throw OpalException("CorrectionBase::addKick()",
                            "Orbit Corrector \"" + name + "\" not found.");
    }

    // Update the CLASSIC corrector.
    elem->update();
}


void CorrectionBase::listCorrectors(bool list, int plane) {
    double rms = 0.0;
    double top = 0.0;
    std::cout << "\nTotal corrector settings for plane "
              << (plane == 0 ? 'X' : 'Y') << '\n';

    for(LocalIter iter = correctorTable[plane].begin();
        iter != correctorTable[plane].end(); ++iter) {
        ElementBase *design = (*iter)->getElement()->removeWrappers();
        Corrector *corr = dynamic_cast<Corrector *>(design);
        BDipoleField &field = corr->getField();
        double setting = (plane == 0) ? - field.getBy() : field.getBx();
        setting *= Physics::c / OPAL.getP0();

        if(list) {
            std::cout << std::setw(20) << corr->getName()
                      << std::setw(12) << setting << '\n';
        }

        rms += setting * setting;
        top = std::max(abs(setting), top);
    }

    rms = sqrt(rms / double(correctorTable[plane].size()));
    std::cout << "R.m.s. : " << std::setw(20) << rms << '\n'
              << "Maximum: " << std::setw(20) << top << endl;
}


void CorrectionBase::listMonitors(bool list, int plane) {
    double rms = 0.0;
    double top = 0.0;
    int q = 2 * plane;
    std::cout << "\nFinal monitor readings for plane "
              << (plane == 0 ? 'X' : 'Y') << '\n';

    for(LocalIter iter = monitorTable[plane].begin();
        iter != monitorTable[plane].end(); ++iter) {
        double reading = (*iter)->orbit[q];
        if(list) {
            std::cout << std::setw(20) << (*iter)->getElement()->getName()
                      << std::setw(12) << reading << '\n';
        }

        rms += reading * reading;
        top = std::max(abs(reading), top);
    }

    rms = sqrt(rms / double(monitorTable[plane].size()));
    std::cout << "R.m.s. : " << std::setw(20) << rms << '\n'
              << "Maximum: " << std::setw(20) << top << endl;
}


void CorrectionBase::setupTables() {
    IdealMapper mapper(itsLine, reference, false, false);
    double arc = 0.0;

    for(TLine::iterator iter = itsLine.begin();
        iter != itsLine.end(); ++iter) {
        // Build ideal matrix.
        iter->accept(mapper);
        mapper.getMatrix(iter->matrix);
        iter->isUsed[0] = iter->isUsed[1] = false;

        // Accumulate arc length.
        ElementBase *elem = iter->getElement();
        arc += elem->getElementLength();
        iter->arc = arc;

        // Put correctors and monitors in relevant table.
        test(elem);
        if(isCorr[0]) correctorTable[0].push_back(&*iter);
        if(isCorr[1]) correctorTable[1].push_back(&*iter);
        if(isMoni[0]) monitorTable[0].push_back(&*iter);
        if(isMoni[1]) monitorTable[1].push_back(&*iter);
    }
}


void CorrectionBase::test(ElementBase *base) {
    ElementBase *design = base->removeWrappers();
    isCorr[0] = isCorr[1] = isMoni[0] = isMoni[1] = false;

    if(Corrector *corr = dynamic_cast<Corrector *>(design)) {
        switch(corr->getPlane()) {

            case Corrector::X:
                isCorr[0] = true;
                break;

            case Corrector::Y:
                isCorr[1] = true;
                break;

            case Corrector::XY:
                isCorr[0] = isCorr[1] = true;
                break;

            default:
                break;
        }
    } else if(Monitor *moni = dynamic_cast<Monitor *>(design)) {
        switch(moni->getPlane()) {

            case Monitor::X:
                isMoni[0] = true;
                break;

            case Monitor::Y:
                isMoni[1] = true;
                break;

            case Monitor::XY:
                isMoni[0] = isMoni[1] = true;
                break;

            default:
                break;
        }
    }
}
