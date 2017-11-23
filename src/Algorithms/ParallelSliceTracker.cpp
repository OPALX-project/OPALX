#include "Algorithms/ParallelSliceTracker.h"

#include <memory>
#include <string>

#include "AbstractObjects/OpalData.h"
#include "Beamlines/Beamline.h"
#include "Distribution/Distribution.h"
#include "Lines/Sequence.h"
#include "Utilities/Timer.h"
#include "AbsBeamline/SpecificElementVisitor.h"

class PartData;

ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
        const PartData &reference, bool revBeam, bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack) {
    itsOpalBeamline_m = std::unique_ptr<OpalBeamline>(new OpalBeamline());
}


ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
        EnvelopeBunch &bunch, DataSink &ds, const PartData &reference,
        bool revBeam, bool revTrack, int maxSTEPS, double zstop):
    Tracker(beamline, reference, revBeam, revTrack),
    maxSteps_m(maxSTEPS),
    zstop_m(zstop) {

    itsBunch_m      = &bunch;
    itsDataSink_m   = &ds;

    itsOpalBeamline_m = std::unique_ptr<OpalBeamline>(new OpalBeamline());

    timeIntegrationTimer1_m  = IpplTimings::getTimer("Time integration1");
    timeIntegrationTimer2_m  = IpplTimings::getTimer("Time integration2");
    timeFieldEvaluation_m    = IpplTimings::getTimer("Field evaluation");
    BinRepartTimer_m         = IpplTimings::getTimer("Time of Binary repart.");
    WakeFieldTimer_m         = IpplTimings::getTimer("Time of Wake Field calc.");

    cavities_m = itsOpalBeamline_m->getElementByType(ElementBase::RFCAVITY);
    auto travelingwaves = itsOpalBeamline_m->getElementByType(ElementBase::TRAVELINGWAVE);
    cavities_m.merge(travelingwaves, ClassicField::SortAsc);
}


ParallelSliceTracker::~ParallelSliceTracker()
{}


/*
 * The maximum phase is added to the nominal phase of
 * the element. This is done on all nodes except node 0 where
 * the Autophase took place.
 */
void ParallelSliceTracker::updateRFElement(std::string elName, double maxPhase) {
    double phase = 0.0;

    for (FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++fit) {
        if ((*fit).getElement()->getName() == elName) {
            if ((*fit).getElement()->getType() == ElementBase::TRAVELINGWAVE) {
                phase  =  static_cast<TravelingWave *>((*fit).getElement().get())->getPhasem();
                static_cast<TravelingWave *>((*fit).getElement().get())->updatePhasem(phase + maxPhase);
            } else {
                phase  = static_cast<RFCavity *>((*fit).getElement().get())->getPhasem();
                static_cast<RFCavity *>((*fit).getElement().get())->updatePhasem(phase + maxPhase);
            }

            break;
        }
    }
}


/*
 * All RF-Elements gets updated, where the phiShift is the
 * global phase shift in units of seconds.
 */
void ParallelSliceTracker::printRFPhases() {
    Inform msg("ParallelSliceTracker ");

    FieldList &cl = cavities_m;

    const double RADDEG = 180.0 / Physics::pi;
    const double globalTimeShift = OpalData::getInstance()->getGlobalPhaseShift();

    msg << "\n-------------------------------------------------------------------------------------\n";

    for (FieldList::iterator it = cl.begin(); it != cl.end(); ++it) {
        std::shared_ptr<Component> element(it->getElement());
        std::string name = element->getName();
        double frequency;
        double phase;

        if (element->getType() == ElementBase::TRAVELINGWAVE) {
            phase = static_cast<TravelingWave *>(element.get())->getPhasem();
	    frequency = static_cast<TravelingWave *>(element.get())->getFrequencym();
        } else {
            phase = static_cast<RFCavity *>(element.get())->getPhasem();
	    frequency = static_cast<RFCavity *>(element.get())->getFrequencym();
        }

        msg << (it == cl.begin()? "": "\n")
            << name
            << ": phi = phi_nom + phi_maxE + global phase shift = " << (phase - globalTimeShift * frequency) * RADDEG << " degree, "
            << "(global phase shift = " << -globalTimeShift *frequency *RADDEG << " degree) \n";
    }

    msg << "-------------------------------------------------------------------------------------\n"
	<< endl;
}

void ParallelSliceTracker::applyEntranceFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {
}


void ParallelSliceTracker::applyExitFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {
}


void ParallelSliceTracker::execute() {

    Inform msg("ParallelSliceTracker");

    currentSimulationTime_m = itsBunch_m->getT();

    setTime();
    setLastStep();
    dumpPhaseSpaceOnScan();

    msg << "Track start at: t= " << itsBunch_m->getT()
        << " zstop@ " << zstop_m << " [m]" << endl;

    unsigned long long step = OpalData::getInstance()->getLastStep();
    msg << "executing ParallelSliceTracker, initial DT " << itsBunch_m->getdT()
        << " [s]; max integration steps " << maxSteps_m << " step= " << step << endl
        << "the mass is: " << itsReference.getM() * 1e-6 << " MeV, its charge: "
        << itsReference.getQ() << endl;

    prepareSections();
    handleAutoPhasing();

    Vector_t rmin, rmax;

    for (; step < maxSteps_m; ++step) {

        globalEOL_m = true;

        switchElements();
        computeExternalFields();

        //reduce(&globalEOL_m, &globalEOL_m, OpBitwiseOrAssign());
        //reduce(&globalEOL_m, &globalEOL_m + 1, &globalEOL_m, OpBitwiseAndAssign());
        MPI_Allreduce(MPI_IN_PLACE, &globalEOL_m, 1, MPI_INT, MPI_LAND, Ippl::getComm());

        computeSpaceChargeFields();
        timeIntegration();

        if (step % 100 == 0)
            msg << " Step " << step << " at " << itsBunch_m->zAvg() << " [m] t= "
                << itsBunch_m->getT() << " [s] E=" << itsBunch_m->Eavg() * 1e-6
                << " [MeV]" << endl;

        currentSimulationTime_m += itsBunch_m->getdT();
        itsBunch_m->setT(currentSimulationTime_m);

        dumpStats(step);

        if (hasEndOfLineReached()) break;
    }

    OpalData::getInstance()->setLastStep(step);
    writeLastStepPhaseSpace(step, itsBunch_m->get_sPos());
    itsOpalBeamline_m->switchElementsOff();
    msg << "done executing ParallelSliceTracker" << endl;
}


void ParallelSliceTracker::timeIntegration() {

    IpplTimings::startTimer(timeIntegrationTimer1_m);
    itsBunch_m->timeStep(itsBunch_m->getdT());
    IpplTimings::stopTimer(timeIntegrationTimer1_m);
}


void ParallelSliceTracker::handleAutoPhasing() {

    if (Options::autoPhase == 0) return;

    if(!OpalData::getInstance()->inRestartRun()) {
        itsDataSink_m->storeCavityInformation();
    }

    auto it = OpalData::getInstance()->getFirstMaxPhases();
    auto end = OpalData::getInstance()->getLastMaxPhases();
    for (; it != end; ++ it) {
        updateRFElement((*it).first, (*it).second);
    }
    printRFPhases();
}


void ParallelSliceTracker::prepareSections() {

    Inform msg("ParallelSliceTracker");
    itsBeamline_m.accept(*this);
    itsOpalBeamline_m->prepareSections();
    itsOpalBeamline_m->print(msg);

    for (int i = 0; i < itsBunch_m->getLocalNum(); i++) {
        Vector_t pos = Vector_t(itsBunch_m->getX(i),
                itsBunch_m->getY(i), itsBunch_m->getZ(i));
        itsOpalBeamline_m->getSectionIndexAt(pos, l);
    }
}


void ParallelSliceTracker::computeExternalFields() {

    IpplTimings::startTimer(timeFieldEvaluation_m);

    Vector_t externalE, externalB, KR, KT;
    bool globalEOL_m;

    for (int i = 0; i < itsBunch_m->getLocalNum(); i++) {

        externalB = Vector_t(0.0);
        externalE = Vector_t(0.0);
        KR        = Vector_t(0.0);
        KT        = Vector_t(0.0);

        //FIXME: why not x=y=0.0?
        Vector_t pos = Vector_t(itsBunch_m->getX(i),
                itsBunch_m->getY(i), itsBunch_m->getZ(i));

        unsigned long rtv = itsOpalBeamline_m->getFieldAt(i, pos, ls,
                currentSimulationTime_m , externalE, externalB);

        globalEOL_m = globalEOL_m && (rtv & BEAMLINE_EOL);

        itsOpalBeamline_m->getKFactors(i, pos, ls, currentSimulationTime_m, KR, KT);

        itsBunch_m->setExternalFields(i, externalE, externalB, KR, KT);
    }

    IpplTimings::stopTimer(timeFieldEvaluation_m);
}


void ParallelSliceTracker::computeSpaceChargeFields() {

    itsBunch_m->computeSpaceCharge();
}


void ParallelSliceTracker::dumpStats(long long step) {

    double sposRef = itsBunch_m->get_sPos();
    if (step != 0 && (step % Options::psDumpFreq == 0 || step % Options::statDumpFreq == 0))
        writePhaseSpace(step, sposRef);
}


void ParallelSliceTracker::switchElements(double scaleMargin) {

    double margin = 1.0;

    itsOpalBeamline_m->resetStatus();
    currentSimulationTime_m = itsBunch_m->getT();
    itsOpalBeamline_m->switchElements(itsBunch_m->zTail() - margin,
                                      itsBunch_m->zHead() + margin,
                                      itsBunch_m->Eavg() * 1e-6);
}


void ParallelSliceTracker::setLastStep() {

    unsigned long long step = 0;

    if (OpalData::getInstance()->inRestartRun()) {

        int prevDumpFreq = OpalData::getInstance()->getRestartDumpFreq();
        step = OpalData::getInstance()->getRestartStep() * prevDumpFreq + 1;
        maxSteps_m += step;
    } else if (OpalData::getInstance()->hasBunchAllocated() && Options::scan) {

        step = 1;
        if (!itsBunch_m->doEmission())
            writePhaseSpace(step - 1, 0.0); // write initial phase space
    } else {

        step = OpalData::getInstance()->getLastStep() + 1;
        maxSteps_m += step;
    }

    OpalData::getInstance()->setLastStep(step);
}


void ParallelSliceTracker::setTime() {

    if (OpalData::getInstance()->inRestartRun())
        currentSimulationTime_m = itsBunch_m->getT();
    else if (OpalData::getInstance()->hasBunchAllocated() && Options::scan)
        itsBunch_m->setT(0.0);
    else
        currentSimulationTime_m = itsBunch_m->getT();
}


void ParallelSliceTracker::dumpPhaseSpaceOnScan() {

    unsigned long long step = OpalData::getInstance()->getLastStep();

    if (OpalData::getInstance()->hasBunchAllocated() && Options::scan) {

        if (!itsBunch_m->doEmission())
            writePhaseSpace(step - 1, 0.0);
    }
}


bool ParallelSliceTracker::hasEndOfLineReached() {
    return (itsBunch_m->get_sPos() > zstop_m || !(itsBunch_m->isValid_m));
}
