// ------------------------------------------------------------------------
// $RCSfile: Monitor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Monitor
//   Defines the abstract interface for a beam position monitor.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------
#include "AbsBeamline/Monitor.h"
#include "Physics/Physics.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include <boost/filesystem.hpp>
#include "AbstractObjects/OpalData.h"
#include "Structure/MonitorStatisticsWriter.h"

#include <fstream>
#include <memory>


// Class Monitor
// ------------------------------------------------------------------------

std::map<double, SetStatistics> Monitor::statFileEntries_sm;
const double Monitor::halfLength_s = 0.005;

Monitor::Monitor():
    Monitor("")
{}


Monitor::Monitor(const Monitor &right):
    Component(right),
    filename_m(right.filename_m),
    plane_m(right.plane_m),
    type_m(right.type_m),
    numPassages_m(0)
{}


Monitor::Monitor(const std::string &name):
    Component(name),
    filename_m(""),
    plane_m(OFF),
    type_m(SPATIAL),
    numPassages_m(0)
{}


Monitor::~Monitor()
{}


void Monitor::accept(BeamlineVisitor &visitor) const {
    visitor.visitMonitor(*this);
}

bool Monitor::apply(const size_t &i, const double &t, Vector_t &/*E*/, Vector_t &/*B*/) {
    const Vector_t &R = RefPartBunch_m->R[i];
    const Vector_t &P = RefPartBunch_m->P[i];
    const double &dt = RefPartBunch_m->dt[i];
    const Vector_t singleStep  = Physics::c * dt * Util::getBeta(P);
    if (online_m && type_m == SPATIAL) {
        if (dt * R(2) < 0.0 &&
            dt * (R(2) + singleStep(2)) > 0.0) {
            double frac = R(2) / singleStep(2);

            lossDs_m->addParticle(R + frac * singleStep,
                                  P,
                                  RefPartBunch_m->ID[i],
                                  t + frac * dt,
                                  0);
        }
    }

    return false;
}

bool Monitor::applyToReferenceParticle(const Vector_t &R,
                                       const Vector_t &P,
                                       const double &t,
                                       Vector_t &,
                                       Vector_t &) {
    if (!OpalData::getInstance()->isInPrepState()) {
        const double dt = RefPartBunch_m->getdT();
        const double cdt = Physics::c * dt;
        const Vector_t singleStep = cdt * Util::getBeta(P);

        if (dt * R(2) < 0.0 &&
            dt * (R(2) + singleStep(2)) > 0.0) {
            double frac = -R(2) / singleStep(2);
            double time = t + frac * dt;
            Vector_t dR = frac * singleStep;
            double ds = euclidean_norm(dR + 0.5 * singleStep);
            lossDs_m->addReferenceParticle(csTrafoGlobal2Local_m.transformFrom(R + dR),
                                           csTrafoGlobal2Local_m.rotateFrom(P),
                                           time,
                                           RefPartBunch_m->get_sPos() + ds,
                                           RefPartBunch_m->getGlobalTrackStep());

            if (type_m == TEMPORAL) {
                const unsigned int localNum = RefPartBunch_m->getLocalNum();

                for (unsigned int i = 0; i < localNum; ++ i) {
                    Vector_t shift = ((frac - 0.5) * cdt * Util::getBeta(RefPartBunch_m->P[i])
                                      - singleStep);
                    lossDs_m->addParticle(RefPartBunch_m->R[i] + shift,
                                          RefPartBunch_m->P[i],
                                          RefPartBunch_m->ID[i],
                                          time,
                                          0);
                }
                OpalData::OPENMODE openMode;
                if (numPassages_m > 0) {
                    openMode = OpalData::OPENMODE::APPEND;
                } else {
                    openMode = OpalData::getInstance()->getOpenMode();
                }
                lossDs_m->save(1, openMode);
            }

            ++ numPassages_m;
        }
    }
    return false;
}

void Monitor::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    endField = startField + halfLength_s;
    startField -= halfLength_s;

    if (filename_m == std::string(""))
        filename_m = getName();
    else
        filename_m = filename_m.substr(0, filename_m.rfind("."));

    const size_t totalNum = bunch->getTotalNum();
    double currentPosition = endField;
    if (totalNum > 0) {
        currentPosition = bunch->get_sPos();
    }

    if (OpalData::getInstance()->getOpenMode() == OpalData::OPENMODE::WRITE ||
        currentPosition < startField) {
        namespace fs = boost::filesystem;

        fs::path lossFileName = fs::path(filename_m + ".h5");
        if (fs::exists(lossFileName)) {
            Ippl::Comm->barrier();
            if (Ippl::myNode() == 0)
                fs::remove(lossFileName);

            Ippl::Comm->barrier();
        }
    }

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m, !Options::asciidump, getType()));
}

void Monitor::finalise() {

}

void Monitor::goOnline(const double &) {
    online_m = true;
}

void Monitor::goOffline() {
    auto stats = lossDs_m->computeStatistics(numPassages_m);
    for (auto &stat: stats) {
        statFileEntries_sm.insert(std::make_pair(stat.spos_m, stat));
    }

    if (type_m != TEMPORAL) {
        lossDs_m->save(numPassages_m);
    }
}

bool Monitor::bends() const {
    return false;
}

void Monitor::setOutputFN(std::string fn) {
    filename_m = fn;
}

void Monitor::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = -halfLength_s;
    zEnd = halfLength_s;
}


ElementBase::ElementType Monitor::getType() const {
    return MONITOR;
}

void Monitor::writeStatistics() {
    if (statFileEntries_sm.size() == 0) return;

    std::string fileName = OpalData::getInstance()->getInputBasename() + std::string("_Monitors.stat");
    auto instance = OpalData::getInstance();
    bool hasPriorTrack = instance->hasPriorTrack();
    bool inRestartRun = instance->inRestartRun();

    auto it = statFileEntries_sm.begin();
    double spos = it->first;
    Util::rewindLinesSDDS(fileName, spos, false);

    MonitorStatisticsWriter writer(fileName, hasPriorTrack || inRestartRun);

    for (const auto &entry: statFileEntries_sm) {
        writer.addRow(entry.second);
    }

    statFileEntries_sm.clear();
}