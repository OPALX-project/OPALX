#ifndef OPAL_MULTI_BUNCH_HANDLER_H
#define OPAL_MULTI_BUNCH_HANDLER_H

#include "Algorithms/PartBunchBase.h"
#include "Structure/DataSink.h"

/* Helper class that stores bunch injection
 * information like azimuth, radius etc. of first
 * bunch in multi-bunch mode of ParallelCyclotronTracker.
 */
class MultiBunchHandler {
public:
    // multi-bunch modes
    enum class MB_MODE {
        FORCE  = 0,
        AUTO   = 1
    };

    // multi-bunch binning type
    enum class MB_BINNING {
        GAMMA = 0,
        BUNCH = 1
    };

    /*
     * @param numBunch > 1 --> multi bunch mode
     * @param eta binning value
     * @param para only for MB_MODE::AUTO
     * @param mode of multi-bunch
     * @param binning type of particle binning
     */
    MultiBunchHandler(PartBunchBase<double, 3> *beam,
                      const int& numBunch,
                      const double& eta,
                      const double& para,
                      const std::string& mode,
                      const std::string& binning,
                      DataSink& ds);

    void saveBunch(PartBunchBase<double, 3> *beam,
                   const double& azimuth);

    bool readBunch(PartBunchBase<double, 3> *beam,
                   const PartData& ref);

    /* Returns:
     *     0 - if nothing happened
     *     1 - if bunch got saved
     *     2 - if bunch got injected
     */
    short injectBunch(PartBunchBase<double, 3> *beam,
                      const PartData& ref,
                      bool& flagTransition,
                      const double& azimuth);

    void updateParticleBins(PartBunchBase<double, 3> *beam);
    
    /// set the working sub-mode for multi-bunch mode: "FORCE" or "AUTO"
    void setMode(const std::string& mbmode);

    // set binning type
    void setBinning(std::string binning);

    void setRadiusTurns(const double& radius);

    /// set total number of tracked bunches
    void setTotalNumBunch(short n);

    /// get total number of tracked bunches
    short getTotalNumBunch() const;

    void setNumBunch(short n);

    short getNumBunch() const;

    bool isForceMode() const;

    void setPathLength(const double& pathlength, short bunchNr);

    const double& getPathLength(short bunchNr) const;

    void setAzimuth(const double& azimuth, short bunchNr);

    const double& getAzimuth(short bunchNr) const;

    void setTime(const double& time, short bunchNr);

    const double& getTime(short bunchNr) const;

private:
    // store the data of the beam which are required for injecting a
    // new bunch for multibunch filename
    std::string onebunch_m;

    /// The number of bunches specified in TURNS of RUN commond
    short numBunch_m;

    // parameter for reset bin in multi-bunch run
    double eta_m;

    // 0 for single bunch (default),
    // 1 for FORCE,
    // 2 for AUTO
    MB_MODE mode_m; //multiBunchMode_m;

    // 0 for GAMMA (default),
    // 1 for BUNCH
    MB_BINNING binning_m; //binningType_m;

    // control parameter for AUTO multi-bunch mode
    double coeffDBunches_m;

    // used for automatic injection in multi-bunch mode
    double radiusLastTurn_m;
    double radiusThisTurn_m;

    DataSink& ds_m;

    // record how many bunches has already been injected.
    short bunchCount_m;
};


inline
void MultiBunchHandler::setTotalNumBunch(short n) {
    numBunch_m = n;
}


inline
short MultiBunchHandler::getTotalNumBunch() const {
    return numBunch_m;
}


inline
void MultiBunchHandler::setNumBunch(short n) {
    bunchCount_m = n;
}


inline
short MultiBunchHandler::getNumBunch() const {
    return bunchCount_m;
}


inline
bool MultiBunchHandler::isForceMode() const {
    return (mode_m == MB_MODE::FORCE);
}


inline
void MultiBunchHandler::setPathLength(const double& pathlength, short bunchNr) {
    MultiBunchDump *mbd = ds_m.getMultiBunchWriter(bunchNr);
    mbd->setPathLength(pathlength);
}


inline
const double& MultiBunchHandler::getPathLength(short bunchNr) const {
    MultiBunchDump *mbd = ds_m.getMultiBunchWriter(bunchNr);
    return mbd->getPathLength();
}


inline
void MultiBunchHandler::setAzimuth(const double& azimuth, short bunchNr) {
    MultiBunchDump *mbd = ds_m.getMultiBunchWriter(bunchNr);
    mbd->setAzimuth(azimuth);
}


inline
const double& MultiBunchHandler::getAzimuth(short bunchNr) const {
    MultiBunchDump *mbd = ds_m.getMultiBunchWriter(bunchNr);
    return mbd->getAzimuth();
}


inline
void MultiBunchHandler::setTime(const double& time, short bunchNr) {
    MultiBunchDump *mbd = ds_m.getMultiBunchWriter(bunchNr);
    mbd->setTime(time);
}


inline
const double& MultiBunchHandler::getTime(short bunchNr) const {
    MultiBunchDump *mbd = ds_m.getMultiBunchWriter(bunchNr);
    return mbd->getTime();
}

#endif
