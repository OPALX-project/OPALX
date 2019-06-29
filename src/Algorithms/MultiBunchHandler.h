#ifndef OPAL_MULTI_BUNCH_HANDLER_H
#define OPAL_MULTI_BUNCH_HANDLER_H

#include "Algorithms/PartBunchBase.h"

/* Helper class that stores bunch injection
 * information like azimuth, radius etc. of first
 * bunch in multi-bunch mode of ParallelCyclotronTracker.
 */
class MultiBunchHandler {
public:
    // multi-bunch modes
    enum class MB_MODE {
        NONE   = 0,
        FORCE  = 1,
        AUTO   = 2
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
                      const std::string& binning);

    void saveBunch(PartBunchBase<double, 3> *beam,
                   const double& azimuth);

    bool readBunch(PartBunchBase<double, 3> *beam,
                   const PartData& ref);

    /* Returns:
     *     0 - if nothing happened
     *     1 - if bunch got saved
     *     2 - if bunch got injected
     */
    int injectBunch(PartBunchBase<double, 3> *beam,
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
    void setNumBunch(int n);

    /// get total number of tracked bunches
    int getNumBunch() const;
    
    bool isMultiBunch() const;

    bool isForceMode() const;

    // record how many bunches has already been injected. ONLY FOR MPM
    int bunchCount_m;

private:
    // store the data of the beam which are required for injecting a
    // new bunch for multibunch filename
    std::string onebunch_m;

    /// The number of bunches specified in TURNS of RUN commond
    int numBunch_m;

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
};


inline
void MultiBunchHandler::setNumBunch(int n) {
    numBunch_m = n;
}


inline
int MultiBunchHandler::getNumBunch() const {
    return numBunch_m;
}


inline
bool MultiBunchHandler::isMultiBunch() const {
    return (mode_m != MB_MODE::NONE);
}


inline
bool MultiBunchHandler::isForceMode() const {
    return (mode_m == MB_MODE::FORCE);
}

#endif
