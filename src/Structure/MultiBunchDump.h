#ifndef OPAL_MULTI_BUNCH_DUMP_H
#define OPAL_MULTI_BUNCH_DUMP_H

#include "SDDSWriter.h"

class MultiBunchDump : public SDDSWriter {

public:
    struct beaminfo_t {
        double time;
        double azimuth;
        double pathlength;
        long unsigned int nParticles;
        double ekin;
        double dEkin;
        double rrms[3];
        double prms[3];
        double emit[3];
        double mean[3];
        double halo[3];
    };

public:
    MultiBunchDump(const std::string& fname, bool restart,
                   const short& bunch);

    void fillHeader();

    void write(PartBunchBase<double, 3>* beam, const double& azimuth);

    bool calcBunchBeamParameters(PartBunchBase<double, 3>* beam);

    /*!
     * @param azimuth in deg
     */
    void setAzimuth(const double& azimuth);

    /*!
     * @returns azimuth in deg
     */
    const double& getAzimuth() const;

    /*!
     * @param time in ns
     */
    void setTime(const double& time);

    /*!
     * @preturns time in ns
     */
    const double& getTime() const;

    /*!
     * @param pathlength in m
     */
    void setPathLength(const double& pathlength);

    /*!
     * @returns pathlength in m
     */
    const double& getPathLength() const;

private:
    bool isNotFirst_m;

    // the bunch number
    short bunch_m;

    beaminfo_t binfo_m;
};


inline
void MultiBunchDump::setAzimuth(const double& azimuth) {
    binfo_m.azimuth = azimuth;
}


inline
const double& MultiBunchDump::getAzimuth() const {
    return binfo_m.azimuth;
}


inline
void MultiBunchDump::setTime(const double& time) {
    binfo_m.time = time;
}


inline
const double& MultiBunchDump::getTime() const {
    return binfo_m.time;
}


inline
void MultiBunchDump::setPathLength(const double& pathlength) {
    binfo_m.pathlength = pathlength;
}


inline
const double& MultiBunchDump::getPathLength() const {
    return binfo_m.pathlength;
}

#endif
