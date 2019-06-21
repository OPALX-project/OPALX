#ifndef OPAL_STAT_WRITER_H
#define OPAL_STAT_WRITER_H

#include "SDDSWriter.h"

class StatWriter : public SDDSWriter {

public:
    /** \brief
     *  delete the last 'numberOfLines' lines of the statistics file
     */
    unsigned int rewindToSpos(double maxSpos) const;
    
private:
    void writeHeader();
    
    
    /** \brief Write statistical data.
     *
     * Writes statistical beam data to proper output file. This is information such as RMS beam parameters
     * etc.
     *
     * Also gathers and writes load balancing data to load balance statistics file.
     * \param beam The beam.
     * \param FDext The external E and B field for the head, reference and tail particles. The vector array
     * has the following layout:
     *  - FDext[0] = B at head particle location (in x, y and z).
     *  - FDext[1] = E at head particle location (in x, y and z).
     *  - FDext[2] = B at reference particle location (in x, y and z).
     *  - FDext[3] = E at reference particle location (in x, y and z).
     *  - FDext[4] = B at tail particle location (in x, y, and z).
     *  - FDext[5] = E at tail particle location (in x, y, and z).
     */
    void writeData(PartBunchBase<double, 3> *beam,
                   PartBunchBase<double, 3> *beam, Vector_t FDext[],
                   const losses_t &losses, const double& azimuth = -1);
    
    
//     void writeStatData(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail);
    
private:
    /// Timer to track statistics write time.
    IpplTimings::TimerRef StatMarkerTimer_m;
};

#endif
