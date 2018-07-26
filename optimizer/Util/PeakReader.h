#ifndef __PEAKREADER_H__
#define __PEAKREADER_H__

#include <string>
#include <map>

#include "Util/OptPilotException.h"

/**
 *  \class PeakReader
 *  \brief Implements a parser and value extractor for peak files (*.peaks)
 */
class PeakReader {
    
public:
    
    PeakReader(std::string filename);
    ~PeakReader();
    
    void parseFile();
    
    /**
     * @param nPeak is the peak number
     * @param radius stores result [mm]
     */
    void getPeak(int nPeak, double& radius);
    
    /**
     * @returns the number of peaks in the file
     */
    std::size_t getNumberOfPeaks();
    
private:
    /// Peak filename
    std::string filename_m;
    
    /// all found peaks < peak number, radius >
    std::map<int, double> peaks_m;
};

#endif
