#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include "Utility/IpplInfo.h"
#include "Algorithms/Vektor.h"

#include <fstream>
#include <string>
#include <vector>

class PeakFinder {
    
public:
    typedef std::vector<double> container_t;
    
public:
    
    PeakFinder();

    PeakFinder(std::string elem);
    
    void addParticle(const Vector_t& R);
    
    void createHistogram();
    
    void save();
    
    inline void setNumBins(unsigned int nBins);
    
    /** 
      * Find peaks of probe - function based on implementation in probe programs
      * @param[in] smoothingNumber   Smooth nr measurements
      * @param[in] minArea           Minimum Area for a single peak
      * @param[in] minFractionalArea Minimum fractional Area
      * @param[in] minAreaAboveNoise Minimum area above noise
      * @param[in] minSlope          Minimum slope
    */
    void findPeaks(int smoothingNumber, double minArea, double minFractionalArea, double minAreaAboveNoise, double minSlope);
    /** 
     * Analyse single peak
     * @param[in]  values     probe values
     * @param[in]  positions  probe positions
     * @param[in]  startIndex start position index
     * @param[in]  endIndex   end   position index
     * @param[out] peakRadius peak radius
     * @param[out] fourSigma  four sigma width
     */
    void analysePeak(const std::vector<double>& values,
		     const std::vector<double>& positions,
		     const int startIndex, const int endIndex,
		     float& peak,
		     float& fourSigma)const;
                     
private:
    
    void open_m();
    
    void append_m();
    
    void close_m();
    
    void saveASCII_m();
    
    
private:
     container_t radius_m;
     container_t globHist_m;
     
     // filename without extension
    std::string fn_m;

    // used to write out the data
    std::ofstream os_m;
    
    std::string element_m;
    
    unsigned int nBins_m;
    
    double binWidth_m;
    
    /// Radial position of peaks
    std::vector<double> peakRadii_m;
    /// Four sigma width of peaks
    std::vector<double> fourSigmaPeaks_m;
};

#endif