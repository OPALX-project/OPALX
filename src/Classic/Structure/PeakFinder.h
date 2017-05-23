#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include "Utility/IpplInfo.h"

class PeakFinder {
    
public:
    
    PeakFinder();

    PeakFinder(std::string elem);
    
    void addParticle(const Vector_t& R);
    
    void createHistogram();
    
    void save();
    
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
    void analysePeak(const std::vector<float>& values,
		     const std::vector<float>& positions, 
		     const int startIndex, const int endIndex,
		     float& peak,
		     float& fourSigma)const;
            
private:
    std::vector<double> radius_m;
    std::vector<double> globHist_m;
    
    std::string element_m;
    int nBins_m;
    int binWidth_m;
    
};

#endif