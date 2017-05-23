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
      * Analyse peaks of profile measurement and save in profile measurement
      *       * @param[in] finger Finger to analyse, default of -1 means default finger
      *          */
     void analysePeaks(int smoothingNumber,
                       double minArea,
                       double minFractionalArea,
                       double minAreaAboveNoise,
                       double minSlope);
     
     /// Analyse single peak
     void analysePeak(const std::vector<float>& values,
                      const std::vector<float>& positions, 
                      const int startIndex, const int endIndex,
                      float& peak,
                      float& rightPeak,
                      float& mean,
                      float& fourSigma,
                      std::pair<float,float>& radiiEnd,
                      std::pair<float,float>& radii4Perc,
                      std::pair<float,float>& radii25Perc) const;
    
    
            
private:
    std::vector<double> radius_m;
    std::vector<double> globHist_m;
    
    std::string element_m;
    int nBins_m;
    int binWidth_m;
    
};

#endif