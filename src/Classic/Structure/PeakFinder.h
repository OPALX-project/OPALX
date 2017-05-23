#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include "Utility/IpplInfo.h"
#include "Algorithms/Vektor.h"

#include <fstream>
#include <string>
#include <vector>

class PeakFinder {
    
public:
    using container_t = std::vector<double>;
    
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

private:
    /***************
     * Output file *
     ***************/
    /// Open output file
    void open_m();
    /// Open output file in append mode
    void append_m();
    /// Close output file
    void close_m();
    /// Write to output file
    void saveASCII_m();
    
    /** 
     * Analyse single peak
     * @param[in]  values     probe values
     * @param[in]  positions  probe positions
     * @param[in]  startIndex start position index
     * @param[in]  endIndex   end   position index
     * @param[out] peakRadius peak radius
     * @param[out] fourSigma  four sigma width
     */
    void analysePeak(const container_t& values,
		     const container_t& positions,
		     const int startIndex, const int endIndex,
		     double & peak,
		     double& fourSigma)const;
                         
private:
    container_t radius_m;
    /// global histogram values
    container_t globHist_m;
     
    /// filename without extension
    std::string fn_m;

    // used to write out the data
    std::ofstream os_m;
    
    /// Element/probe name, for name output file
    std::string element_m;
    
    // Histogram details
    /// Number of bins
    unsigned int nBins_m;
    /// Bin width
    double binWidth_m;

    ///@{ Peak analysis parameters (copied from RRI2 probe program for now, need to be tuned a bit)
    const int    smoothingNumber_m   = 5;
    const double minArea_m           = 0.025;
    const double minFractionalArea_m = 0.6;
    const double minAreaAboveNoise_m = 5e-5;
    const double minSlope_m          = 1e-7;
    ///@}
    /// Radial position of peaks
    container_t peakRadii_m;
    /// Four sigma width of peaks
    container_t fourSigmaPeaks_m;
};

#endif