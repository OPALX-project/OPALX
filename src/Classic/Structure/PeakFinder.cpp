#include "PeakFinder.h"

#include <algorithm>
#include <cmath>
#include <iterator>

#include "AbstractObjects/OpalData.h"
#include "Ippl.h"

PeakFinder::PeakFinder(std::string elem, double min, double max, double binWidth, bool singlemode):
    element_m(elem), binWidth_m(binWidth),
    min_m(min), max_m(max), singlemode_m(singlemode)
{
    if (min_m > max_m) {
        std::swap(min_m, max_m);
    }
    // calculate bins, round up so that histogram is large enough (add one for safety)
    nBins_m = static_cast<unsigned int>(std::ceil(( max_m - min_m ) / binWidth_m)) + 1;
}


void PeakFinder::addParticle(const Vector_t& R) {
    double radius = std::hypot(R(0),R(1));
    radius_m.push_back(radius);
}


void PeakFinder::save() {
    
    createHistogram_m();

    bool found = false;
    
    if ( Ippl::getNodes() == 1 && singlemode_m == true) {
        found = findPeaks();
    } else {
        found = findPeaks(smoothingNumber_m,
                          minArea_m,
                          minFractionalArea_m,
                          minAreaAboveNoise_m,
                          minSlope_m);
    }
    
    if ( found ) {
        fn_m   = element_m + std::string(".peaks");
        hist_m = element_m + std::string(".hist");
        
        INFOMSG("Save " << fn_m << " and " << hist_m << endl);
        
        if(OpalData::getInstance()->inRestartRun())
            this->append_m();
        else
            this->open_m();
        
        this->saveASCII_m();
        
        this->close_m();
        Ippl::Comm->barrier();
        
    }
    
    radius_m.clear();
    globHist_m.clear();
    peakRadii_m.clear();
    fourSigmaPeaks_m.clear();
}


bool PeakFinder::findPeaks() {
    for (const auto &radius : radius_m)
        peakRadii_m.push_back(radius);
    return !peakRadii_m.empty();
}


bool PeakFinder::findPeaks(int smoothingNumber,
                           double minAreaFactor,
                           double minFractionalAreaFactor,
                           double minAreaAboveNoise,
                           double minSlope)
{
    /* adapted from subroutine SEPAPR
     * Die Routine waehlt einen Beobachtungsindex.
     * Von diesem Aus wird fortlaufend die Peakflaeche
     * FTP integriert und mit dem aus dem letzten Messwert
     * und dessen Abstand vom Beobachtungspunkt gebildete
     * Dreieck ZPT verglichen. Ist FTP > ZPT ist ein neuer
     * Peak identifiziert. Der Beobachtungspunkt
     * verschiebt sich zum letzten Messwertindex und ab da
     * weiter, solange der Messwert abnimmt. Parallel wird
     * die Gesamtmessung aufintegriert, als 
     * zusaetzliches Kriterium zur Unterscheidung von echten
     * Peaks und Rauscheffekten.
     * 
     * smoothingNumber          Startindex in VAL für die Peakidentifikation
     * minAreaFactor            Zulässiger minimaler Anteil eines Einzelpeaks
     *                          am Messdatenintegral = Gewichtsfaktor für die
     *                          Elimination von Rauschpeaks
     * minFractionalAreaFactor  Gewichtsfaktor für die Gegenueberstellung FTP - ZPT
     *                          smoothen the data by summing neighbouring bins
     */
  
    const int size = static_cast<int>(globHist_m.size());
    if (size < smoothingNumber) {
        // no peaks can be found
        return false;
    }
    
    container_t smoothValues, sumSmoothValues;
    smoothValues.resize(size);
    sumSmoothValues.resize(size);
    double totalSum = 0.0;
    
    for (int i = smoothingNumber; i < size-smoothingNumber; i++) {
        double sum = 0.0;
        for (int j = -smoothingNumber; j <= smoothingNumber; j++) {
            sum += globHist_m[i+j];
        }
        sum /= smoothingNumber * 2 + 1;
        totalSum += sum;
        smoothValues[i] = sum;
        sumSmoothValues[i] = totalSum;
    }
    
    // set first and last values to value of smoothingNumber
    for (int i=0; i<smoothingNumber; i++) {
        smoothValues[i]        = smoothValues[smoothingNumber];
        smoothValues[size-1-i] = smoothValues[size-1-smoothingNumber];
    }
  
    std::vector<int> peakSeparatingIndices; // indices at minima (one more than number of peaks(!))
    peakSeparatingIndices.push_back(0);
  
    int nrPeaks            = 0;
    const double minArea   = minAreaFactor * totalSum; // minimum area for a peak
#ifdef PEAK_DEBUG
    INFOMSG("minArea " << minArea << endl);
#endif
    // number of indices corresponding to 10 mm
    const int maxIndex     = static_cast<int> (10 * size / binWidth_m);
    bool upwards           = false;
    bool newPeak           = false;
    
    for (int i = 1; i < size; i++) {
        int startIndex = std::max(i-maxIndex, peakSeparatingIndices.back());
        double ftp     = sumSmoothValues[i] - sumSmoothValues[startIndex];
        double ftpPeak = ftp - (i - startIndex)*smoothValues[startIndex]; // peak - noiselevel
        double slope   = (smoothValues[i] - smoothValues[startIndex]) / (i-startIndex);
        double zpt     = minFractionalAreaFactor * (smoothValues[i] - smoothValues[startIndex]) * (i - startIndex);
        
        if (ftpPeak >= zpt && ftp > minArea && ftpPeak > minAreaAboveNoise && slope > minSlope) {
#ifdef PEAK_DEBUG
            if (newPeak == false) {
                INFOMSG("Peak "     << peakSeparatingIndices.size() << endl);
                INFOMSG("Fraction " << ftpPeak << " " << zpt << endl);
                INFOMSG("Area "     << ftp     << " " << minArea << endl);
                INFOMSG("Noise "    << ftpPeak << " " << minAreaAboveNoise << endl);
                INFOMSG("Slope "    << slope   << " " << minSlope << endl);
            }
#endif
            newPeak = true;
        }
        
        if (smoothValues[i] > smoothValues [i-1] || i == size-1) {
            if (upwards == false || i == size-1) {
                upwards = true;
                if (newPeak == true) {
                    nrPeaks++;
                    peakSeparatingIndices.push_back(i-1);
                    newPeak = false;
                } else if (smoothValues[peakSeparatingIndices.back()] >= smoothValues[i]) {
                    peakSeparatingIndices.back() = i;
                }
            }
        } else {
            upwards = false;
        }
    }
    
    // debug
#ifdef PEAK_DEBUG
    INFOMSG("Number of peaks found: " << nrPeaks << endl);
#endif
    peakRadii_m.resize(nrPeaks);
    fourSigmaPeaks_m.resize(nrPeaks);
    
    // the position of the peak is bin centered
    container_t positions;
    positions.reserve(nBins_m);
    
    for (unsigned int i = 0; i < nBins_m; i++) {
        positions.push_back(min_m + (i + 0.5) * binWidth_m);
    }
    
    for (int i = 1; i < (int)(peakSeparatingIndices.size()); i++) {
        int startIndex = peakSeparatingIndices[i-1];
        int endIndex   = peakSeparatingIndices[i];
        analysePeak(globHist_m, positions, startIndex, endIndex,
                    peakRadii_m[i-1], fourSigmaPeaks_m[i-1]);
    }
    
    return !peakRadii_m.empty();
}


void PeakFinder::analysePeak(const container_t& values,
                             const container_t& positions,
                             const int startIndex, const int endIndex,
                             double& peak,
                             double& fourSigma)const
{
    // original subroutine ANALPR
    int range      = endIndex - startIndex;
    // find maximum
    double maximum   = -1;
    int maximumIndex = -1;
    int relMaxIndex  = -1;
    
    for (int j = startIndex; j <= endIndex; j++) {
        if (values[j] > maximum) {
            maximum = values[j];
            maximumIndex = j;
            relMaxIndex  = j - startIndex; // count from peak separation
        }
    }
    
    peak = positions[maximumIndex];
    
    // left limits, go down from peak to separation index
    int index20 = -1;
    int indexLeftEnd = 0; // left limit of peak
    
    for (int j = relMaxIndex; j >= 0; j--) {
        int index = j + startIndex;
        double value = values[index];
        if (value > 0.2 * maximum)
            index20 = j; // original code had i-1
        
        // if too far out, then break (not sure where formula comes from)
        if ( j < (3 * index20 - 2 * relMaxIndex) ) {
            indexLeftEnd = j;
            break;
        }
    }
    
    // right limits
    index20 = -1;
    int indexRightEnd = range; // right limit of peak
    
    // loop on right side of peak
    for (int j = relMaxIndex; j <= range; j++) {
        int index = j + startIndex;
        double value = values[index];
        if (value > 0.2 *maximum) {index20    = j;}
        // if too far out, then break (not sure where formula comes from)
        if ( j > (3 * index20 - 2 * relMaxIndex) ) {
            indexRightEnd = j;
            break;
        }
    }
    
    if (indexRightEnd - indexLeftEnd == 0) { // no peak
        fourSigma = 0.0;
        return; // return zeros for sigma
    }
    
    double sum=0.0, radialSum=0.0;
    
    for (int j = indexLeftEnd; j <= indexRightEnd; j++) {
        int index = j + startIndex;
        sum       += values[index];
        radialSum += values[index] * positions[index];
    }
    
    double mean = radialSum / sum;
    double variance = 0.0;
    
    for (int j = indexLeftEnd; j <= indexRightEnd; j++) {
        int index = j + startIndex;
        double value = values[index];
        double dx = positions[index] - mean;
        variance += value * dx * dx;
    }
    fourSigma = 4 * std::sqrt(variance / sum);
}


void PeakFinder::createHistogram_m() {
    /*
     * create local histograms
     */
    
    globHist_m.resize(nBins_m);
    container_t locHist(nBins_m,0.0);

    double invBinWidth = 1.0 / binWidth_m;
    for(container_t::iterator it = radius_m.begin(); it != radius_m.end(); ++it) {
        int bin = static_cast<int>(std::abs(*it - min_m ) * invBinWidth);
        if (bin < 0 || (unsigned int)bin >= nBins_m) continue; // Probe might save particles outside its boundary
        ++locHist[bin];
    }
    
    /*
     * create global histograms
     */
    reduce(&(locHist[0]), &(locHist[0]) + locHist.size(),
           &(globHist_m[0]), OpAddAssign());
}


void PeakFinder::open_m() {
    if ( Ippl::myNode() == 0 ) {
        os_m.open(fn_m.c_str(), std::ios::out);
        hos_m.open(hist_m.c_str(), std::ios::out);
    }
}


void PeakFinder::append_m() {
    if ( Ippl::myNode() == 0 ) {
        os_m.open(fn_m.c_str(), std::ios::app);
        hos_m.open(hist_m.c_str(), std::ios::app);
    }
}


void PeakFinder::close_m() {
    if ( Ippl::myNode() == 0 ) {
        os_m.close();
        hos_m.close();
    }
}


void PeakFinder::saveASCII_m() {
    if ( Ippl::myNode() == 0 )  {
        os_m << "# Peak Radii (mm)" << std::endl;
        for (auto &radius : peakRadii_m)
            os_m << radius << std::endl;
        
        hos_m << "# Histogram bin counts (min, max, nbins, binsize) "
              << min_m << " mm "
              << max_m << " mm "
              << nBins_m << " "
              << binWidth_m << " mm" << std::endl;
        for (auto binCount : globHist_m)
            hos_m << binCount << std::endl;
    }
}

