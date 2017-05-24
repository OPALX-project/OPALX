#include "PeakFinder.h"

#include <algorithm>
#include <cmath>
#include <iterator>

#include "AbstractObjects/OpalData.h"
#include "Ippl.h"

PeakFinder::PeakFinder(std::string elem):
    element_m(elem), nBins_m(10), binWidth_m(1),
    globMin_m(0.0),globMax_m(5000.0)
{}


PeakFinder::PeakFinder() : PeakFinder(std::string("NULL"))
{ }


void PeakFinder::addParticle(const Vector_t& R) {
    double radius = std::sqrt( dot(R, R) );
    radius_m.push_back(radius);
}


void PeakFinder::save() {
    fn_m = element_m + std::string(".peaks");
    
    createHistogram_m();

    findPeaks(smoothingNumber_m,
	      minArea_m,
	      minFractionalArea_m,
	      minAreaAboveNoise_m,
	      minSlope_m);

    INFOMSG("Save " << fn_m << endl);
    
    if(OpalData::getInstance()->inRestartRun())
	this->append_m();
    else {
        this->open_m();
    }
    this->saveASCII_m();
    
    this->close_m();
    Ippl::Comm->barrier();
    
    radius_m.clear();
    globHist_m.clear();
}


void PeakFinder::findPeaks(int smoothingNumber,
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
    container_t& values = globHist_m;
  
    const int size = static_cast<int>(values.size());
    if (size < smoothingNumber) {
	// no peaks can be found
	return;
    }
    container_t smoothValues;
    container_t sumSmoothValues;
    smoothValues.resize(size);
    sumSmoothValues.resize(size);
    double totalSum = 0.0;
    for (int i = smoothingNumber; i<size-smoothingNumber; i++) {
	double sum = 0.0;
	for (int j = -smoothingNumber; j<=smoothingNumber; j++) {
	    sum += values[i+j];
	}
	sum /= smoothingNumber*2+1;
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
    INFOMSG("minArea " << minArea << endl);
    // number of indices corresponding to 10 mm
    const int maxIndex     = static_cast<int> (10 * size / binWidth_m);
    bool upwards           = false;
    bool newPeak           = false;
    for (int i=1; i<size; i++) {
	int startIndex = std::max(i-maxIndex, peakSeparatingIndices.back());
	double ftp     = sumSmoothValues[i] - sumSmoothValues[startIndex];
	double ftpPeak = ftp - (i - startIndex)*smoothValues[startIndex]; // peak - noiselevel
	double slope   = (smoothValues[i] - smoothValues[startIndex]) / (i-startIndex);
	double zpt     = minFractionalAreaFactor * (smoothValues[i] - smoothValues[startIndex]) * (i - startIndex);
	if (ftpPeak >= zpt && ftp > minArea && ftpPeak > minAreaAboveNoise && slope > minSlope) {
	    if (newPeak == false) {
		INFOMSG("Peak "     << peakSeparatingIndices.size() << endl);
		// 	INFOMSG("Position " << histogram->getPosition(i) << endl);
		INFOMSG("Fraction " << ftpPeak << " " << zpt << endl);
		INFOMSG("Area "     << ftp     << " " << minArea << endl);
		INFOMSG("Noise "    << ftpPeak << " " << minAreaAboveNoise << endl);
		INFOMSG("Slope "    << slope   << " " << minSlope << endl);
	    }
	    newPeak = true;
	}
	if (smoothValues[i] > smoothValues [i-1] || i == size-1) {
	    if (upwards == false || i == size-1) {
		upwards = true;
		if (newPeak == true) {
		    nrPeaks++;
		    // 	  INFOMSG("Separating position " << histogram->getPosition(i) << endl);
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
    INFOMSG("Number of peaks found: " << nrPeaks << endl);
    peakRadii_m.resize(nrPeaks);
    fourSigmaPeaks_m.resize(nrPeaks);
    
    container_t positions;
    positions.reserve(nBins_m);
    for (unsigned int i=0; i<nBins_m; i++) {
	positions.push_back(globMin_m + (i+0.5)*binWidth_m);
    }
    
    for (int i=1; i<(int)(peakSeparatingIndices.size()); i++) {
	int startIndex = peakSeparatingIndices[i-1];
	int endIndex   = peakSeparatingIndices[i];
	analysePeak(values,positions,startIndex,endIndex,
		    peakRadii_m[i-1],fourSigmaPeaks_m[i-1]);
    }
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
    for (int j=startIndex; j<=endIndex; j++) {
	if (values[j] > maximum) {
	    maximum = values[j];
	    maximumIndex = j;
	    relMaxIndex  = j - startIndex; // count from peak separation
	}
    }
    peak = positions[maximumIndex];
    // qDebug() << "Peak " << i << " at " << positions[maximumIndex] << " mm";
    
    // left limits, go down from peak to separation index
    int index20 = -1;
    int indexLeftEnd = 0; // left limit of peak
    for (int j=relMaxIndex; j>=0; j--) {
	int index = j + startIndex;
	double value = values[index];
	if (value > 0.2 *maximum) {index20 = j;} // original code had i-1
	// if too far out, then break (not sure where formula comes from)
	if (j < (3*index20 - 2*relMaxIndex)) {
	    indexLeftEnd   = j;
	    break;
	}
    }
    // right limits
    index20 = -1;
    int indexRightEnd = range; // right limit of peak
    // loop on right side of peak
    for (int j=relMaxIndex; j<=range; j++) {
	int index = j + startIndex;
	double value = values[index];
	if (value > 0.2 *maximum) {index20    = j;}
	// if too far out, then break (not sure where formula comes from)
	if (j > (3*index20 - 2*relMaxIndex)) {
	    indexRightEnd   = j;
	    break;
	}
    }
    // INFOMSG("width of Peak " << indexRightEnd - indexLeftEnd << "steps " << endl);
    
    if (indexRightEnd - indexLeftEnd == 0) { // no peak
	fourSigma = 0.0;
	return; // return zeros for sigma
    }
    double sum=0.0, radialSum=0.0;
    for (int j=indexLeftEnd; j<=indexRightEnd; j++) {
	int index = j + startIndex;
	sum       += values[index];
	radialSum += values[index] * positions[index];
    }
    double mean = radialSum / sum;
    double variance = 0.0;
    for (int j=indexLeftEnd; j<=indexRightEnd; j++) {
	int index = j + startIndex;
	double value = values[index];
	double dx = positions[index] - mean;
	variance += value * dx * dx;
    }
    fourSigma = 4 * std::sqrt(variance / sum);
}


void PeakFinder::createHistogram_m() {

    double locMin=1e10, locMax=-1e10;
    if (!radius_m.empty()) {
	// compute global minimum and maximum radius
	auto result = std::minmax_element(radius_m.begin(), radius_m.end());
    
	locMin = *result.first;
	locMax = *result.second;
    }
    MPI_Allreduce(&locMin, &globMin_m, 1, MPI_DOUBLE, MPI_MIN, Ippl::getComm());
    MPI_Allreduce(&locMax, &globMax_m, 1, MPI_DOUBLE, MPI_MAX, Ippl::getComm());
    
    /*
     * create local histograms
     */
    
    binWidth_m = 1.0; // mm

    if (globMax_m < -1e9) nBins_m = 10; // no particles in probe
    // calculate bins, round up so that histogram is large enough
    else {nBins_m = static_cast<unsigned int>(std::ceil( globMax_m - globMin_m ) / binWidth_m);}

    // std::cout << "number of bins:      " << nBins_m << std::endl;
    // std::cout << "number of particles: " << radius_m.size() << std::endl;
    // std::cout << "global min , max     " << globMin_m << " " << globMax_m << std::endl;
    globHist_m.resize(nBins_m);
    
    container_t locHist(nBins_m);

    double invBinWidth = 1.0 / binWidth_m;
    for(container_t::iterator it = radius_m.begin(); it != radius_m.end(); ++it) {
        int bin = (*it - globMin_m ) * invBinWidth;
        ++locHist[bin];
    }
    
    /*
     * create global histograms
     */
    // reduce(&(locHist[0]), &(locHist[locHist.size()-1]),
    //        &(globHist_m[0]), OpAddAssign());
    MPI_Reduce(&locHist[0], &globHist_m[0], locHist.size(), MPI_DOUBLE, MPI_SUM, 0, Ippl::getComm());

    if (Ippl::myNode() == 0) {
        std::string histFilename = element_m + ".hist";
        std::ofstream os(histFilename);
        if (!os) {std::cout << "cant open histogram file" << std::endl; return;}
        for (auto value : globHist_m) {
            os << value << std::endl;
        }
        os.close();
    }
}


void PeakFinder::open_m() {
    if ( Ippl::myNode() == 0 ) {
        os_m.open(fn_m.c_str(), std::ios::out);
    }
}


void PeakFinder::append_m() {
    if ( Ippl::myNode() == 0 ) {
        os_m.open(fn_m.c_str(), std::ios::app);
    }
}


void PeakFinder::close_m() {
    if ( Ippl::myNode() == 0 )
        os_m.close();
}


void PeakFinder::saveASCII_m() {
    if ( Ippl::myNode() == 0 )  {
	os_m << "#Peak Radii " << std::endl;
        for (auto &radius : peakRadii_m) {
	    os_m << radius << std::endl;
	}
    }
}

