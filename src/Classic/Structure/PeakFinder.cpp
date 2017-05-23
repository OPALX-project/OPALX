#include "PeakFinder.h"

#include <algorithm>
#include <cmath>
#include <iterator>

extern Inform *gmsg;

PeakFinder::PeakFinder(std::string elem):
    element_m(elem), nBins_m(0), binWidth_m(0) {
}

PeakFinder::PeakFinder() {
    PeakFinder(std::string("NULL"));
}


void PeakFinder::addParticle(const Vector_t& R) {
    double radius = std::sqrt( R(0) * R(0) + R(2) * R(2) );
    radius_m.push_back(radius);
}

void PeakFinder::createHistogram() {
    
    // compute global minimum and maximum radius
    auto result = std::minmax_element(radius_m.begin(), radius_m.end());
    
    double locMin = *result.first;
    double locMax = *result.second;
    
    double globMin = 0;
    double globMax = 0;
    
    reduce(locMin, globMin, OpAddAssign());
    
    reduce(locMax, globMax, OpAddAssign());
    
    
    /*
     * create local histograms
     */
    
    std::vector<double> locHist(nBins_m);
    double invBinWidth = 1.0 / binWidth_m;
    for(std::vector<double>::iterator it = radius_m.begin(); it != radius_m.end(); ++it) {
        int bin = (*it - globMin ) * invBinWidth;
        ++locHist[bin];
    }
    
    /*
     * create global histograms
     */
    globHist_m.resize(nBins_m);
    reduce(&(locHist[0]), &(locHist[0]) + locHist.size(),
               &(globHist[0]), OpAddAssign());
    
    
    radius_m.clear();
}


void PeakFinder::save() {
    
    
}

void PeakFinder::findPeaks(int smoothingNumber, double minAreaFactor, double minFractionalAreaFactor, double minAreaAboveNoise, double minSlope)
{
  // adapted from subroutine SEPAPR
  // Die Routine waehlt einen Beobachtungsindex. Von diesem Aus wird fortlaufend die Peakflaeche FTP integriert und mit dem aus dem letzten Messwert
  // und dessen Abstand vom Beobachtungspunkt gebildete Dreieck ZPT verglichen. Ist FTP > ZPT ist ein neuer Peak identifiziert. Der Beobachtungspunkt
  // verschiebt sich zum letzten Messwertindex und ab da weiter, solange der Messwert abnimmt. Parallel wird die Gesamtmessung aufintegriert, als 
  // zusaetzliches Kriterium zur Unterscheidung von echten Peaks und Rauscheffekten.
  // smoothingNumber            Startindex in VAL für die Peakidentifikation
  // minAreaFactor              Zulässiger minimaler Anteil eines Einzelpeaks am Messdatenintegral = Gewichtsfaktor für die Elimination von Rauschpeaks
  // minFractionalAreaFactor    Gewichtsfaktor für die Gegenueberstellung FTP - ZPT
  // smoothen the data by summing neighbouring bins
  std::vector<float> values = histogram->values();
  
  const int size = static_cast<int>(values.size());
  if (size < smoothingNumber) {
    // no peaks can be found
    return;
  }
  std::vector<double> smoothValues;
  std::vector<double> sumSmoothValues;
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
  *gmsg << "minArea " << minArea << endl;
  // number of indices corresponding to 10 mm
  const int maxIndex     = static_cast<int> (10 * size / (histogram->width()));
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
	*gmsg << "Peak "     << peakSeparatingIndices.size() << endl;
	*gmsg << "Position " << histogram->getPosition(i) << endl;
        *gmsg << "Fraction " << ftpPeak << " " << zpt << endl;
	*gmsg << "Area "     << ftp     << " " << minArea << endl;
	*gmsg << "Noise "    << ftpPeak << " " << minAreaAboveNoise << endl;
	*gmsg << "Slope "    << slope   << " " << minSlope << endl;
      }
      newPeak = true;
    }
    if (smoothValues[i] > smoothValues [i-1] || i == size-1) {
      if (upwards == false || i == size-1) {
	upwards = true;
	if (newPeak == true) {
	  nrPeaks++;
	  *gmsg << "Separating position " << histogram->getPosition(i) << endl;
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
  *gmsg << "Number of peaks found: " << nrPeaks << endl;
  peakRadii.resize(nrPeaks);
  fourSigmaPeaks.resize(nrPeaks);
  const std::vector<float>& positions = histogram->getPositions();
  for (int i=1; i<(int)(peakSeparatingIndices.size()); i++) {
    int startIndex = peakSeparatingIndices[i-1];
    int endIndex   = peakSeparatingIndices[i];
    analysePeak(values,positions,startIndex,endIndex,
		peakRadii[i-1],fourSigmaPeaks[i-1]);
  }
}


void PeakFinder::analysePeak(const std::vector<float>& values,
			     const std::vector<float>& positions, 
			     const int startIndex, const int endIndex,
			     float& peak,
			     float& fourSigma)const
{
  // original subroutine ANALPR
  int range      = endIndex - startIndex;
  // find maximum
  float maximum   = -1;
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
    float value = values[index];
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
    float value = values[index];
    if (value > 0.2 *maximum) {index20    = j;}
    // if too far out, then break (not sure where formula comes from)
    if (j > (3*index20 - 2*relMaxIndex)) {
      indexRightEnd   = j;
      break;
    }
  }
  // *gmsg << "width of Peak " << indexRightEnd - indexLeftEnd << "steps " << endl;
    
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
    float value = values[index];
    double dx = positions[index] - mean;
    variance += value * dx * dx;
  }
  fourSigma = 4 * std::sqrt(variance / sum);
}
