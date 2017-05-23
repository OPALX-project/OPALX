// ------------------------------------------------------------------------
// $RCSfile: Probe.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Probe
//   Defines the abstract interface for a septum magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2009/10/07 09:32:32 $
// $Author: Bi, Yang $
// 2012/03/01: fix bugs and change the algorithm in the checkProbe()
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Probe.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include <iostream>
#include <fstream>
using Physics::pi;

extern Inform *gmsg;

// Class Probe
// ------------------------------------------------------------------------

Probe::Probe():
    Component(),
    filename_m(""),
    position_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0),
    step_m(0){
    A_m = yend_m - ystart_m;
    B_m = xstart_m - xend_m;
    R_m = sqrt(A_m*A_m+B_m*B_m);
    C_m = ystart_m*xend_m - xstart_m*yend_m;
}


Probe::Probe(const Probe &right):
    Component(right),
    filename_m(right.filename_m),
    position_m(right.position_m),
    xstart_m(right.xstart_m),
    xend_m(right.xend_m),
    ystart_m(right.ystart_m),
    yend_m(right.yend_m),
    width_m(right.width_m),
    step_m(right.step_m){
    A_m = yend_m - ystart_m;
    B_m = xstart_m - xend_m;
    R_m = sqrt(A_m*A_m+B_m*B_m);
    C_m = ystart_m*xend_m - xstart_m*yend_m;
}


Probe::Probe(const std::string &name):
    Component(name),
    filename_m(""),
    position_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0),
    step_m(0){
    A_m = yend_m - ystart_m;
    B_m = xstart_m - xend_m;
    R_m = sqrt(A_m*A_m+B_m*B_m);
    C_m = ystart_m*xend_m - xstart_m*yend_m;
}


Probe::~Probe() {
    idrec_m.clear();
}


void Probe::accept(BeamlineVisitor &visitor) const {
    visitor.visitProbe(*this);
}

bool Probe::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    return apply(i, t, Ev, Bv);
}

bool Probe::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool Probe::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

void Probe::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    position_m = startField;
    startField -= 0.005;
    endField = position_m + 0.005;
}

void Probe::initialise(PartBunch *bunch, const double &scaleFactor) {
    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
}

void Probe::finalise() {
    lossDs_m->save();
    *gmsg << "* Finalize probe " << getName() << endl;
}

bool Probe::bends() const {
    return false;
}

void Probe::goOffline() {
    online_m = false;
    lossDs_m->save();
}

void  Probe::setXstart(double xstart) {
    xstart_m = xstart;
}

void  Probe::setXend(double xend) {
    xend_m = xend;
}

void  Probe::setYstart(double ystart) {
    ystart_m = ystart;
}


void  Probe::setYend(double yend) {
    yend_m = yend;
}
void  Probe::setWidth(double width) {
    width_m = width;
}


double  Probe::getXstart() const {
    return xstart_m;
}

double  Probe::getXend() const {
    return xend_m;
}

double  Probe::getYstart() const {
    return ystart_m;
}

double  Probe::getYend() const {
    return yend_m;
}
double  Probe::getWidth() const {
    return width_m;
}

void Probe::setGeom(const double dist) {

    double slope;
    if (xend_m == xstart_m)
      slope = 1.0e12;
    else
      slope = (yend_m - ystart_m) / (xend_m - xstart_m);

    double coeff2 = sqrt(1 + slope * slope);
    double coeff1 = slope / coeff2;
    double halfdist = dist / 2.0;
    geom_m[0].x = xstart_m - halfdist * coeff1;
    geom_m[0].y = ystart_m + halfdist / coeff2;

    geom_m[1].x = xstart_m + halfdist * coeff1;
    geom_m[1].y = ystart_m - halfdist / coeff2;

    geom_m[2].x = xend_m + halfdist * coeff1;
    geom_m[2].y = yend_m - halfdist  / coeff2;

    geom_m[3].x = xend_m - halfdist * coeff1;
    geom_m[3].y = yend_m + halfdist / coeff2;

    geom_m[4].x = geom_m[0].x;
    geom_m[4].y = geom_m[0].y;
}

bool  Probe::checkProbe(PartBunch &bunch, const int turnnumber, const double t, const double tstep) {

    bool flagprobed = false;
    Vector_t rmin, rmax, probepoint;
    bunch.get_bounds(rmin, rmax);
    double r_start = sqrt(xstart_m * xstart_m + ystart_m * ystart_m);
    double r_end = sqrt(xend_m * xend_m + yend_m * yend_m);
    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));
    double r2 = sqrt(rmin(0) * rmin(0) + rmin(1) * rmin(1));

    if( r1 > r_start - 10.0 && r2 < r_end + 10.0 ){
        size_t tempnum = bunch.getLocalNum();
        int pflag = 0;

	Vector_t meanP(0.0, 0.0, 0.0);
	for(unsigned int i = 0; i < bunch.getLocalNum(); ++i) {
	  for(int d = 0; d < 3; ++d) {
            meanP(d) += bunch.P[i](d);
	  }
	}
	reduce(meanP, meanP, OpAddAssign());
	meanP = meanP / Vector_t(bunch.getTotalNum());

	double sk1, sk2, stangle = 0.0;
	if ( B_m == 0.0 ){
	  sk1 = meanP(1)/meanP(0);
	  if (sk1 == 0.0)
	    stangle = 1.0e12;
	  else
	    stangle = std::abs(1/sk1);
	}else if (meanP(0) == 0.0 ){
	  sk2 = - A_m/B_m;
	  if ( sk2 == 0.0 )
	    stangle = 1.0e12;
	  else
	    stangle = std::abs(1/sk2);
	}else {
	  sk1 = meanP(1)/meanP(0);
	  sk2 = - A_m/B_m;
	  stangle = std::abs(( sk1-sk2 )/(1 + sk1*sk2));
	}
	double lstep = (sqrt(1.0-1.0/(1.0+dot(meanP, meanP))) * Physics::c) * tstep*1.0e-6; // [mm]
	double Swidth = lstep / sqrt( 1 + 1/stangle/stangle );
	setGeom(Swidth);

    for(unsigned int i = 0; i < tempnum; ++i) {
	  pflag = checkPoint(bunch.R[i](0), bunch.R[i](1));
	  if(pflag != 0) {
	     // dist1 > 0, right hand, dt > 0; dist1 < 0, left hand, dt < 0
		  double dist1 = (A_m*bunch.R[i](0)+B_m*bunch.R[i](1)+C_m)/R_m/1000.0;
		  double k1, k2, tangle = 0.0;
		  if ( B_m == 0.0 ){
		    k1 = bunch.P[i](1)/bunch.P[i](0);
		    if (k1 == 0.0)
		      tangle = 1.0e12;
		    else
                      tangle = std::abs(1/k1);
		  }else if (bunch.P[i](0) == 0.0 ){
		    k2 = -A_m/B_m;
		    if (k2 == 0.0)
		      tangle = 1.0e12;
		    else
                      tangle = std::abs(1/k2);
		  }else {
		    k1 = bunch.P[i](1)/bunch.P[i](0);
		    k2 = -A_m/B_m;
		    tangle = std::abs(( k1-k2 )/(1 + k1*k2));
		  }
		  double dist2 = dist1 * sqrt( 1+1/tangle/tangle );
		  double dt = dist2/(sqrt(1.0-1.0/(1.0 + dot(bunch.P[i], bunch.P[i]))) * Physics::c)*1.0e9;

	    probepoint(0) = (B_m*B_m*bunch.R[i](0) - A_m*B_m*bunch.R[i](1)-A_m*C_m)/(R_m*R_m);
	    probepoint(1) = (A_m*A_m*bunch.R[i](1) - A_m*B_m*bunch.R[i](0)-B_m*C_m)/(R_m*R_m);
	    probepoint(2) = bunch.R[i](2);
	    lossDs_m->addParticle(probepoint, bunch.P[i], bunch.ID[i], t+dt, turnnumber);
	    flagprobed = true;
	  }
	}
    }

    reduce(&flagprobed, &flagprobed + 1, &flagprobed, OpBitwiseOrAssign());
    return flagprobed;
}

int Probe::checkPoint(const double &x, const double &y) {
    int    cn = 0;

    for(int i = 0; i < 4; i++) {
        if(((geom_m[i].y <= y) && (geom_m[i+1].y > y))
           || ((geom_m[i].y > y) && (geom_m[i+1].y <= y))) {

            float vt = (float)(y - geom_m[i].y) / (geom_m[i+1].y - geom_m[i].y);
            if(x < geom_m[i].x + vt * (geom_m[i+1].x - geom_m[i].x))
                ++cn;
        }
    }
    return (cn & 1);  // 0 if even (out), and 1 if odd (in)

}

void Probe::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m - 0.005;
    zEnd = position_m + 0.005;
}

ElementBase::ElementType Probe::getType() const {
    return PROBE;
}

void Probe::analysePeaks(int smoothingNumber, double minAreaFactor, double minFractionalAreaFactor, double minAreaAboveNoise, double minSlope)
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
  std::cout << "minArea " << minArea << std::endl;;
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
	std::cout << "Peak "     << peakSeparatingIndices.size();
	std::cout << "Position " << histogram->getPosition(i);
        std::cout << "Fraction " << ftpPeak << zpt;
	std::cout << "Area "     << ftp     << minArea;
	std::cout << "Noise "    << ftpPeak << minAreaAboveNoise;
	std::cout << "Slope "    << slope   << minSlope;
      }
      newPeak = true;
    }
    if (smoothValues[i] > smoothValues [i-1] || i == size-1) {
      if (upwards == false || i == size-1) {
	upwards = true;
	if (newPeak == true) {
	  nrPeaks++;
	  std::cout << "Separating position " << histogram->getPosition(i) << std::endl;
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
  std::cout << "Number of peaks found: " << nrPeaks << std::endl;;
  // get turn number and mid peak radius for display
//  QVector<float> peakRadii(nrPeaks), rightPeakRadii(nrPeaks), peakMean(nrPeaks), peakFourSigma(nrPeaks);
  //QVector<QPair<float,float>> radiiEnd(nrPeaks),radii4Perc(nrPeaks),radii25Perc(nrPeaks);
  const std::vector<float>& positions = histogram->getPositions();
  for (int i=1; i<(int)(peakSeparatingIndices.size()); i++) {
    int startIndex = peakSeparatingIndices[i-1];
    int endIndex   = peakSeparatingIndices[i];
    analysePeak(values,positions,startIndex,endIndex,
		peakRadii[i-1],rightPeakRadii[i-1],peakMean[i-1],peakFourSigma[i-1],
                radiiEnd[i-1],radii4Perc[i-1],radii25Perc[i-1]);
  }
}

void Probe::analysePeak(const std::vector<float>& values,
			     const std::vector<float>& positions, 
			     const int startIndex, const int endIndex,
			     float& peak,
			     float& rightPeak,
			     float& mean,
			     float& fourSigma,
			     std::pair<float,float>& radiiEnd,
			     std::pair<float,float>& radii4Perc,
			     std::pair<float,float>& radii25Perc)const
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
    if (value > 0.25*maximum) {radii25Perc.first = positions[index];}
    if (value > 0.04*maximum) { radii4Perc.first = positions[index];}
    // if too far out, then break (not sure where formula comes from)
    if (j < (3*index20 - 2*relMaxIndex)) {
      indexLeftEnd   = j;
      radiiEnd.first = positions[index]; 
      break;
    }
  }
  // right limits
  double radiusHigh=0.0, radiusLow=0.0;
  index20 = -1;
  int indexRightEnd = range; // right limit of peak
  // loop on right side of peak
  for (int j=relMaxIndex; j<=range; j++) {
    int index = j + startIndex;
    float value = values[index];
    if (value > 0.2 *maximum) {index20    = j;}
    if (value > 0.85*maximum) {radiusHigh = positions[index];}
    if (value > 0.15*maximum) {radiusLow  = positions[index];}
    if (value > 0.25*maximum) {radii25Perc.second = positions[index];}
    if (value > 0.04*maximum) { radii4Perc.second = positions[index];}
    // if too far out, then break (not sure where formula comes from)
    if (j > (3*index20 - 2*relMaxIndex)) {
      indexRightEnd   = j;
      radiiEnd.second = positions[index];
      break;
    }
  }
  // qDebug() << "width of Peak" << indexRightEnd - indexLeftEnd << "steps";
  rightPeak = (radiusHigh + radiusLow) / 2.;
    
  if (indexRightEnd - indexLeftEnd == 0) { // no peak
    return; // return zeros for the mean and sigma
  }
  double sum=0.0, radialSum=0.0;
  for (int j=indexLeftEnd; j<=indexRightEnd; j++) {
    int index = j + startIndex;
    sum       += values[index];
    radialSum += values[index] * positions[index];
  }
  mean = radialSum / sum;
  double variance = 0.0;
  for (int j=indexLeftEnd; j<=indexRightEnd; j++) {
    int index = j + startIndex;
    float value = values[index];
    double dx = positions[index] - mean;
    variance += value * dx * dx;
  }
  fourSigma = 4 * std::sqrt(variance / sum);
  // qDebug() << "Peak" << i << "maximum at" << peak << "and right peak" << rightPeak;
  // qDebug() << "Mean" << mean << "four sigma" << fourSigma;
}
