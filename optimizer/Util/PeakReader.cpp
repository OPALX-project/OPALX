#include <fstream>
#include <iterator>

#include "Util/PeakReader.h"
#include "Util/OptPilotException.h"

PeakReader::PeakReader(std::string filename)
    : filename_m(filename)
{ }


PeakReader::~PeakReader() { }


void PeakReader::parseFile() {
    
    std::ifstream peak_file;
    
    peak_file.open(filename_m.c_str(), std::ios::in);
    
    if ( !peak_file ) {
        throw OptPilotException("PeakReader::parseFile()",
                                "Error opening file " + filename_m);
    }
    
    // skip header
    std::string header;
    std::getline(peak_file, header);
    
    if ( header.find("# Peak") == std::string::npos ) {
        throw OptPilotException("PeakReader::parseFile()",
                                "Error reading file " + filename_m);
    }
    
    
    int nPeaks = 0;
    peaks_m.clear();
    std::istream_iterator<double> it(peak_file);
    while ( it != std::istream_iterator<double>() ) {
        peaks_m[++nPeaks] = *it;
        ++it;
    }
    
    peak_file.close();
}


void PeakReader::getPeak(int nPeak, double& radius) {
    
    if ( peaks_m.count(nPeak) > 0 ) {
        radius = peaks_m[nPeak];
    } else {
        throw OptPilotException("PeakReader::getPeak",
                                "peak not found!");
    }
}


std::size_t PeakReader::getNumberOfPeaks() {
    return peaks_m.size();
}