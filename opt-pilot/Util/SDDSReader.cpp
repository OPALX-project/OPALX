#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "Util/SDDSReader.h"
#include "Util/OptPilotException.h"

SDDSReader::SDDSReader(std::string filename) {

    filename_m = filename;
    numParams_m = 0;
    numColumns_m = 0;
}

SDDSReader::~SDDSReader()
{}


void SDDSReader::parseFile() {
    std::ifstream sdds;
    sdds.open(filename_m.c_str(), std::ios::in);
    if(!sdds) {
        throw OptPilotException("SDDSReader::parseFile()",
                                "Error opening file " + filename_m);
    }

    std::string header;
    std::getline(sdds, header, '\n');
    if( !header.compare("SDDS1\n") ) {
        throw OptPilotException("SDDSReader::parseFile()",
                                "Error parsing SDDS header!");
    }

    // parse header
    while(sdds) {
        std::string linestart;
        std::getline(sdds, linestart, ' ');

        if( linestart.compare(0, 12, "&description") == 0) {
            std::getline(sdds, linestart, '&');
            std::getline(sdds, linestart, '\n');
        } else if( linestart.compare(0, 10, "&parameter") == 0) {

            std::getline(sdds, linestart, '&');
            std::vector<std::string> tmp = split(linestart, ",");
            std::string ts = tmp[0];
            tmp.clear();
            tmp = split(ts, "=");

            addParameter(tmp[1]);

            std::getline(sdds, linestart, '\n');

        } else if( linestart.compare(0, 7, "&column") == 0 ) {

            std::getline(sdds, linestart, '&');
            std::vector<std::string> tmp = split(linestart, ",");
            std::string ts = tmp[0];
            tmp.clear();
            tmp = split(ts, "=");

            addColumn(tmp[1]);

            std::getline(sdds, linestart, '\n');

        } else if( linestart.compare(0, 5, "&data") == 0 ) {
            std::getline(sdds, linestart, '&');
            std::getline(sdds, linestart, '\n');
            break;
        } else
            throw OptPilotException("SDDSReader::parseFile()",
                            "Header line starts with unrecognized field: '" + linestart + "'!");
    }

    // read params
    for(int i=0; i<numParams_m; i++) {
        std::string param;
        std::getline(sdds, param, '\n');
        params_m.push_back(param);
    }

    // read all dumps
    while(sdds) {
        std::string line;
        std::vector<std::string> tmp;
        for(int i=0; i<numColumns_m-1; i++) {
            std::getline(sdds, line, '\t');
            tmp.push_back(line);
        }
        std::getline(sdds, line, '\n');
        tmp.push_back(line);
        columns_m.push_back(tmp);
    }
}


std::vector<std::string> SDDSReader::split(
        std::string str, std::string delimiter) {

    std::string tmp;
    std::vector<std::string> res;

    size_t fpos = str.find(delimiter);
    while(fpos != std::string::npos) {
        tmp = str.substr(0,fpos);
        str = str.substr(fpos+1, std::string::npos);
        res.push_back(tmp);
        fpos = str.find(delimiter);
    }
    res.push_back(str);

    return res;
}


//XXX use either all upper, or all lower case chars
void SDDSReader::fixCaseSensitivity(std::string &for_string) {

    boost::to_lower(for_string);
}


void SDDSReader::addParameter(std::string param_name) {

    fixCaseSensitivity(param_name);
    paramNameToID_m.insert(std::pair<std::string, int>(param_name, numParams_m));
    numParams_m++;
}


void SDDSReader::addColumn(std::string column_name) {

    fixCaseSensitivity(column_name);
    columnNameToID_m.insert(std::pair<std::string, int>(column_name, numColumns_m));
    numColumns_m++;
}