#ifndef __PROBEREADER_H__
#define __PROBEREADER_H__

#include <string>
#include <vector>
#include <map>
// #include <locale>

/**
 *  \class SDDSReader
 *  \brief Implements a parser and value extractor for Probe loss files
 */
class ProbeReader {
    
public:
    explicit ProbeReader(std::string filename);
    
    ~ProbeReader();

    void parseFile();
    
    void getVariableValue(int id, std::string varname, double& sim_value);
    
private:
    /// Probe loss filename
    std::string filename_m;
    
    /// Number of variables
    int nColumns_m;
    
    /// Number of values per variable
    int nRows_m;
    
    std::map<std::string, int> columnNamesToID_m;
    std::vector< std::vector<double> > data_m;
    
};

#endif
