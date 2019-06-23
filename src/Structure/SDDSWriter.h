#ifndef OPAL_SDDS_WRITER_H
#define OPAL_SDDS_WRITER_H

#include <fstream>
#include <string>
#include <queue>
#include <tuple>
#include <utility>
#include <iomanip>

#include <boost/filesystem.hpp>

#include "Algorithms/PartBunchBase.h"

class SDDSWriter {

public:
    // order: text, content
    typedef std::pair<std::string, std::string> desc_t;
    
    // order: name, type, description
    typedef std::tuple<std::string,
                       std::string,
                       std::string> param_t;
    
    // order: mode, no row counts
    typedef std::pair<std::string, size_t> data_t;

    // order: name, type, unit, description
    typedef std::tuple<std::string,
                       std::string,
                       std::string,
                       std::string> cols_t;
    
    
    SDDSWriter(const std::string& fname, bool restart);
    
    virtual void write(PartBunchBase<double, 3>* beam) { };
    
    /** \brief
     *  delete the last 'numberOfLines' lines of the file 'fileName'
     */
    void rewindLines(size_t numberOfLines);
    
    void replaceVersionString();
    
    
    bool exists();
    
protected:
    
    void addDescription(const std::string& text,
                        const std::string& content);
    
    void addParameter(const std::string& name,
                      const std::string& type,
                      const std::string& desc);
    
    void addColumn(const std::string& name,
                   const std::string& type,
                   const std::string& unit,
                   const std::string& desc);
    
    void addInfo(const std::string& mode,
                 const size_t& no_row_counts);

    
    template<typename T>
    void writeValue(const T& value);
    
    void newline();
    
    void open();
    
    void close();
    
    /** \brief Write SDDS header.
     *
     * Writes the appropriate SDDS format header information, The SDDS tools can be used
     * for plotting data.
     */
    void writeHeader();
    
    std::string fname_m;
    
    /** \brief First write to the statistics output file.
     *
     * Initially set to std::ios::out so that SDDS format header information is written to file
     * during the first write call to the statistics output file. Variable is then
     * reset to std::ios::app so that header information is only written once.
     */
    std::ios_base::openmode mode_m;

private:
    
    void writeDescription_m();
    
    void writeParameters_m();
    
    void writeColumns_m();
    
    void writeInfo_m();
    
    std::ofstream os_m;
    
    std::string indent_m;
    
    desc_t desc_m;
    std::queue<param_t> params_m;
    std::queue<cols_t>  cols_m;
    data_t info_m;
};




inline
bool SDDSWriter::exists() {
    return boost::filesystem::exists(fname_m);
}


template<typename T>
void SDDSWriter::writeValue(const T& value) {
    os_m << value << std::setw(10) << "\t";
}


inline
void SDDSWriter::newline() {
    os_m << std::endl;
}


inline
void SDDSWriter::addDescription(const std::string& text,
                                const std::string& content)
{
    desc_m = std::make_pair(text, content);
}


inline
void SDDSWriter::addParameter(const std::string& name,
                              const std::string& type,
                              const std::string& desc)
{
    params_m.push(std::make_tuple(name, type, desc));
}


inline
void SDDSWriter::addInfo(const std::string& mode,
                         const size_t& no_row_counts) {
    info_m = std::make_pair(mode, no_row_counts);
}


inline
void SDDSWriter::addColumn(const std::string& name,
                           const std::string& type,
                           const std::string& unit,
                           const std::string& desc)
{
    cols_m.push(std::make_tuple(name, type, unit, desc));
}

#endif
