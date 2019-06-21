#ifndef OPAL_SDDS_WRITER_H
#define OPAL_SDDS_WRITER_H

#include <fstream>
#include <string>

#include "Structure/Writer.h"


class SDDSWriter : public Writer {

public:
    SDDSWriter(const std::string& fname, const double& spos);
    
    /** \brief
     *  delete the last 'numberOfLines' lines of the file 'fileName'
     */
    void rewindLines(size_t numberOfLines);
    
    void replaceVersionString() const;
    
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

    void addData(const std::string& mode,
                 const short& no_row_counts);

    /** \brief Write SDDS header.
     *
     * Writes the appropriate SDDS format header information, The SDDS tools can be used
     * for plotting data.
     */
    virtual void writeHeader() = 0;
    
    virtual void writeData(PartBunchBase<double, 3> *beam) = 0;
    
    template<typename T>
    void writeValue(const T& value);
    
    
    void open_m();
    
    void close_m();
    
protected:
    /** \brief First write to the statistics output file.
     *
     * Initially set to std::ios::out so that SDDS format header information is written to file
     * during the first write call to the statistics output file. Variable is then
     * reset to std::ios::app so that header information is only written once.
     */
    std::ios_base::openmode mode_m;
    
    std::ofstream os_m;
};

#endif
