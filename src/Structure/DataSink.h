//
//  Copyright & License: See Copyright.readme in src directory
//

/**
   \brief Class: DataSink

   This class acts as an observer during the calculation. It generates diagnostic
   output of the accelerated beam such as statistical beam descriptors of particle
   positions, momenta, beam phase space (emittance) etc. These are written to file
   at periodic time steps during the calculation.

   This class also writes the full beam phase space to an H5 file at periodic time
   steps in the calculation (this period is different from that of the statistical
   numbers).

   Class also writes processor load balancing data to file to track parallel
   calculation efficiency.
*/

#ifndef _OPAL_DATA_SINK_H
#define _OPAL_DATA_SINK_H

#include <fstream>

#include "Algorithms/PBunchDefs.h"

#include "Writer.h"

template <class T, unsigned Dim>
class PartBunchBase;
class BoundaryGeometry;
class H5PartWrapper;

class DataSink {
public:
    
    enum Format {
        STAT = 0,
        LBAL,
        MEMORY,
#ifdef ENABLE_AMR
        GRID,
#endif
        H5,
        SIZE
    };

    typedef std::vector<std::pair<std::string, unsigned int> > losses_t;
    
    typedef std::unique_ptr<Writer> writer_t;

    /** \brief Default constructor.
     *
     * The default constructor is called at the start of a new calculation (as
     * opposed to a calculation restart).
     */
    DataSink();
    DataSink(H5PartWrapper *h5wrapper, bool restart = false);

    template<typename ... Arguments>
    void dump(const Format& format, PartBunchBase<double, 3> *beam, Arguments& ... args) const;

    /** \brief Write cavity information from  H5 file
     *
     *
     *
     */
    void storeCavityInformation();
    
private:
    void rewindLines_m();
    
    void initWriters_m(bool restart = false);
    
    

public:
    /**
     * Write particle loss data to an ASCII fille for histogram
     * @param fn specifies the name of ASCII file
     * @param beam
     */

    void writePartlossZASCII(PartBunchBase<double, 3> *beam, BoundaryGeometry &bg, std::string fn);

    /**
     * Write geometry points and surface triangles to vtk file
     *
     * @param fn specifies the name of vtk file contains the geometry
     *
     */
    void writeGeomToVtk(BoundaryGeometry &bg, std::string fn);
    //void writeGeoContourToVtk(BoundaryGeometry &bg, std::string fn);

    /**
     * Write impact number and outgoing secondaries in each time step
     *
     * @param fn specifies the name of vtk file contains the geometry
     *
     */
    void writeImpactStatistics(PartBunchBase<double, 3> *beam, long long int &step, size_t &impact, double &sey_num,
                               size_t numberOfFieldEmittedParticles, bool nEmissionMode, std::string fn);

    void writeSurfaceInteraction(PartBunchBase<double, 3> *beam, long long int &step, BoundaryGeometry &bg, std::string fn);

#ifdef ENABLE_AMR
    bool writeAmrStatistics(PartBunchBase<double, 3> *beam);
    
    void noAmrDump(PartBunchBase<double, 3> *beam);
#endif


private:

    DataSink(const DataSink &) { }
    DataSink &operator = (const DataSink &) { return *this; }
    
    
    std::vector<writer_t> writer_m;
    
    static std::string convertToString(int number);

    /** \brief First write to the H5 surface loss file.
     *
     * If true, file name will be assigned and file will be prepared to write.
     * Variable is then reset to false so that H5 file is only initialized once.
     */
    bool firstWriteH5Surface_m;

    /// Name of output file for surface loss data.
    std::string surfaceLossFileName_m;

    /// H5 file for surface loss data.
    h5_file_t H5fileS_m;

    

    /// needed to create index for vtk file
    unsigned int lossWrCounter_m;
};

// inline
// void DataSink::reset() {
//     H5call_m = 0;
// }


inline
std::string DataSink::convertToString(int number) {
    std::stringstream ss;
    ss << std::setw(5) << std::setfill('0') <<  number;
    return ss.str();
}

#endif // DataSink_H_

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End: