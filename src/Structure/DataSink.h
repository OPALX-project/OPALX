//
// Class DataSink
//   This class acts as an observer during the calculation. It generates diagnostic
//   output of the accelerated beam such as statistical beam descriptors of particle
//   positions, momenta, beam phase space (emittance) etc. These are written to file
//   at periodic time steps during the calculation.
//
//   This class also writes the full beam phase space to an H5 file at periodic time
//   steps in the calculation (this period is different from that of the statistical
//   numbers).

//   Class also writes processor load balancing data to file to track parallel
//   calculation efficiency.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef _OPAL_DATA_SINK_H
#define _OPAL_DATA_SINK_H

#include "Algorithms/Vektor.h"
#include "Structure/H5Writer.h"
#include "Structure/SDDSWriter.h"
#include "Structure/StatWriter.h"

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

template <class T, unsigned Dim>
class PartBunch;
class BoundaryGeometry;
class H5PartWrapper;

class DataSink {
private:
    typedef StatWriter::losses_t            losses_t;
    typedef std::unique_ptr<StatWriter>     statWriter_t;
    typedef std::unique_ptr<SDDSWriter>     sddsWriter_t;
    typedef std::unique_ptr<H5Writer>       h5Writer_t;
    
public:
    /** \brief Default constructor.
     *
     * The default constructor is called at the start of a new calculation (as
     * opposed to a calculation restart).
     */
    DataSink();
    DataSink(H5PartWrapper *h5wrapper, bool restart);
    DataSink(H5PartWrapper *h5wrapper);

    void dumpH5(PartBunch<double, 3> *beam, Vector_t<double, 3> FDext[]) const;
    
    int dumpH5(PartBunch<double, 3> *beam, Vector_t<double, 3> FDext[], double meanEnergy,
               double refPr, double refPt, double refPz,
               double refR, double refTheta, double refZ,
               double azimuth, double elevation, bool local) const;
    
    void dumpSDDS(PartBunch<double, 3> *beam, Vector_t<double, 3> FDext[],
                  const double& azimuth = -1) const;
    
    void dumpSDDS(PartBunch<double, 3> *beam, Vector_t<double, 3> FDext[],
                  const losses_t &losses = losses_t(), const double& azimuth = -1) const;
    
    /** \brief Write cavity information from  H5 file
     */
    void storeCavityInformation();
    
    void changeH5Wrapper(H5PartWrapper *h5wrapper);
    
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
    void writeImpactStatistics(const PartBunch<double, 3> *beam,
                               long long int &step,
                               size_t &impact,
                               double &sey_num,
                               size_t numberOfFieldEmittedParticles,
                               bool nEmissionMode,
                               std::string fn);

private:
    DataSink(const DataSink& ds) = delete;
    DataSink &operator = (const DataSink &) = delete;

    void rewindLines();

    void init(bool restart = false,
              H5PartWrapper* h5wrapper = nullptr);


    h5Writer_t      h5Writer_m;
    statWriter_t    statWriter_m;
    std::vector<sddsWriter_t> sddsWriter_m;

    static std::string convertToString(int number, int setw = 5);

    /// needed to create index for vtk file
    unsigned int lossWrCounter_m;

    /// Timer to track statistics write time.
    IpplTimings::TimerRef StatMarkerTimer_m;

};


inline
std::string DataSink::convertToString(int number, int setw) {
    std::stringstream ss;
    ss << std::setw(setw) << std::setfill('0') <<  number;
    return ss.str();
}


#endif // DataSink_H_
