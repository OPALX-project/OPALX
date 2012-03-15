// ------------------------------------------------------------------------
// $RCSfile: DataSink.hh,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DataSink
//   Original, Observer in the Linac code written by  Tim Cleland,
//             Julian Cummings, William Humphrey, and Graham Mark
//             Salman Habib and Robert Ryne
//             Los Alamos National Laboratory
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2003/01/23 13:29:44 $
// $Author: adelmann $
// $Log: DataSink.hh,v $
// Revision 1.1.1.1  2003/01/23 13:29:44  adelmann
// Classic
//
// Revision 3.1  2001/08/30 11:27:22  adelmann
// Clean up, remobe .madprogress stuff
//
// Revision 3.0  2001/08/22 14:41:33  adelmann
// The stable Version
//
// Revision 2.17  2001/08/11 05:30:22  adelmann
// Production Version
//
// Revision 1.1.1.1  2000/11/30 20:29:52  adelmann
// g++ and KCC
//
// Revision 1.3  2000/10/26 05:17:59  adelmann
// Remove DX stuff and add timestam and title
//
// Revision 1.2  2000/08/10 10:56:35  adelmann
// Some cleanup and add a option dx (data explorer) !
// The new printall skript extracts data from this stat file format
//
// Revision 1.1.1.1  2000/07/14 07:20:54  adelmann
// linux version Fri Jul 14 09:15:27 CEST 2000
//
// Revision 1.1.1.1  2000/05/20 11:13:58  adelmann
// Initial working version, thick elements, without: field dump and collimators
//
// Revision 1.3  2000/01/28 07:22:42  adelmann
// Fixrd some bugs with the Mad9pOutput
//
// Revision 1.2  2000/01/27 14:13:29  adelmann
// - Add  bunch->dataSink_m.saveStatDataGnuplotFormat( . )
//        DiscParticle write out
//        updateDotProgress
//
// Revision 1.1.1.1  2000/01/06 07:33:27  adelmann
// linux version works with gcc 991007 V2.96
//
// Revision 2.2  1999/10/29 05:02:07  adelmann
// *** empty log message ***
//
// Revision 2.1  1999/10/27 06:37:00  adelmann
// SGI-LINUX g++, with RF-Gap, REVSCATTER and read distribution from file
//
// Revision 1.1.1.1  1999/10/26 04:22:18  adelmann
// Classic 2.1 (with p)
//
// Revision 1.1.1.1  1999/10/26 04:14:36  adelmann
// Classic 2.1 (with p)
//
//
// ------------------------------------------------------------------------

#ifndef DataSink_H_
#define DataSink_H_

#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunch.h"
#include "Ippl.h"

#include <hdf5.h>
#include "H5Part.h"

/** \brief Class: DataSink
 *
 * This class acts as an observer during the calculation. It generates diagnostic
 * output of the accelerated beam such as statistical beam descriptors of particle
 * positions, momenta, beam phase space (emittance) etc. These are written to file
 * at periodic time steps during the calculation.
 *
 * This class also writes the full beam phase space to an H5 file at periodic time
 * steps in the calculation (this period is different from that of the statistical
 * numbers).
 *
 * Class also writes processor load balancing data to file to track parallel
 * calculation efficiency.
 */

class DataSink
{
public:
    /** \brief Default constructor.
     *
     * The default constructor is called at the start of a new calculation (as
     * opposed to a calculation restart).
     */
    DataSink();

    /** \brief Restart constructor.
     *
     * This constructor is called when a calculation is restarted using data from
     * an existing H5 file.
     */
    DataSink(int restartStep);

    ~DataSink();

private:

    DataSink( const DataSink & ) { }
    DataSink & operator = ( const DataSink & ) { return *this; }

public:

    /** \brief Write H5 file attributes.
     *
     * Called when new H5 is created to initialize file. Write file attributes
     * to describe stored data.
     */
    void writeH5FileAttributes();

    /** \brief Write statistical data.
     *
     * Writes statistical beam data to proper output file. This is information such as RMS beam parameters
     * etc.
     *
     * Also gathers and writes load balancing data to load balance statistics file.
     * \param beam The beam.
     * \param FDext The external E and B field for the head, reference and tail particles. The vector array
     * has the following layout:
     *  - FDext[0] = B at head particle location (in x, y and z).
     *  - FDext[1] = E at head particle location (in x, y and z).
     *  - FDext[2] = B at reference particle location (in x, y and z).
     *  - FDext[3] = E at reference particle location (in x, y and z).
     *  - FDext[4] = B at tail particle location (in x, y, and z).
     *  - FDext[5] = E at tail particle location (in x, y, and z).
     *  \param sposHead Longitudinal position of the head particle.
     *  \param sposRef Longitudinal position of the reference particle.
     *  \param sposTail Longitudinal position of the tail particles.
     */
    void writeStatData(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail);

    /** \brief Write SDDS header.
     *
     * Writes the appropriate SDDS format header information to beam statistics file so the SDDS tools can be used
     * for plotting data.
     * \param outputFile Name of file to write to.
     *
     */
    void writeSDDSHeader(ofstream &outputFile);

    /** \brief Dumps Phase Space to H5 file.
     *
     * \param beam The beam.
     * \param FDext The external E and B field for the head, reference and tail particles. The vector array
     * has the following layout:
     *  - FDext[0] = B at head particle location (in x, y and z).
     *  - FDext[1] = E at head particle location (in x, y and z).
     *  - FDext[2] = B at reference particle location (in x, y and z).
     *  - FDext[3] = E at reference particle location (in x, y and z).
     *  - FDext[4] = B at tail particle location (in x, y, and z).
     *  - FDext[5] = E at tail particle location (in x, y, and z).
     *  \param sposHead Longitudinal position of the head particle.
     *  \param sposRef Longitudinal position of the reference particle.
     *  \param sposTail Longitudinal position of the tail particles.
     */
    void writePhaseSpace(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail);

    /** \brief Dumps phase space to H5 file in OPAL cyclotron calculation.
     *
     * \param beam The beam.
     * \param FDext The external E and B field for the head, reference and tail particles. The vector array
     * has the following layout:
     *  - FDext[0] = B at head particle location (in x, y and z).
     *  - FDext[1] = E at head particle location (in x, y and z).
     *  - FDext[2] = B at reference particle location (in x, y and z).
     *  - FDext[3] = E at reference particle location (in x, y and z).
     *  - FDext[4] = B at tail particle location (in x, y, and z).
     *  - FDext[5] = E at tail particle location (in x, y, and z).
     *  \return Returns the number of the time step just written.
     */
    int writePhaseSpace_cycl(PartBunch &beam, Vector_t FDext[]);
private:

    /** \brief First write to the statistics output file.
     *
     * Initially set to true so that SDDS format header information is written to file
     * during the first write call to the statistics output file. Variable is then
     * reset to false so that header information is only written once.
     */
    bool firstWriteToStat_m;

    /** \brief First write to the H5 file.
     *
     * If true, file attributes and other initialization information are written to file.
     * Variable is then reset to false so that H5 file is only initialized once.
     */
    bool firstWriteH5part_m;

    /// Name of output file for beam statistics.
    string statFileName_m;

    /// Name of output file for processor load balancing information.
    string lBalFileName_m;

    /// %Pointer to H5 file for particle data.
    H5PartFile *H5file_m;

    /// Current record, or time step, of H5 file.
    h5part_int64_t H5call_m;

    /// Timer to track statistics write time.
    IpplTimings::TimerRef StatMarkerTimer_m;

    /// Timer to track particle data/H5 file write time.
    IpplTimings::TimerRef H5PartTimer_m;

    /** \brief Longitudinal shift of reference particle position.
     *
     * Can be used to shift data. For instance, when comparing to IMPACT-T results, which
     * uses a different convention for the reference particle position.
     */
    double sshift_m;

};

#endif // DataSink_H_

/***************************************************************************
 * $RCSfile: DataSink.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 13:29:44 $
 ***************************************************************************/

