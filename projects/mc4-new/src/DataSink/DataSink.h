// ------------------------------------------------------------------------
// $RCSfile: DataSink.hh,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DataSink
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2003/01/23 13:29:44 $
// $Author: adelmann $
// $Log: DataSink.hh,v $
// Revision 1.1.1.1  2003/01/23 13:29:44  adelmann
//
// ------------------------------------------------------------------------

#ifndef DataSink_H_
#define DataSink_H_

#ifdef IPPL_USE_SINGLE_PRECISION
 #define TT float
#else
 #define TT double
#endif

#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

#include "Ippl.h"

template<class T, unsigned int Dim> class ChargedParticles;

#include <hdf5.h>
#include "H5hut.h"



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
    DataSink(string fn, ChargedParticles<TT,3> *univ);

    /** \brief Restart constructor.
     *
     * This constructor is called when a calculation is restarted using data from
     * an existing H5 file or start at step 0 when reading initial conditions from file.
     */
    DataSink(string fn, ChargedParticles<TT,3> *univ, int restartStep);

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

    /** \brief Dumps Phase Space to H5 file.
     *
     *  \param univ The universe to be dumped.
     */
    void writePhaseSpace(TT time, TT z, int step);


    /** \brief Read Phase Space from H5 file.
     *
     *  \param univ The universe to be filled.
     */
    void readPhaseSpace();


private:

    /** \brief First write to the H5 file.
     *
     * If true, file attributes and other initialization information are written to file.
     * Variable is then reset to false so that H5 file is only initialized once.
     */
    bool firstWriteH5part_m;

    /// Name of output file for processor load balancing information.
    string lBalFileName_m;

    /// %Pointer to H5 file for particle data.
    h5_file_t *H5file_m;

    /// Current record, or time step, of H5 file.
    h5_int64_t H5call_m;

    /// Timer to track particle data/H5 file write time.
    IpplTimings::TimerRef H5PartTimer_m;
    
    ChargedParticles<TT,3> *univ_m;

};

#endif // DataSink_H_

/***************************************************************************
 * $RCSfile: DataSink.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 13:29:44 $
 ***************************************************************************/

