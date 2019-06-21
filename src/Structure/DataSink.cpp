//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/DataSink.h"

#include "OPALconfig.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Fields/Fieldmap.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPS.h"
#include "Utilities/Timer.h"
#include "Util/SDDSParser.h"

#ifdef __linux__
    #include "MemoryProfiler.h"
#endif


#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif

#include "H5hut.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include <queue>
#include <sstream>

extern Inform *gmsg;

DataSink::DataSink()
    : lossWrCounter_m(0)
{
    this->initWriters_m();
}

DataSink::DataSink(H5PartWrapper *h5wrapper, bool restart)
    : DataSink()
{
    if (restart && !Options::enableHDF5) {
        throw OpalException("DataSink::DataSink()",
                            "Can not restart when HDF5 is disabled");
    }
    
    this->initWriters_m(restart);
    
    H5Writer* h5writer   = dynamic_cast<H5Writer*>(writer_m[Format::H5].get());
    h5writer->changeH5Wrapper(h5wrapper);
    
    if ( restart )
        rewindLines_m();
}

// DataSink::DataSink(H5PartWrapper *h5wrapper):
//     lossWrCounter_m(0)
// {
//     /// Set file write flags to true. These will be set to false after first
//     /// write operation.
//     firstWriteH5Surface_m = true;
//     /// Define file names.
//     std::string fn = OpalData::getInstance()->getInputBasename();
//     surfaceLossFileName_m = fn + std::string(".SurfaceLoss.h5");
//     h5wrapper_m->writeHeader();
//     h5wrapper_m->close();
// }


template<typename ... Arguments>
void DataSink::dump(const Format& format, PartBunchBase<double, 3> *beam, Arguments& ... args) const {
    writer_m[format]->writeData(beam, args ...);
}


void DataSink::storeCavityInformation() {
    if (!Options::enableHDF5) return;
    
    H5Writer* h5writer   = dynamic_cast<H5Writer*>(writer_m[Format::H5].get());
    h5writer->storeCavityInformatio();
}


void DataSink::rewindLines_m() {
    
    if ( Ippl::myNode() != 0 )
        return;
    
    unsigned int linesToRewind = 0;
    
    StatWriter* statwriter = dynamic_cast<StatWriter*>(writer_m[Format::STAT].get());
    H5Writer* h5writer   = dynamic_cast<H5Writer*>(writer_m[Format::H5].get());
    
    // use stat file to get position
    if ( statwriter->exists() ) {
        double spos = h5writer->getLastPosition();
        linesToRewind = statwriter->rewindToSpos(spos);
        statwriter->replaceVersionString();
        
        h5wrapper_m->close();
    }
    
    // rewind all others
    if ( linesToRewind > 0 ) {
        writer_m[Format::LBAL]->rewindLines(linesToRewind);
        writer_m[Format::MEMORY]->rewindLines(linesToRewind);
#ifdef ENABLE_AMR
        writer_m[Format::GRID]->rewindLines(linesToRewind);
#endif
    }
}


void initWriters_m(bool restart) {
    std::string fn = OpalData::getInstance()->getInputBasename();
    writer_m.resize(Format::SIZE);
    writer_m[Format::STAT]      = writer_t(new StatWriter(fn + std::string(".stat"), restart));
    writer_m[Format::LBAL]      = writer_t(new LBalWriter(fn + std::string(".lbal"), restart));
    #ifdef __linux__
    writer_m[Format::MEMORY]    = writer_t(new MemoryProfiler());
#else
    writer_m[Format::MEMORY]    = writer_t(new MemoryWriter(fn + std::string(".mem"), restart));
#endif

#ifdef ENABLE_AMR
    writer_m[Format::GRID]      = writer_t(new GridLBalWriter(fn + std::string(".grid"), restart));
#endif
    writer_m[Format::H5]        = writer_t(new H5Writer());
    
    if ( Options::enableHDF5 ) {
        
    } else {
        writer_m.pop_back();
    }
}


// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End: