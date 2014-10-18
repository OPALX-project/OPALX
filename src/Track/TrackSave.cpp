// ------------------------------------------------------------------------
// $RCSfile: TrackSave.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackSave
//   The class for the OPAL TSAVE command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:47 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Track/TrackSave.h"
#include "Algorithms/PartBunch.h"
#include "Attributes/Attributes.h"
#include "Track/Track.h"
#include "Utilities/OpalException.h"
#include <fstream>
#include <iomanip>


// Class TrackSave
// ------------------------------------------------------------------------


TrackSave::TrackSave():
    Action(1, "TSAVE",
           "The \"TSAVE\" sub-command saves the defined particles "
           "on the given file.") {
    itsAttr[0] = Attributes::makeString
                 ("FILE", "Name of file to be written", "TRACKSAVE");
}


TrackSave::TrackSave(const std::string &name, TrackSave *parent):
    Action(name, parent)
{}


TrackSave::~TrackSave()
{}


TrackSave *TrackSave::clone(const std::string &name) {
    return new TrackSave(name, this);
}


void TrackSave::execute() {
    // open output file.
    std::string file = Attributes::getString(itsAttr[0]);
    std::ofstream os(file.c_str());

    if(os.bad()) {
        throw OpalException("TrackSave::execute()",
                            "Unable to open output file \"" + file + "\".");
    }

    os << "\nSaved particle positions:\n";
    std::streamsize old_prec = os.precision(8);
    os.setf(std::ios::fixed, std::ios::floatfield);
    PartBunch *bunch = Track::block->bunch;

    for(unsigned int i = 0; i < bunch->getLocalNum(); i++) {
        Particle part = bunch->get_part(i);
        os << part.x() << ' ' << part.px() << ' '
           << part.y() << ' ' << part.py() << ' '
           << part.t() << ' ' << part.pt() << '\n';
    }

    os << std::flush;
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}
