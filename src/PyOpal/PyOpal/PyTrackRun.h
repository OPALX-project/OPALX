#ifndef PYTRACKRUN_h
#define PYTRACKRUN_h

#include "PyOpal/Globals.h"
#include "PyOpal/PyOpalObject.h"
#include "Track/TrackRun.h"

namespace PyOpal {
namespace PyTrackRunNS {

boost::python::object execute(PyOpalObjectNS::PyOpalObject<TrackRun>& trackRun);
}
}

#endif
