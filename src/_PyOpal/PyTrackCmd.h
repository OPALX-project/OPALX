#ifndef PYTRACKCMD_h
#define PYTRACKCMD_h

#include "PyOpal/Globals.h"
#include "PyOpal/PyOpalObject.h"
#include "Track/TrackCmd.h"

namespace PyOpal {
namespace PyTrackCmdNS {

void executeWrapper(PyOpalObject<TrackCmd>& cmd);

}
}

#endif
