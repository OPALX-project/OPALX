// ------------------------------------------------------------------------
// $RCSfile: Stop.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Stop
//   The class for the OPAL STOP command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Stop.h"


// Class Stop
// ------------------------------------------------------------------------

Stop::Stop(): Action(0, "STOP",
                         "The \"STOP\" statement terminates program execution "
                         "or reading of a called file.")
{}


Stop::Stop(const std::string &name, Stop *parent):
    Action(name, parent)
{}


Stop::~Stop()
{}


Stop *Stop::clone(const std::string &name) {
    return new Stop(name, this);
}


void Stop::execute()
{}
