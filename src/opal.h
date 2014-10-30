#ifndef __OPAL_H__
#define __OPAL_H__

#include "Ippl.h"
#include "H5hut.h"

#include "AbstractObjects/OpalData.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Parser/TerminalStream.h"
#include "Utilities/ClassicException.h"
#include "Utilities/ParseError.h"
#include "Utilities/Timer.h"
#include "Fields/Fieldmap.hh"

#include <iostream>
#include <vector>
#include <new>
#include <exception>
#include <string>

#include "config.h"

Ippl *ippl;
Inform *gmsg;

int run_opal(char *arg[], std::string inputfile, int restartStep = -1, MPI_Comm comm = MPI_COMM_WORLD);

#endif
