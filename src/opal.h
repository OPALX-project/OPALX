#ifndef __OPAL_H__
#define __OPAL_H__

#include "Ippl.h"


Ippl *ippl;
Inform *gmsg;

int run_opal(char *arg[], std::string inputfile, int restartStep = -1, MPI_Comm comm = MPI_COMM_WORLD);

#endif
