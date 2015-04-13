// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "Message/CommCreator.h"
#include "Message/Communicate.h"


// include files for specific communication libraries


#if defined(IPPL_MPI)
#include "Message/CommMPI.h"
#endif

#if defined(IPPL_SHMEMPI)
#include "Message/CommSHMEMPI.h"
#endif

#if defined(IPPL_PM)
#include "Message/CommPM.h"
#endif

#include <cstring>

// static data for this file
static const char *CommLibraryNames[CommCreator::COMMLIBRARIES] =
  { "pm", "mpi", "shmempi","serial" };

static const char *CommLibraryList = "pm, mpi, shmempi, or serial";


/////////////////////////////////////////////////////////////////////////
// return the name of the Nth library
const char *CommCreator::getLibraryName(int n)
{
    
    if (n >= 0 && n < COMMLIBRARIES)
        return CommLibraryNames[n];
    else
        return 0;
}


/////////////////////////////////////////////////////////////////////////
// return a list of all the libraries, as a single string
const char *CommCreator::getAllLibraryNames()
{
    
    return CommLibraryList;
}


/////////////////////////////////////////////////////////////////////////
bool CommCreator::supported(int cm)
{
    
    if (cm == PM)
    {
#ifdef IPPL_PM
        return true;
#endif
    }
    else if (cm == MPI)
    {
#ifdef IPPL_MPI
        return true;
#endif
    }
    else if (cm == SHMEMPI)
    {
#ifdef IPPL_SHMEMPI
        return true;
#endif
  } 
  else if (cm == SERIAL) {
    return true;
  }
 return false;
}


/////////////////////////////////////////////////////////////////////////
// return the index of the given named library, or (-1) if not found
int CommCreator::libindex(const char *nm)
{
    
    for (int i=0; i < COMMLIBRARIES; ++i)
    {
        if (strcmp(nm, getLibraryName(i)) == 0)
            return i;
    }

    // if here, it was not found
    return (-1);
}


/////////////////////////////////////////////////////////////////////////
// create a new Communicate object.  Arguments = type, argc, argv, num nodes,
// whether to do initialization or not (ignored by some libs).
// If the type code is illegal, or the comm library is not supported,
// return 0.
Communicate *CommCreator::create(int cm, int& argc, char**& argv, int nodes,
                                 bool doinit, MPI_Comm mpicomm)
{

    Communicate *comm = 0;

    // to avoid warning message
    if (doinit) { }

    if (cm == PM)
    {
#ifdef IPPL_PM
        comm = new CommPM(argc, argv, nodes);
#endif
    }
    else if (cm == MPI)
    {
#ifdef IPPL_MPI
        comm = new CommMPI(argc, argv, nodes, doinit, mpicomm);
#endif
    }
    else if (cm == SHMEMPI)
    {
#ifdef IPPL_SHMEMPI
        comm = new CommSHMEMPI(argc, argv, nodes);
#endif
    }
    else if (cm == SERIAL)
    {
        // just make a dummy comm object, which does nothing
        comm = new Communicate(argc, argv, nodes);
    }
    // return the Communicate object
    return comm;
}
