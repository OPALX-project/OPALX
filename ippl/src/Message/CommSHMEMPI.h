// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef COMM_SHMEMPI_H
#define COMM_SHMEMPI_H

/***************************************************************************
 * CommSHMEMPI.h - MPI-specific communications object for use with the
 * Ippl framework.  This differs from CommMPI in that it will, for on-box
 * communication on symmetric multiprocessors, use SHMEM instead of MPI
 * for greater efficiency.  It is basically "MPI + SHMEM on-box".  John
 * Hall is the person to praise or blame for the name SHMEMPI.
 *
 * Allows user to establish id's for available nodes, establish connections,
 * and send/receive data.
 *
 ***************************************************************************/

// include files
#include "Message/Communicate.h"


class CommSHMEMPI : public Communicate
{

public:
    // constructor and destructor
    // constructor arguments: command-line args, and number of processes
    // to start (if < 0, start the 'default' number ... for MPI, this value
    // will be ignored, since the number of nodes is determined by the args
    // to mpirun.
    CommSHMEMPI(int& argc, char**& argv, int procs = (-1));
    virtual ~CommSHMEMPI(void);

    // return the name of this item
    virtual const char *name() const
    {
        return "SHMEMPI";
    }

    //
    //    virtual routines to deal with memory management
    //

    // clean up after a Message has been used (called by Message).
    virtual void cleanupMessage(void *);

protected:

    //
    // implementation-specific routines (which begin with 'my')
    //	these should be provided in a derived class, and contain the
    //	comm-library-specific code
    //

    // send a message ... arguments are the Message itself, the
    // destination node, the 'user' tag, and the 'encoding' tag.
    // Messages should be sent via the underlying mechanism by using the
    // encoding tag (one of the COMM_ tags),
    // and should embed the information about what the user
    // tag is in the data sent between nodes.  Return success.
    virtual bool mysend(Message *, int node, int utag, int etag);

    // receive a message from the given node and user tag.  Return a NEW
    // Message object if a message arrives, or NULL if no message available.
    // node will be set to the node from which the message was sent.
    // tag will be set to the 'user tag' for that message.
    // etag is the 'encoding' tag, and must be one of the COMM_ tags.
    // Only message sent via the underlying mechanism with the
    // given etag are checked.  When one is found, the user tag and sending
    // node are extracted from the sent data.
    // If node = COMM_ANY_NODE, checks for messages from any node.
    // If tag = COMM_ANY_TAG, checks for messages with any user tag.
    virtual Message *myreceive(int& node, int& tag, int etag);

    // Synchronize all processors (everybody waits for everybody
    // else to get here before returning to calling function).
    virtual void mybarrier(void);

private:
    // take data from the given Message, and pack it into the current send buf
    void *pack_message(Message *msg, int tag, int &buffsize);
};

#endif // COMM_SHMEMPI_H
