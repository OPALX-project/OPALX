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
#include "Message/CommSHMEMPI.h"
#include "Message/Message.h"
#include "Utility/IpplInfo.h"
#include "Profile/Profiler.h"

// include mpi header file.
#include <mpi.h>

#include <cstring>
#include <cstdlib>
#include <unistd.h>


// if an error occurs during myreceive more times than this, CommSHMEMPI
// will just exit.  Make it negative to totally disable checking for a
// maximum number of errors
#define MAX_SHMEMPI_ERRS	500


// static data to keep track of errors
static int numErrors = 0;
static int size_of_SHMEMPI_INT; /* needed for tracing */

// temporary buffer used for speed
#define PSIZE 1024*16
#define PACKSIZE ((PSIZE)*sizeof(long))
static long shmempipackbuf[PSIZE];



////////////////////////////////////////////////////////////////////////////
// constructor.   arguments: command-line args, and number of processes
// to start (if < 0, start the 'default' number ... for MPI, this value
// will be ignored, since the number of nodes is determined by the args
// to mpirun.
CommSHMEMPI::CommSHMEMPI(int& argc , char**& argv, int procs)
        : Communicate(argc, argv, procs)
{
    TAU_PROFILE("CommSHMEMPI::CommSHMEMPI()", "void (int, char **, int)",
                TAU_MESSAGE);

    int i, reported, rep_host, ierror, result_len;
    MPI_Status stat;
    char *currtok, *nexttok, *execname;

    // a little "string magic" to strip the absolute pathname off the executable
    currtok = strstr(argv[0],"/");
    if (!currtok)
    {
        execname = strdup(argv[0]);
    }
    else
    {
        currtok++;
        nexttok = strstr(currtok,"/");
        while (nexttok)
        {
            currtok = nexttok+1;
            nexttok = strstr(currtok,"/");
        }
        execname = strdup(currtok);
    }

    // initialize mpi
    MPI_Init(&argc, &argv);

    // restore original executable name without absolute path
    strcpy(argv[0],execname);

    // determine the number of nodes running and my node number
    MPI_Comm_size(MPI_COMM_WORLD,&TotalNodes);
    MPI_Comm_rank(MPI_COMM_WORLD,&myHost);

    // make sure we do not have too many processes running
    if (procs > 0 && procs < TotalNodes)
    {
        // if this is a process that is beyond what we had requested, just exit
        if (myHost >= procs)
            Ippl::abort();
        TotalNodes = procs;
    }

    MPI_Type_size ( MPI_INT, &size_of_SHMEMPI_INT );
    if (myHost == 0)      // this code is run by the master process
    {
        // send a messages to each child node
        for (i = 1; i < TotalNodes; i++)
        {
            MPI_Send(&myHost, 1, MPI_INT, i, COMM_HOSTS_TAG, MPI_COMM_WORLD);
            TAU_TRACE_SENDMSG(COMM_HOSTS_TAG, i, 1*size_of_SHMEMPI_INT);
        }

        // wait for the spawned processes to report back that they're ready
        std::vector<int> child_ready(TotalNodes);
        for (i = 0; i < TotalNodes; child_ready[i++] = 0);
        INFOMSG("CommSHMEMPI: Parent process waiting for children ..." << endl);
        reported = 1;		// since the parent is already ready
        while (reported < TotalNodes)
        {
            ierror = MPI_Recv(&rep_host, 1, MPI_INT, MPI_ANY_SOURCE,
                              COMM_HOSTS_TAG, MPI_COMM_WORLD, &stat);
            TAU_TRACE_RECVMSG(COMM_HOSTS_TAG, stat.MPI_SOURCE,
                              1*size_of_SHMEMPI_INT);
            if (rep_host >= 0 && rep_host < TotalNodes && !(child_ready[rep_host]))
            {
                child_ready[rep_host] = 1;
                reported++;
                INFOMSG("CommSHMEMPI: Child " << rep_host << " ready." << endl);
            }
            else
            {
                ERRORMSG("CommSHMEMPI: Error with child reporting to parent.  ");
                ERRORMSG("rep_host = " << rep_host);
                ERRORMSG(", child_ready[] = " << child_ready[rep_host] << endl);
            }
        }

        INFOMSG("CommSHMEMPI: Initialization complete." << endl);

    }
    else  			// this is a child process; get data from pops
    {
        char host_name[MPI_MAX_PROCESSOR_NAME];
        ierror = MPI_Get_processor_name(host_name, &result_len);
        if (ierror >= 0)
        {
            INFOMSG("CommSHMEMPI: Started job " << myHost << " on host `");
            INFOMSG(host_name <<  "'." << endl);
        }
        else
        {
            ERRORMSG("CommSHMEMPI: failed" << endl);
        }

        // receive message from the master node
        int checknode;
        MPI_Recv(&checknode, 1, MPI_INT, 0, COMM_HOSTS_TAG, MPI_COMM_WORLD,
                 &stat);
        TAU_TRACE_RECVMSG(COMM_HOSTS_TAG, 0, 1*size_of_SHMEMPI_INT);
        if (checknode != 0)
            WARNMSG("CommSHMEMPI: Child received bad message during startup." << endl);

        // send back an acknowledgement
        MPI_Send(&myHost, 1, MPI_INT, 0, COMM_HOSTS_TAG, MPI_COMM_WORLD);
        TAU_TRACE_SENDMSG(COMM_HOSTS_TAG, 0, 1*size_of_SHMEMPI_INT);
    }

    // set up the contexts and processes arrays properly
    if (TotalNodes > 1)
    {
        vector<int> proccount;
        proccount.push_back(1);
        for (i = 1; i < TotalNodes; i++)
        {
            Contexts.push_back(1);
            Processes.push_back(proccount);
        }
    }
}


////////////////////////////////////////////////////////////////////////////
// class destructor
CommSHMEMPI::~CommSHMEMPI(void)
{
    TAU_PROFILE("CommSHMEMPI::~CommSHMEMPI()", "void()", TAU_MESSAGE);
    int i, dieCode = 0;
    MPI_Status stat;

    // on all nodes, when running in parallel, get any extra messages not
    // yet received
    if (TotalNodes > 1)
    {
        int trial, node, tag;
        Message *msg;
        for (trial = 0; trial < 50000; ++trial)
        {
            do
            {
                node = COMM_ANY_NODE;
                tag = COMM_ANY_TAG;
                msg = myreceive(node, tag, COMM_SEND_TAG);
                if (msg != 0 && tag != IPPL_ABORT_TAG && tag != IPPL_EXIT_TAG)
                {
                    WARNMSG("CommSHMEMPI: Found extra message from node " << node);
                    WARNMSG(", tag " << tag << ": msg = " << *msg << endl);
                }
            }
            while (msg != 0);
        }
    }

    // broadcast a message to all other nodes to tell them to quit
    if (myNode() == 0)
    {
        // on master node, send out messages
        for (i = 1; i < TotalNodes; i++)
        {
            MPI_Send(&dieCode, 1, MPI_INT, i, COMM_DIE_TAG, MPI_COMM_WORLD);
            TAU_TRACE_SENDMSG(COMM_DIE_TAG, i, 1*size_of_SHMEMPI_INT);
        }
    }
    else
    {
        // on client nodes, receive message
        MPI_Recv(&dieCode, 1, MPI_INT, 0, COMM_DIE_TAG, MPI_COMM_WORLD, &stat);
        TAU_TRACE_RECVMSG(COMM_DIE_TAG, 0, 1*size_of_SHMEMPI_INT);
    }

    MPI_Finalize();
}


////////////////////////////////////////////////////////////////////////////
// take the data from a Message object and pack it into the current send buf.
// each message is packed in this order:
//      tag, sending node, number of items             (3-int array)
//              type of item 1  (short)
//              size of item 1, in number of elements   (int)
//              item 1 data     (various)
//              ...
//              type of item N  (short)
//              size of item N, in number of elements   (int)
//              item N data     (various)
void *CommSHMEMPI::pack_message(Message *msg, int tag, int &buffsize)
{
    TAU_PROFILE("CommSHMEMPI::pack_message()", "(Message *, int, int)",
                TAU_MESSAGE);
    // calculate size of buffer
    buffsize = find_msg_length(*msg);

    // allocate storage for buffer
    void *pos = (buffsize > PACKSIZE) ? makebuffer(buffsize) : shmempipackbuf;

    // pack message data and return the necessary pointer
    fill_msg_buffer(pos, *msg, tag, buffsize);
    return pos;
}


////////////////////////////////////////////////////////////////////////////
// send a message ... arguments are the Message itself, the
// destination node, the 'user' tag, and the 'encoding' tag.
// Messages should be sent via the underlying mechanism by using the
// encoding tag (one of the COMM_ tags),
// and should embed the information about what the user
// tag is in the data sent between nodes.  Return success.
bool CommSHMEMPI::mysend(Message *msg, int node, int tag, int etag)
{
    TAU_PROFILE("CommSHMEMPI::mysend()", "bool (Message *, int, int, int)",
                TAU_MESSAGE);

    int nitems = msg->size();
    int errstat = (-1);
    int flag = false;
    MPI_Request request;
    MPI_Status status;

    MPI_Status rec_status;
    int src_node, rec_node, rec_tag, rec_size, rec_utag, bufid, rec_flag = 0;
    Message* newmsg = NULL;

    // pack the message data into the buffer
    int size;
    void *outbuffer = pack_message(msg, tag, size);

    // send the message (non-blocking)
    errstat = MPI_Isend(outbuffer, size, MPI_BYTE, node, etag,
                        MPI_COMM_WORLD, &request);
    TAU_TRACE_SENDMSG(etag, node, size);

    while (!flag)
    {
        // get info about messages to be received
        bufid = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                           &rec_flag, &rec_status);
        if ( (bufid >= 0) && (rec_flag != 0) )
        {
            // a message is available to be received
            src_node = rec_status.MPI_SOURCE;
            rec_tag = rec_status.MPI_TAG;
            MPI_Get_count(&rec_status, MPI_BYTE, &rec_size);
            if ( (rec_size >= 0) && (rec_tag >= 0) && (src_node >= 0) )
            {
                // message is a valid one, so malloc the output buffer
                void *rec_buff = makebuffer(rec_size);

                // blocking receive, unpack message
                MPI_Recv(rec_buff, rec_size, MPI_BYTE, src_node, rec_tag,
                         MPI_COMM_WORLD, &rec_status);
                TAU_TRACE_RECVMSG(rec_tag, src_node, rec_size);
                newmsg = unpack_message(rec_node, rec_utag, rec_buff);

                // tell this new Message that we were the one that created its
                // storage buffer, so that when the Messageis deleted, we can
                // be told about it in order to free the storage.
                newmsg->useCommunicate(this, rec_buff);

                // put message in my message queue
                if (add_msg(newmsg,rec_node,rec_utag))
                {
                    newmsg = NULL; // reset message pointer
                    rec_flag = 0; // reset receive flag
                }
            }
        }

        // check for completion of send
        MPI_Test(&request, &flag, &status);
    }

    //  free up the send buffer
    if (size > PACKSIZE)
        freebuffer(outbuffer);

    // return the success of the operation
    return (errstat == 0);
}


////////////////////////////////////////////////////////////////////////////
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
Message *CommSHMEMPI::myreceive(int& node, int& tag, int etag)
{
    TAU_PROFILE("CommSHMEMPI::myreceive()", "Message *(int, int, int)",
                TAU_MESSAGE);

    int bufid, size, checknode, checktag, flag = false;
    Message *newmsg = 0;
    MPI_Status stat;

    checknode = (node < 0 || node >= TotalNodes ? MPI_ANY_SOURCE : node);
    checktag = etag;

    // get info about message
    bufid = MPI_Iprobe(checknode, checktag, MPI_COMM_WORLD, &flag, &stat);
    if (bufid < 0)
    {
        // an error has occurred
        ERRORMSG("CommSHMEMPI: cannot receive msg from node " << checknode);
        ERRORMSG(", tag " << checktag << endl);

        if (MAX_SHMEMPI_ERRS > 0 && ++numErrors > MAX_SHMEMPI_ERRS)
        {
            ERRORMSG("Maximum number of MPI receive errors (" << numErrors);
            ERRORMSG(") exceeded. MPI is hosed!!" << endl);
            Ippl::abort();
        }
    }

    // if the message is actually available, see if we can get it now
    if (flag == true)
    {
        MPI_Get_count(&stat,MPI_BYTE,&size);
        if (size < 0)
        {
            ERRORMSG("CommSHMEMPI: received message has size " << size << endl);
        }
        else if ((stat.MPI_TAG != checktag) || (stat.MPI_TAG < 0))
        {
            ERRORMSG("CommSHMEMPI: received message with invalid tag ");
            ERRORMSG(stat.MPI_TAG << endl);
        }
        else if (stat.MPI_SOURCE < 0)
        {
            ERRORMSG("CommSHMEMPI: received message from invalid source ");
            ERRORMSG(stat.MPI_SOURCE << endl);
        }
        else
        {
            checknode = stat.MPI_SOURCE;
            checktag = stat.MPI_TAG;

            // malloc the receive buffer
            void *outbuff = makebuffer(size);

            // blocking receive
            MPI_Recv(outbuff, size, MPI_BYTE, checknode, checktag,
                     MPI_COMM_WORLD, &stat);
            TAU_TRACE_RECVMSG(checktag, checknode, size);
            newmsg = unpack_message(node, tag, outbuff);

            // tell this new Message that we were the one that created its
            // storage buffer, so that when the Messageis deleted, we can
            // be told about it in order to free the storage.
            newmsg->useCommunicate(this, outbuff);
            numErrors = 0;
        }
    }
    else
    {
        // no message is available
        DEBUGMSG(level2 << "CommSHMEMPI: No Message Received to Match Request");
        DEBUGMSG(endl);
    }

    // return the new Message, or NULL if no message available
    return newmsg;
}


////////////////////////////////////////////////////////////////////////////
// Synchronize all processors (everybody waits for everybody
// else to get here before returning to calling function).
// Uses MPI barrier for all procs
void CommSHMEMPI::mybarrier(void)
{
    TAU_PROFILE("CommSHMEMPI::mybarrier()", "void ()", TAU_MESSAGE);

    MPI_Barrier(MPI_COMM_WORLD);
}


////////////////////////////////////////////////////////////////////////////
// clean up after a Message has been used (called by Message).
void CommSHMEMPI::cleanupMessage(void *d)
{
    TAU_PROFILE("CommSHMEMPI::cleanupMessage()", "void (void *)", TAU_MESSAGE);

    // need to free the allocated storage
    freebuffer(d);
}
