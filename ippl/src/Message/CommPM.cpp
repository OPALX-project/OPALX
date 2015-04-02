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
#include "Message/CommPM.h"
#include "Message/Message.h"
#include "Utility/IpplInfo.h"


// include mpi header file
#include <mpi.h>

// include PM header file. need to invoke PM library.
// #include <Pm.h>
// define interface directly, instead of using Pm.h.
// (because KCC couldn't parse Pm.h)
// The following definitions are written in Pm.h
typedef void* pmCtx;
extern "C"
{
    int _pmGetSendBuf(pmCtx pmc, caddr_t *bufp, size_t length);
    int _pmSend(pmCtx pmc, int dst_node);
    int _pmSendDone(pmCtx pmc);
    int _pmReceive(pmCtx pmc, caddr_t *bufp, size_t *length);
    int _pmPutReceiveBuf(pmCtx pmc);
    extern pmCtx _pm_subnet[];
    extern int _pm_subnet_count;
};
#define PM_MTU                  (8192 + 32)

// include score header file. need to invoke sub-network functions.
#include <score.h>

#include <cstring>
#include <cstdlib>
#include <unistd.h>

// if an error occurs during myreceive more times than this, CommPM
// will just exit.  Make it negative to totally disable checking for a
// maximum number of errors
#define MAX_MPI_ERRS	500


// static data to keep track of errors
static int size_of_MPI_INT; /* needed for tracing */

// sub-network for message passing using PM directly.
static pmCtx pm_network;

////////////////////////////////////////////////////////////////////////////
// constructor.  arguments: command-line args, and number of processes
// to start (if < 0, start the 'default' number, i.e. the number of
// hosts in a MPI virtual machine, the number of nodes in an O2K, etc)
// Note: The base-class constructor does not need the argument info or
// the number of nodes, it just by default sets the number of nodes=1
CommPM::CommPM(int& argc , char**& argv, int procs)
        : Communicate(argc, argv, procs)
{
    
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

    // allocate a sub-network for message passing using PM directly.
    if (_score_alloc_subnet())
    {
        ERRORMSG("CommPM: Error with allocating sub-network.");
    }
    pm_network = _pm_subnet[_pm_subnet_count - 1];

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

    MPI_Type_size ( MPI_INT, (MPI_Aint*)&size_of_MPI_INT );
    if (myHost == 0)      // this code is run by the master process
    {
        // send a messages to each child node
        for (i = 1; i < TotalNodes; i++)
        {
            MPI_Send(&myHost, 1, MPI_INT, i, COMM_HOSTS_TAG, MPI_COMM_WORLD);
            
        }

        // wait for the spawned processes to report back that they're ready
        //~ int *child_ready = new int[TotalNodes];
        std::vector<int> child_ready(TotalNodes);
        for (i = 0; i < TotalNodes; child_ready[i++] = 0);
        INFOMSG("CommPM: Parent process waiting for children ..." << endl);
        reported = 1;		// since the parent is already ready
        while (reported < TotalNodes)
        {
            ierror = MPI_Recv(&rep_host, 1, MPI_INT, MPI_ANY_SOURCE,
                              COMM_HOSTS_TAG, MPI_COMM_WORLD, &stat);
            
            if (rep_host >= 0 && rep_host < TotalNodes && !(child_ready[rep_host]))
            {
                child_ready[rep_host] = 1;
                reported++;
                INFOMSG("CommPM: Child " << rep_host << " ready." << endl);
            }
            else
            {
                ERRORMSG("CommPM: Error with child reporting to parent.  ");
                ERRORMSG("rep_host = " << rep_host);
                ERRORMSG(", child_ready[] = " << child_ready[rep_host] << endl);
            }
        }

        //~ delete [] child_ready;
        INFOMSG("CommPM: Initialization complete." << endl);

    }
    else  			// this is a child process; get data from pops
    {
        char host_name[MPI_MAX_PROCESSOR_NAME];
        ierror = MPI_Get_processor_name(host_name, &result_len);
        if (ierror >= 0)
        {
            INFOMSG("CommPM: Started job " << myHost << " on host `");
            INFOMSG(host_name <<  "'." << endl);
        }
        else
        {
            ERRORMSG("CommPM: failed" << endl);
        }

        // receive message from the master node
        int checknode;
        MPI_Recv(&checknode, 1, MPI_INT, 0, COMM_HOSTS_TAG, MPI_COMM_WORLD,
                 &stat);
        
        if (checknode != 0)
            WARNMSG("CommPM: Child received bad message during startup." << endl);

        // send back an acknowledgement
        MPI_Send(&myHost, 1, MPI_INT, 0, COMM_HOSTS_TAG, MPI_COMM_WORLD);
        
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
CommPM::~CommPM(void)
{
    
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
                if (msg != 0)
                {
                    WARNMSG("CommPM: Found extra message from node " << node);
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
            
        }
    }
    else
    {
        // on client nodes, receive message
        MPI_Recv(&dieCode, 1, MPI_INT, 0, COMM_DIE_TAG, MPI_COMM_WORLD, &stat);
        
    }

    // release sub-network
    _score_free_subnet();

    MPI_Finalize();
}


struct PM_Message
{
    int node;
    int tag;
};

#if 0
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
void *CommPM::pack_message(Message *msg, int tag, int &buffsize)
{
    // calculate size of buffer
    buffsize = find_msg_length(*msg);

    // allocate storage for buffer
    void *pos = (buffsize > PACKSIZE) ? makebuffer(buffsize) : mpipackbuf;

    // pack message data and return the necessary pointer
    fill_msg_buffer(pos, *msg, tag);
    return pos;
}
#endif

////////////////////////////////////////////////////////////////////////////
// send a message ... arguments are the Message itself, the
// destination node, the 'user' tag, and the 'encoding' tag.
// Messages should be sent via the underlying mechanism by using the
// encoding tag (one of the COMM_ tags),
// and should embed the information about what the user
// tag is in the data sent between nodes.  Return success.
bool CommPM::mysend(Message *msg, int node, int tag, int etag)
{

    int nerr = 0;

    //  printf("mysend\n");

    // calculate size of buffer
    int length;
    length = find_msg_length(*msg) + sizeof(PM_Message);

    //  printf("length = %d\n", length);

    // allocate storage for buffer
    PM_Message* msgbuf;
    int timeout_counter = 0;
    while (_pmGetSendBuf(pm_network, (caddr_t*)&msgbuf, length) == ENOBUFS)
    {
        timeout_counter++;
        if (timeout_counter > 10000000)
        {
            ERRORMSG("CommPM: _pmGetSendBuf TIMEOUT");
            timeout_counter = 0;
            nerr++;
        }
        // receive a message when the buffer couldn't be allocated.
        nerr += pickup_message();
    }

    //  printf("buff %x is allocated\n", msgbuf);
    // pack message
    msgbuf->node = myHost;
    msgbuf->tag = tag;
    //  printf("invoke fill_msg_buffer\n");
    fill_msg_buffer((void*) (&(msgbuf->tag) + 1), *msg, tag);
    //  printf("fill_msg_buffer done\n");

    // send message
    if (_pmSend(pm_network, node))
    {
        ERRORMSG("CommPM: _pmSend Error");
        nerr++;
    }
    // receive message waiting the sending done.
    while (_pmSendDone(pm_network) == EBUSY)
    {
        nerr += pickup_message();
    }

    //  printf("mysend done\n");
    // return the succsess of the operation
    return (nerr == 0);
}

int CommPM::pickup_message(void)
{
    int nerr = 0;
    int length;
    PM_Message* msgbuf;
    int error;
    Message* newmsg = 0;
    // int rec_tag;
    int src_node, rec_size, rec_utag;
    void* rec_buff;

    // pickup message
    if (error = _pmReceive(pm_network, (caddr_t*)&msgbuf, (size_t*)&length))
    {
        // no message is received
        if (error != ENOBUFS)
        {
            ERRORMSG("CommPM: _pmReceive Error (in pickup_message)");
            nerr++;
        }
    }
    else
    {
        // a message is received, unpack it
        
        src_node = msgbuf->node;
        // rec_tag = msgbuf->tag;
        rec_size = length - sizeof(PM_Message);
        rec_buff = makebuffer(rec_size);
        //    bcopy((void*) (&(msgbuf->tag) + 1), rec_buff, rec_size);
        memcpy(rec_buff, (void*) (&(msgbuf->tag) + 1), rec_size);
        _pmPutReceiveBuf(pm_network);
        newmsg = unpack_message(src_node, rec_utag, rec_buff);
        newmsg->useCommunicate(this, rec_buff);

        // put message in my message queue
        if (add_msg(newmsg, src_node, rec_utag))
        {
            newmsg = NULL; // reset message pointer
        }
    }

    // return the number of errors;
    return nerr;
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
Message *CommPM::myreceive(int& node, int& tag, int etag)
{
    
    int error;
    PM_Message* msgbuf;
    int length;
    int rec_tag, src_node, rec_size, rec_utag;
    void* rec_buff;
    Message* newmsg = 0;

    // pickup message
    if (error = _pmReceive(pm_network, (caddr_t*)&msgbuf, (size_t*)&length))
    {
        // no message is received
        if (error != ENOBUFS)
        {
            ERRORMSG("CommPM: _pmReceive Error (in myreceive)");
        }
        else
        {
            // no message is available
            DEBUGMSG(level2<<"CommPM: No Message Received to Match Request"<<endl);
        }
    }
    else
    {
        // a message is received, unpack it
        
        src_node = msgbuf->node;
        rec_tag = msgbuf->tag;
        rec_size = length - sizeof(PM_Message);
        rec_buff = makebuffer(rec_size);
        //    bcopy((void*) (&(msgbuf->tag) + 1), rec_buff, rec_size);
        memcpy(rec_buff, (void*) (&(msgbuf->tag) + 1), rec_size);
        _pmPutReceiveBuf(pm_network);
        newmsg = unpack_message(src_node, rec_utag, rec_buff);
        newmsg->useCommunicate(this, rec_buff);
    }

    // Retrun the new Message, or NULL if no message available
    node = src_node;
    tag = rec_tag;
    return newmsg;
}


////////////////////////////////////////////////////////////////////////////
// Synchronize all processors (everybody waits for everybody
// else to get here before returning to calling function).
// Uses MPI barrier for all procs
void CommPM::mybarrier(void)
{
    

    MPI_Barrier(MPI_COMM_WORLD);
}


////////////////////////////////////////////////////////////////////////////
// clean up after a Message has been used (called by Message).
void CommPM::cleanupMessage(void *d)
{
    

    // need to free the allocated storage
    freebuffer(d);
}
