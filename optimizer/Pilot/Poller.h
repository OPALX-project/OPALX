#ifndef __POLLER_H__
#define __POLLER_H__

#include <iostream>

#include "Util/MPIHelper.h"

/**
 *  \class Poller
 *  \brief An interface implementing the basics of a poll loop, posting an
 *         MPI_Irecv and waiting for new requests on a specific communicator.
 *
 *  @see Pilot
 *  @see Worker
 *  @see Optimizer
 *  @see MPIHelper.h
 */
class Poller {

public:

    Poller(MPI_Comm comm, double delay = 0.1)
        : comm_m(comm)
        , is_running_(false)
        , poll_delay_(delay)
    {
        last_polled_ = MPI_Wtime();
    }

    virtual ~Poller()
    {}

protected:
    /// communicator the poller listens to requests
    MPI_Comm comm_m;

    bool is_running_;

    /// time of last MPI_Test
    double last_polled_;
    /// delay in seconds between polls
    double poll_delay_;

    /**
     *  User specific behavior on receiving a message.
     *  \return boolean indicating if Irecv has to be re-posted
     */
    virtual bool onMessage(MPI_Status status, size_t recv_value) = 0;
    /// enable implementation to react to STOP tag
    virtual void onStop() = 0;

    /// executed before starting polling loop
    virtual void setupPoll() = 0;
    /// executed before checking for new request
    virtual void prePoll()  = 0;
    /// executed after handling (if any) new request
    virtual void postPoll() = 0;

    /** The poll loop stops when receiving a 'MPI_STOP_TAG' otherwise passes
     *  message to user.
     */
    virtual void run() {
        MPI_Status status;
        MPI_Request req;
        size_t recv_value = 0;
        int flag = 0;

        setupPoll();

        MPI_Irecv(&recv_value, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
                  comm_m, &req);

        is_running_ = true;

        while(true) {

            // regulate the amount of MPI_Test calls (expensive)
            double tnow = MPI_Wtime();
            if(tnow - last_polled_ > poll_delay_)
                last_polled_ = tnow;
            else
                continue;

            prePoll();

            if(req != MPI_REQUEST_NULL) {
                MPI_Test(&req, &flag, &status);
                if(flag) {
                    if(status.MPI_TAG == MPI_STOP_TAG) {
                        is_running_ = false;
                        onStop();
                        return;
                    } else {
                        if(onMessage(status, recv_value))
                            MPI_Irecv(&recv_value, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE,
                                      MPI_ANY_TAG, comm_m, &req);
                        else
                            break;
                    }
                }
            }

            postPoll();
        }
    }
};

#endif