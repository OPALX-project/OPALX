#include <vector>
#include <algorithm>
#include <string.h>
#include "mpi.h"
#include "boost/smart_ptr.hpp"
#include "boost/bind.hpp"

class AsyncSendBuffer {

public:

    AsyncSendBuffer(std::ostringstream &os) {
        this->size_req   = new MPI_Request();
        this->buffer_req = new MPI_Request();
        this->buf_size   = os.str().length();
        buffer           = new char[buf_size];
        memcpy(buffer, os.str().c_str(), buf_size);
    }

    ~AsyncSendBuffer() {
        delete   size_req;
        delete   buffer_req;
        delete[] buffer;
    }

    bool hasCompleted() {
        int bufferflag = 0;
        MPI_Test(this->buffer_req, &bufferflag, MPI_STATUS_IGNORE);
        if(bufferflag) {
            int sizeflag = 0;
            MPI_Test(this->buffer_req, &sizeflag, MPI_STATUS_IGNORE);
            if(sizeflag) {
                return true;
            }
        }
        return false;
    }

    void send(int recv_rank, int size_tag, int data_tag, MPI_Comm comm) {
        MPI_Isend(&buf_size, 1, MPI_LONG, recv_rank, size_tag, comm, size_req);
        MPI_Isend(buffer, buf_size, MPI_CHAR, recv_rank, data_tag, comm, buffer_req);
    }


private:

    // can't use smart pointers because MPI will hold last valid reference to
    // pointer
    MPI_Request *size_req;
    MPI_Request *buffer_req;
    char        *buffer;

    size_t buf_size;

};


class AsyncSendBuffers {

public:
    AsyncSendBuffers() {}

    void insert(boost::shared_ptr<AsyncSendBuffer> buf) {
        collection_.push_back(buf);
    }

    void cleanup() {
        collection_.erase(std::remove_if(
            collection_.begin(), collection_.end(), boost::bind(&AsyncSendBuffer::hasCompleted, _1)),
            collection_.end());
    }

    size_t size() {
        return collection_.size();
    }

private:
    std::vector< boost::shared_ptr<AsyncSendBuffer> > collection_;
};

