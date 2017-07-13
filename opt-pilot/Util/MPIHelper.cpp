#include <string.h>

#include <boost/static_assert.hpp>
#include <boost/serialization/map.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Util/MPIHelper.h"

void serialize(Param_t params, std::ostringstream &os) {

    boost::archive::text_oarchive oa(os);
    oa << params;
}

void serialize(reqVarContainer_t reqvars, std::ostringstream &os) {

    boost::archive::text_oarchive oa(os);
    oa << reqvars;
}

void deserialize(char *buffer, Param_t &params) {

    params.clear();
    std::istringstream is(buffer);
    boost::archive::text_iarchive ia(is);
    ia >> params;
}

void deserialize(char *buffer, reqVarContainer_t &reqvars) {

    reqvars.clear();
    std::istringstream is(buffer);
    boost::archive::text_iarchive ia(is);
    ia >> reqvars;
}

void MPI_Bcast_params(Param_t &params, size_t root, MPI_Comm comm) {

    int pid = 0;
    size_t my_pid = 0;
    MPI_Comm_rank(comm, &pid);
    my_pid = static_cast<size_t>(pid);

    size_t buf_size = 0;
    std::ostringstream os;

    if(my_pid == root) {
        serialize(params, os);
        buf_size = os.str().length();
    }

    MPI_Bcast(&buf_size, 1, MPI_UNSIGNED_LONG, root, comm);

    char *buffer = new char[buf_size];
    if(my_pid == root) memcpy(buffer, os.str().c_str(), buf_size);

    MPI_Bcast(buffer, buf_size, MPI_CHAR, root, comm);
    if(my_pid != root) deserialize(buffer, params);

    delete[] buffer;
}


void MPI_Send_params(Param_t params, size_t pid, MPI_Comm comm) {

    std::ostringstream os;
    serialize(params, os);
    size_t buf_size = os.str().length();

    MPI_Send(&buf_size, 1, MPI_UNSIGNED_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    char *buffer = new char[buf_size];
    memcpy(buffer, os.str().c_str(), buf_size);

    MPI_Send(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    delete[] buffer;
}

std::pair<size_t*, char*> MPI_ISend_params(Param_t params, size_t pid,
                                           MPI_Comm comm, MPI_Request *req) {

    std::ostringstream os;
    serialize(params, os);
    size_t* buf_size = new size_t();
    *buf_size = os.str().length();

    MPI_Isend(buf_size, 1, MPI_UNSIGNED_LONG, pid,
              MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, req);

    char *buffer = new char[*buf_size];
    memcpy(buffer, os.str().c_str(), *buf_size);

    MPI_Isend(buffer, *buf_size, MPI_CHAR, pid,
              MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, req);

    std::pair<size_t*, char*> p(buf_size, buffer);

    return p;
}


void MPI_Recv_params(Param_t &params, size_t pid, MPI_Comm comm) {

    MPI_Status status;
    size_t buf_size = 0;
    MPI_Recv(&buf_size, 1, MPI_UNSIGNED_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    char *buffer = new char[buf_size];

    MPI_Recv(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    deserialize(buffer, params);

    delete[] buffer;
}


void MPI_Send_reqvars(reqVarContainer_t reqvars, size_t pid, MPI_Comm comm) {

    std::ostringstream os;
    serialize(reqvars, os);
    size_t buf_size = os.str().length();

    MPI_Send(&buf_size, 1, MPI_UNSIGNED_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    char *buffer = new char[buf_size];
    memcpy(buffer, os.str().c_str(), buf_size);

    MPI_Send(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm);

    delete[] buffer;
}


void MPI_Recv_reqvars(reqVarContainer_t &reqvars, size_t pid, MPI_Comm comm) {

    MPI_Status status;
    size_t buf_size = 0;
    MPI_Recv(&buf_size, 1, MPI_UNSIGNED_LONG, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    char *buffer = new char[buf_size];

    MPI_Recv(buffer, buf_size, MPI_CHAR, pid,
             MPI_EXCHANGE_SERIALIZED_DATA_TAG, comm, &status);

    deserialize(buffer, reqvars);

    delete[] buffer;
}

