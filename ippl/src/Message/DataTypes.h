#ifndef IPPL_MPI_DATATYPES_H
#define IPPL_MPI_DATATYPES_H

/*
 * Implementation like Boost.MPI
 */

#include <mpi.h>

template <typename T> MPI_Datatype get_mpi_datatype(const T& x)
{
    return get_mpi_datatype(T());
}


#define IPPL_MPI_DATATYPE(CppType, MPIType)                       \
template<>                                                        \
inline MPI_Datatype                                               \
get_mpi_datatype< CppType >(const CppType&) { return MPIType; }


IPPL_MPI_DATATYPE(char, MPI_CHAR);

IPPL_MPI_DATATYPE(short, MPI_SHORT);

IPPL_MPI_DATATYPE(int, MPI_INT);

IPPL_MPI_DATATYPE(long, MPI_LONG);

IPPL_MPI_DATATYPE(float, MPI_FLOAT);

IPPL_MPI_DATATYPE(double, MPI_DOUBLE);

IPPL_MPI_DATATYPE(long double, MPI_LONG_DOUBLE);

IPPL_MPI_DATATYPE(unsigned char, MPI_UNSIGNED_CHAR);

IPPL_MPI_DATATYPE(unsigned short, MPI_UNSIGNED_SHORT);

IPPL_MPI_DATATYPE(unsigned, MPI_UNSIGNED);

IPPL_MPI_DATATYPE(unsigned long, MPI_UNSIGNED_LONG);

IPPL_MPI_DATATYPE(bool, MPI_CXX_BOOL);

#endif
