#include "Util/Trace/TraceComponent.h"

#include "mpi.h"

class MyPid : public TraceComponent {

public:


    MyPid(std::string name, MPI_Comm comm)
        : TraceComponent(name)
    {
        mypid_ = 0;
        MPI_Comm_rank(comm, &mypid_);
    }

    void execute(std::ostringstream &dump) {
        dump << mypid_;
    }

private:

    int mypid_;

};
