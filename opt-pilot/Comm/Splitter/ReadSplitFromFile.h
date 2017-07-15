#ifndef __READ_SPLIT_FROM_FILE__
#define __READ_SPLIT_FROM_FILE__

#include "mpi.h"

#include "Util/OptPilotException.h"
#include "Comm/Splitter.h"

class ReadSplitFromFile : public Splitter {

public:

    ReadSplitFromFile(MPI_Comm comm = MPI_COMM_WORLD) : comm_m(comm) {
        MPI_Comm_rank(comm_m, &rank_m);
        MPI_Comm_size(comm_m, &size_m);
        iAmOptimizer_m = false;
        iAmWorker_m = false;
        iAmPilot_m = false;

        if(size_m < 3)
            throw OptPilotException("CommSplitter::CommSplitter", "Need 3 or more cores to split!");

        split();
    }

    virtual ~ReadSplitFromFile()
    {}

    void freeComm() {
        MPI_Group_free(&gpw);
        MPI_Group_free(&gpo);
        MPI_Group_free(&all);
        //MPI_Comm_free(&pilot_opt);
        //MPI_Comm_free(&pilot_worker);
    }

    bool isOptimizer() { return iAmOptimizer_m; }
    bool isWorker()    { return iAmWorker_m; }
    bool isPilot()     { return iAmPilot_m; }

    MPI_Comm myInternalComm() {
        return myInternalComm_m;
    }

    MPI_Comm myExternalComm() {
        return myExternalComm_m;
    }

    MPI_Comm PilotListenerComm() {
        if(isPilot())
            return PilotListener_m;
        else
            return MPI_COMM_NULL;
    }

    MPI_Comm PilotOptimizerComm() {
        if(isPilot())
            return myInternalComm_m;
        else
            return MPI_COMM_NULL;
    }

    MPI_Comm PilotWorkerComm() {
        if(isPilot())
            return myExternalComm_m;
        else
            return MPI_COMM_NULL;
    }

    int PilotRank() {
        return 0;
    }

    int OptMaster() {
        return 1;
    }



private:

    //typedef std::vector<MPI_Comm> Communicators;

    MPI_Comm comm_m;
    MPI_Comm PilotListener_m;
    MPI_Comm myInternalComm_m;
    MPI_Comm myExternalComm_m;

    MPI_Group all;
    MPI_Group gpo;
    MPI_Group gpw;
    MPI_Comm pilot_opt;
    MPI_Comm pilot_worker;

    int rank_m;
    int size_m;

    bool iAmOptimizer_m;
    bool iAmWorker_m;
    bool iAmPilot_m;


    void split() {

        MPI_Comm_group(comm_m, &all);

        int rpilot_opt[2] = { 0, 1 };
        MPI_Group_incl(all, 2, rpilot_opt, &gpo);
        MPI_Comm_create(comm_m, gpo, &pilot_opt);

        //XXX: for now every worker has only one core
        int ropt[1] = { 1 };
        MPI_Group_excl(all, 1, ropt, &gpw );
        MPI_Comm_create(comm_m, gpw, &pilot_worker);


        PilotListener_m = MPI_COMM_WORLD;

        if(rank_m == 0) {

            myInternalComm_m = pilot_opt;
            myExternalComm_m = pilot_worker;
            iAmPilot_m = true;

        } else if(rank_m == 1) {

            myInternalComm_m = MPI_COMM_SELF;
            myExternalComm_m = pilot_opt;
            iAmOptimizer_m = true;

        } else {

            myInternalComm_m = MPI_COMM_SELF;
            myExternalComm_m = pilot_worker;
            iAmWorker_m = true;
        }
    }

};

#endif
