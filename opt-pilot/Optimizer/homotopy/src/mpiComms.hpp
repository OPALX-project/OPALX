#ifndef mpiComms_hpp
#define mpiComms_hpp

#include "CommInterface.hpp"

    //typedef MPI_Status Status_t;

    class MPIimpl : public Homotopy::Communication::CommImpl::Iface {
    private:
        struct details {
            MPI_Status status;
            MPI_Request req;
            size_t recv_value;
            int flag;
        };
        std::vector< details > State;
        int active_ind;

        void AsyncRecv( int i ) {
            MPI_Irecv(&State[i].recv_value, 1, MPI_LONG, MPI_ANY_SOURCE,
                        MPI_ANY_TAG, *comm_m[i], &State[i].req);
        }

    public:
        bool is_running_;
        typedef MPI_Status Status_t;
        std::vector< MPI_Comm *> comm_m;

        MPIimpl() : is_running_( false ) {}

        void Init() {
            is_running_ = true;
            active_ind = -1;

            unsigned i;
            for(i=0; i<State.size(); i++) AsyncRecv( i ) ;
        }
        void AsyncRecv() {
            if( active_ind >= 0 ) AsyncRecv( active_ind );
            active_ind = -1;
        }
        bool HasMessage() {
	  unsigned i;
            for(i=0; i<State.size(); i++){
                if (State[i].req != MPI_REQUEST_NULL){
                    MPI_Test(&State[i].req, &State[i].flag, &State[i].status);
                    if(State[i].flag) { 
                        active_ind = i;
                        return true;
                    }
                }
            }
            active_ind = -1;
            return false;
        }
        bool ShouldStop() {
            return State[active_ind].status.MPI_TAG == MPI_STOP_TAG;
        }
        void Shutdown() {

            is_running_ = false;

	    unsigned i;
            for(i=0; i<State.size(); i++) if(i != (unsigned)this->active_ind ) MPI_Cancel( &State[i].req ) ;
        }

        void AddComm( MPI_Comm* c ) {
            comm_m.push_back( c );
            
            details d;
            d.recv_value = 0;
            d.flag = 0;
            State.push_back( d );
        }

        MPI_Status GetStatus(){ return State[active_ind].status; }
        size_t GetValue(){ return State[active_ind].recv_value; }
    };

#endif
