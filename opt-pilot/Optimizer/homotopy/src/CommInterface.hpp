#ifndef CommInterface_hpp
#define CommInterface_hpp

#include "HomotopyTypes.h"

#include <queue>

namespace Homotopy {
    namespace Communication {

        class Interface {
        public:
            virtual void dispatch( const designVars_t &, size_t ) =0;
            virtual void collect(std::queue< std::pair< size_t, objVars_t> >& ) =0;
            virtual void exchange(std::queue< typename Homotopy::NeighborMessage_t > & toSend, std::queue< typename Homotopy::nodeEnvelope_t > & received ) =0;
        };

        template< typename T>
        class SimulatedRemote : public Interface {
            struct req {
                designVars_t x;
                size_t id;
            };
            const int EVALS_AT_ONCE;
            T& P;
            std::vector< req > requests;
        public:
            typedef T BaseType;
            SimulatedRemote( T & p ) : P(p), EVALS_AT_ONCE( 100 ) {
                srand ( unsigned ( time (NULL) ) );
            }
            void dispatch( const designVars_t &x, size_t id ) {
                req r = { x, id };
                requests.push_back( r );
            }
            void collect(std::queue< std::pair< size_t, objVars_t> >& results ){
                int evals = ( (requests.size() >= EVALS_AT_ONCE) ? (EVALS_AT_ONCE) : requests.size() );

                std::vector< req > v( requests.begin(), requests.begin() + evals);
                requests.erase( requests.begin(), requests.begin() + evals );

                //std::random_shuffle( v.begin(), v.end() );

                typename std::vector< req >::iterator it;
                for(it=v.begin(); it!=v.end(); it++){

                    objVars_t r( this->P.size() );
                    int i;
                    for(i=0; i<this->P.size(); i++){
                        r[i] = (this->P[i])( it->x );
                    }

                    std::pair< size_t, objVars_t > p( it->id, r );
                    results.push( p );
                }
            }

            void exchange(std::queue< typename Homotopy::NeighborMessage_t > & toSend, std::queue< typename Homotopy::nodeEnvelope_t > & received ){
                //Do nothing for now, since all simulated mesh subsets are local
            }
        };

        //Communication layer implementations
        //that should allow writing simulated interfaces for testing
        namespace CommImpl {
            
            //Interface to implement
            class Iface {
            public:
                virtual void Init() {}
                virtual void AsyncRecv() {}
                virtual bool HasMessage() { return false;}
                virtual bool ShouldStop() { return false; }
                virtual void Shutdown() {}
                //virtual Status_t GetStatus(){ return 0; }
                virtual size_t GetValue(){ return 0; }
            };

            template< typename T >
            class Simulator : public Iface {
            private:
                struct req {
                    designVars_t x;
                    size_t id;
                };
                std::vector< req > requests;
            public:
                typedef int Status_t;
                virtual void Init() {}
                virtual void AsyncRecv() {}
                virtual bool HasMessage() { return !requests.empty() ;}
                virtual bool ShouldStop() { return false; }
                virtual void Shutdown() {}
                virtual Status_t GetStatus(){ return 0; }
                virtual size_t GetValue(){ return 0; }

                T P;
                void dispatcher( size_t id, const designVars_t &x ) {
                    req r = { x, id };
                    requests.push_back( r );
                }
                bool handler(Status_t status, size_t size, std::queue< std::pair< size_t, objVars_t> >&q, std::queue< typename Homotopy::nodeEnvelope_t > & received){
                    typename std::vector< req >::iterator it;
                    for(it=requests.begin(); it!=requests.end(); it++){

                        objVars_t r( this->P.size() );
                        int i;
                        for(i=0; i<this->P.size(); i++){
                            r[i] = (this->P[i])( it->x );
                        }

                        std::pair< size_t, objVars_t > p( it->id, r );
                        q.push( p );
                    }
                    requests.clear();
                    return true;
                }
                void exchange(std::queue< typename Homotopy::NeighborMessage_t > & toSend ){
                    //Do nothing for now, since all simulated mesh subsets are local
                }
            };
        }

        template< typename T >
        class AdHoc : public Interface {
            //Received evaluations
            std::queue< std::pair< size_t, objVars_t> > q;

            //Received ghost nodes
            std::queue< typename Homotopy::nodeEnvelope_t > g;

            bool HandleMessage() {
                //std::cout << "Handling message..." << std::endl;
                //this->count--;
                return (this->Handler )( comm_.GetStatus(), comm_.GetValue(), this->q, this->g );
            }

            void poll(){
                //std::cout << "Called poll()" << std::endl;
                //comm_.AsyncRecv();
                //Receive all incomming messages
                while( comm_.HasMessage() ) {
                    //std::cout << "hasMessage" << std::endl;
                    if( comm_.ShouldStop() ) {
                        comm_.Shutdown();
                        //TODO: Figure out shutdown protocol
                        return;
                    }
                    else {
                        //std::cout << "ShouldStop = false" << std::endl;
                        if( this->HandleMessage() )
                            comm_.AsyncRecv();
                        else
                            break;
                    }
                }
            }

        public:
            typedef int BaseType;
            T comm_;
            int count;

            void initialize() {
                comm_.Init();
                this->count = 0;
                comm_.AsyncRecv();
            }

            AdHoc() {
                //XXX (IFF): here comm_m is not set!!
                //this->initialize();
            }

            typedef boost::function<void (size_t, const designVars_t &)> dispatcher_t;
            typedef boost::function<bool (typename T::Status_t, size_t, std::queue< std::pair< size_t, objVars_t> >&, std::queue< typename Homotopy::nodeEnvelope_t > &)> handler_t;
            typedef boost::function<void (std::queue< typename Homotopy::NeighborMessage_t > & )> exchanger_t;

            dispatcher_t Dispatcher;
            handler_t Handler;
            exchanger_t Exchanger;

            void shutdown(){
                while( comm_.is_running_ ) this->poll();
            }

            void dispatch( const designVars_t &x , size_t id ) {
                //std::cout << getpid() << ": Called dispatch on ID " << id << std::endl;
                this->poll();
                //while( this->count > 5 ) this->poll();
                //this->count++;
                (this->Dispatcher) (id, x);
                //std::cout << "Completed dispatch on ID " << id << std::endl;
            }

            void collect(std::queue< std::pair< size_t, objVars_t> >& Q ) {
                this->poll();

                while( ! this->q.empty() ) {
                    Q.push( this->q.front() );
                    q.pop();
                }
            }
            void exchange(std::queue< typename Homotopy::NeighborMessage_t > & toSend, std::queue< typename Homotopy::nodeEnvelope_t > & Rec ){
                //Send the local nodes, if we have some
                if( ! toSend.empty() ) (this->Exchanger)( toSend );

                //Poll for other results
                this->poll();

                //Accept the ghost nodes
                while( ! this->g.empty() ) {
                    Rec.push( this->g.front() );
                    g.pop();
                }
            }
        };
    }
}
#endif
