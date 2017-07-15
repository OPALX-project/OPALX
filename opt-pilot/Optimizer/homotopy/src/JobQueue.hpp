#ifndef JobQueue_hpp
#define JobQueue_hpp

#include "HomotopyTypes.h"

#include <queue>
#include <map>
#include <utility>

//For testing
#include <ctime>
#include <cstdlib>
#include <algorithm>

namespace Homotopy {
    
    template< typename T>
    class JobQueue {
        

        typedef size_t jobid_t;
        typedef size_t groupid_t;

        jobid_t job_id;
        groupid_t group_id;

        std::map< designVars_t, jobid_t > JobList;
        std::map< jobid_t, objVars_t > Complete;
            
        std::map< jobid_t, std::vector<groupid_t> > JobsToGroups;
        std::map< groupid_t, int > GroupCounts;
        std::map< groupid_t, int > GroupToInd;
            
        void pushEvaluation(const designVars_t &x, jobid_t j){
            RemoteEvaluator.dispatch( x, j );
        }

        //Adjust the outstanding jobs
        void AcceptJob( jobid_t j, objVars_t &objs ){
            Complete[ j ] = objs;
            std::vector< groupid_t > & g (this->JobsToGroups[ j ]);
            std::vector< groupid_t >::iterator i;
            for(i=g.begin(); i!=g.end(); i++){
                this->GroupCounts[ (*i) ]-=1;
            }
        }

        //Get a completed group to work on
        groupid_t PopGroup(){
            std::map< groupid_t, int >::iterator i;
            for( i = GroupCounts.begin(); i != GroupCounts.end(); i++){
                if( i->second == 0 ) {
                    groupid_t retval = i->first;
                    GroupCounts.erase( i );
                    return retval;
                }
            }
            return 0;
        }

        //Get next index
        int NextInd(){
            groupid_t next_g = PopGroup();
            if( next_g ) {
                int ind = GroupToInd[ next_g ];
                GroupToInd.erase( next_g );
                return ind;
            }
            return -1;
        }

        //Request an evaluation
        objVars_t Request(const designVars_t &x) {
            //Have we requested this job yet?
            std::map< designVars_t, jobid_t >::iterator j = JobList.find( x );
            if( j != JobList.end() ){
                    
                jobid_t this_job = j->second;

                //If so, is it done?
                std::map< jobid_t, objVars_t >::iterator c = Complete.find( this_job );
                if( c != Complete.end() ) {
                    //Get the result
                    objVars_t result( c->second );

                    //Delete job
                    JobList.erase( j );
                    Complete.erase( c );
                    JobsToGroups.erase( this_job );

                    return result;
                }
                else {
                    //If not done, add this job to the group
                    JobsToGroups[ this_job ].push_back( this->group_id );
                    GroupCounts[ this->group_id ] += 1;
                }
            }
            else {
                //Job needs to be requested
                this->job_id++;

                //register with groups list
                std::vector< groupid_t > v;
                v.push_back( this->group_id );
                JobsToGroups[ this->job_id ] = v;
                //register with group counts
                GroupCounts[ this->group_id ] += 1;
                //register with jobs list
                JobList[ x ] = this->job_id;
                //request the evaluation
                this->pushEvaluation( x, this->job_id );
            }

            //return null result
            objVars_t empty_result;
            return empty_result;
        }
 
    public:
        T& RemoteEvaluator;
        typedef JobQueue<T>& BaseType;

        JobQueue( T & c ) : RemoteEvaluator( c ), job_id( 0 ), group_id( 1 ) {}
        //JobQueue( typename T::BaseType& c ) : RemoteEvaluator( c ), job_id( 0 ), group_id( 1 ) {}
        JobQueue( ) : job_id( 0 ), group_id( 1 ) {}

        void NewGroup(int ind){
            this->group_id++;
            GroupToInd[ this->group_id ] = ind;
            GroupCounts[ this->group_id ] = 0;
        }

        objVars_t eval(const designVars_t &x ){ return this->Request( x ); }
            
        int Poll() { 

            int nextInd = -1;
            while( nextInd < 0 ){
                std::queue< std::pair< jobid_t, objVars_t> > results;
                RemoteEvaluator.collect( results );

                while( !results.empty() ){
                    std::pair< jobid_t, objVars_t> &p( results.front() );
                    this->AcceptJob( p.first, p.second );
                    results.pop();
                }

                nextInd = this->NextInd();
            }

            return nextInd;
        }

        void Clear() {
            job_id = 0;
            group_id = 0;

            JobList.clear();
            Complete.clear();
            
            JobsToGroups.clear();
            GroupCounts.clear();
            GroupToInd.clear();
        }

        size_t Size(){ return job_id; }
    };

}
#endif
