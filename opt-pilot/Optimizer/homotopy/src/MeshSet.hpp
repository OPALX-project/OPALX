#ifndef MeshSet_hpp
#define MeshSet_hpp

#include "HomotopyTypes.h"
#include <map>
#include <queue>
#include <algorithm>
#include "mesh.hpp"

namespace Homotopy {
    
    typedef Mesh::ind_t ind_t;
    typedef Mesh::MeshBase::MeshPoint Point_t;
    /*
    typedef std::pair< ind_t, objVars_t > nodeEnvelope_t;
    typedef std::pair< ind_t, std::vector< nodeEnvelope_t > > NeighborMessage_t;
    */

    template< typename T>
    class GhostManager {

        T& Exchanger;
        std::vector< ind_t > IDs;
        std::map< ind_t, std::vector< Point_t* > > GlobalToGhosts;
        std::map< ind_t, std::vector< Point_t* > > SubsToLocals;
 
    public:
        GhostManager( T& exchanger ) : Exchanger( exchanger ) {}

        void Add( ind_t ID, std::map< ind_t, Point_t* > &globaltoghost, std::map< ind_t, std::vector< Point_t* > > &neighborsubsetstolocals ) {
            this->IDs.push_back( ID );

            std::map< ind_t, Point_t* >::iterator gg;
            for(gg=globaltoghost.begin(); gg!=globaltoghost.end(); gg++){
                if( this->GlobalToGhosts.find( gg->first ) == this->GlobalToGhosts.end() ) {
                    this->GlobalToGhosts[ gg->first ] = std::vector< Point_t* > ();
                }
                this->GlobalToGhosts[ gg->first ].push_back( gg->second );
            }
            std::map< ind_t, std::vector< Point_t* > >::iterator nstl;
            for(nstl=neighborsubsetstolocals.begin(); nstl!=neighborsubsetstolocals.end(); nstl++){
                if( this->SubsToLocals.find( nstl->first ) == this->SubsToLocals.end() ) {
                    this->SubsToLocals[ nstl->first ] = std::vector< Point_t* > ();
                }
                this->SubsToLocals[ nstl->first ].insert( this->SubsToLocals[ nstl->first ].end(), nstl->second.begin(), nstl->second.end() );
            }
        }          

        void Receive() {
            //Empty send queue
            std::queue< NeighborMessage_t > ToSend;
            //Container for results
            std::queue< nodeEnvelope_t > newNeighbors;

            Exchanger.exchange( ToSend, newNeighbors );

            //Process received node updates
            while( !newNeighbors.empty() ){
                nodeEnvelope_t &neighbor = newNeighbors.front();

                //Do we have this registered as a ghost
                std::map< ind_t, std::vector< Point_t* > >::iterator ghost_pos = GlobalToGhosts.find( neighbor.first );
                if( ghost_pos != GlobalToGhosts.end() ) {
                    std::vector< Point_t* >::iterator it;
                    for( it=ghost_pos->second.begin(); it!=ghost_pos->second.end(); it++ ){
                        (*it)->ObjectiveCoords = neighbor.second;
                    }
                }
                else {
                    //Something is wrong...
                    printf("GhostManager received erroneous node!\n");
                }
                newNeighbors.pop();
            }
        }

        void Exchange() {
            //Container of send requests format: pair( To, vector Content )
            //Complicated, but should be automatically serializable...
            std::queue< NeighborMessage_t > ToSend;

            //iterate over the neighbors
            std::map< ind_t, std::vector< Point_t* > >::iterator n;
            for(n=SubsToLocals.begin(); n!=SubsToLocals.end(); n++){
                
                //If the node is not managed by this instance,
                //then prepare a message to send
                if( std::find( this->IDs.begin(), this->IDs.end(), n->first ) == this->IDs.end() ){
                    std::vector< nodeEnvelope_t > toNode;
                    //gather the points to send to these neighbors
                    std::vector< Point_t* >::iterator v;
                    for(v=n->second.begin(); v!=n->second.end(); v++){
                        nodeEnvelope_t node( (int) ((*v)->ID), (*v)->ObjectiveCoords );
                        toNode.push_back( node );
                    }
                    ToSend.push( std::pair< ind_t, std::vector< nodeEnvelope_t > > ( n->first, toNode ) );
                }
                //Process locally managed updates
                else {
                    std::vector< Point_t* >::iterator p;
                    for(p=n->second.begin(); p!=n->second.end(); p++){
                        //Maybe we should check this?
                        std::vector< Point_t* >& v = this->GlobalToGhosts[ (*p)->ID ];

                        std::vector< Point_t* >::iterator it;
                        for( it=v.begin(); it!=v.end(); it++ ){
                            (*it)->ObjectiveCoords = (*p)->ObjectiveCoords;
                        }
                    }
                }
            }

            //Container for results
            std::queue< nodeEnvelope_t > newNeighbors;
            Exchanger.exchange( ToSend, newNeighbors );

            //Process received node updates
            while( !newNeighbors.empty() ){
                nodeEnvelope_t &neighbor = newNeighbors.front();

                //Do we have this registered as a ghost
                std::map< ind_t, std::vector< Point_t* > >::iterator ghost_pos = GlobalToGhosts.find( neighbor.first );
                if( ghost_pos != GlobalToGhosts.end() ) {
                    std::vector< Point_t* >::iterator it;
                    for( it=ghost_pos->second.begin(); it!=ghost_pos->second.end(); it++ ){
                        (*it)->ObjectiveCoords = neighbor.second;
                    }
                }
                else {
                    //Something is wrong...
                    printf("GhostManager received erroneous node!\n");
                }
                newNeighbors.pop();
            }
        }

        void Print() {
            printf("Managed IDs: ");
            std::vector< ind_t >::iterator id;
            for(id=IDs.begin(); id!=IDs.end(); id++) printf("%d ", *id);

            printf("\nIncoming Ghosts: ");
            std::map< ind_t, std::vector< Point_t* > >::iterator in;
            for(in=GlobalToGhosts.begin(); in!=GlobalToGhosts.end(); in++) printf("%d ", in->first );
            
            printf("\nOutgoing Meshes: ");
            std::map< ind_t, std::vector< Point_t* > >::iterator out;
            for(out=SubsToLocals.begin(); out!=SubsToLocals.end(); out++) printf("%d ", out->first );
            printf("\n");
        }
    };


    template< typename T, typename E >
    class MeshSet : public Mesh::MeshBase {
    private:
        //TODO: Could encapsulate various T constructor args in struct...
        std::vector< T* > meshes;
        GhostManager< E > gman;
    public:
        MeshSet(    std::vector< Mesh::point_t > DesignSpace,
                    std::vector< Mesh::point_t > ObjectiveSpace,
					std::vector< Mesh::point_t > Lambdas, 
					int NumberOfPoints, 
                    std::vector<Mesh::ind_t> IDs, 
                    Mesh::ind_t SubsetsPerSide,
                    E & exchanger) : 
                                            MeshBase( DesignSpace, ObjectiveSpace, Lambdas ),
                                            gman( exchanger ) {
            std::vector< ind_t >::iterator i;
            for(i=IDs.begin(); i!=IDs.end(); i++){
                T* m;
                try{ 
                    m = new T(    DesignSpace, ObjectiveSpace, Lambdas, 
                                            NumberOfPoints, *i, SubsetsPerSide ) ;
                }
                catch ( std::bad_alloc& e)
                {
                    m = NULL;
                    printf("Error allocating memory for mesh %d\n", *i );
                }

                meshes.push_back( m );
            }
        }

        virtual size_t GetSize() { return meshes.size(); }
        virtual T* Get( ind_t id ){ return meshes[id]; }

        virtual void Generate() {
            typename std::vector< T* >::iterator i;
            for(i=meshes.begin(); i!=meshes.end(); i++){
                //Generate all the meshes
                (*i)->Generate();
                //Place them under local management
                gman.Add( (*i)->ID, (*i)->GlobalToGhost, (*i)->NeighborSubsetsToLocals );
            }
            //Send local nodes to remote locations, and 
            //receive local/remote ghost nodes
            this->gman.Exchange();

            //Make sure to recieve all remote ghost nodes before proceeding
            bool valid = false;
            while( !valid ){
                valid = true;
                for(i=meshes.begin(); i!=meshes.end(); i++){
                    //Probe mesh validity
                    valid &= (*i)->Valid();
                }

                this->gman.Receive();
            }
        }
        virtual void Refresh() {
            typename std::vector< T* >::iterator i;
            for(i=meshes.begin(); i!=meshes.end(); i++){
                (*i)->Refresh();
            }
            this->gman.Exchange();
        }
        virtual void Print() {
            typename std::vector< T* >::iterator i;
            for(i=meshes.begin(); i!=meshes.end(); i++){
                (*i)->Print();
            }
            this->gman.Print();
        }
        virtual void WriteOut() {
            typename std::vector< T* >::iterator i;
            for(i=meshes.begin(); i!=meshes.end(); i++){
                (*i)->WriteOut();
            }
        }
    };
}

#endif
