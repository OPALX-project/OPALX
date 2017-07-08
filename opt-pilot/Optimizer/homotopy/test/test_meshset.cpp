//Validate the evaluator

#include "gtest/gtest.h"

#include "boost/bind.hpp"
#include "boost/function.hpp"

#include "HomotopyTypes.h"
#include "mesh.hpp"
#include "MeshSet.hpp"

#include <stdio.h>

using namespace Homotopy;

namespace {

    void Print(const std::vector< double > & x ){
        int i;
	    for(i=0; i<x.size(); i++) { 
		    printf("%lf ", x[i]);
	    }
	    printf("\n");
    }

    class TestExchanger{
    public:
        ind_t rec_ind;
        objVars_t obj;

        ind_t last_to;
        ind_t last_addr;
        objVars_t last_msg;

        TestExchanger( ind_t ind ) : rec_ind( ind ) {
            this->obj = objVars_t (1, 0.0);
        }
        void exchange(  std::queue< NeighborMessage_t > & toSend, std::queue< nodeEnvelope_t > & received ){
            //Mock send
            while( !toSend.empty() ){
                NeighborMessage_t& m = toSend.front();
                /*
                printf("To %d: \n", m.first );
                std::vector< nodeEnvelope_t >::iterator i;
                for( i=m.second.begin(); i!=m.second.end(); i++){
                    Print( i->second );
                }
                */
                this->last_to = m.first;
                this->last_addr = m.second[ 0 ].first;
                this->last_msg = m.second[ 0].second;
                toSend.pop();
            }

            //Mock receive
            if( this->rec_ind >= 0){
                nodeEnvelope_t e ( this->rec_ind, this->obj );
                received.push( e );
            }
        }
    };
    

	TEST(GHOSTMANAGER, Alive){
        TestExchanger te( 1 );

        std::vector< double > v(3, 1.0);
        Point_t point0( 0, v, v, v);
        Point_t point1( 1, v, v, v);

        std::vector< Point_t* > pointvec;
        pointvec.push_back( &point0 );

        std::map< ind_t, std::vector< Point_t* > > ns;
        ns[0] = pointvec;

        std::map< ind_t, Point_t* > gg;
        gg[1] = &point1;

        GhostManager< TestExchanger > gm( te );
        gm.Add( 1, gg, ns );
        gm.Exchange();

        EXPECT_EQ( point0.ObjectiveCoords, te.last_msg );
        EXPECT_EQ( point1.ObjectiveCoords, te.obj );
	}

    TEST(MESHSET, Alive){
        //Test parameters
        int dim = 3;
        int subsetsperside = 2;
        int pointsperside = 4;

        //Setup the corners
		srand( 0 );
		std::vector< std::vector < double > > D;
		int i;
		for(i = 0; i < dim; i++){
			std::vector< double > v;
			int j;
			for(j = 0; j < dim; j++){ v.push_back( (i==j ? 1.0 : 0.0) ); }
			D.push_back( v );
		}

        std::vector< ind_t > IDs;
        IDs.push_back( 0 );
        IDs.push_back( 1 );
        IDs.push_back( 2 );

        TestExchanger te( -5 );

        MeshSet< Mesh::SimplexSubset, TestExchanger > m( D, D, D, pointsperside, IDs , subsetsperside, te);
        m.Generate();
        m.Refresh();
        m.Print();
    }

    TEST(MESHSET, AllGhosts){
        //Test parameters
        int dim = 2;
        int pointsperside = 100;

        //Setup the corners
		srand( 0 );
		std::vector< std::vector < double > > D;
		int i;
		for(i = 0; i < dim; i++){
			std::vector< double > v;
			int j;
			for(j = 0; j < dim; j++){ v.push_back( (i==j ? 1.0 : 0.0) ); }
			D.push_back( v );
		}

        std::vector< ind_t > IDs;
        for(i=0; i<Mesh::Simplex::eta( dim-1, pointsperside ); i++){
            IDs.push_back( i );
        }

        TestExchanger te( -1 );

        MeshSet< Mesh::SimplexSubset, TestExchanger > m( D, D, D, pointsperside, IDs , pointsperside, te);
        m.Generate();
        m.Refresh();
        m.Refresh();

        for(i=0; i<IDs.size(); i++){
            EXPECT_EQ( m.Get(i)->Points.size(), 1 );
            EXPECT_EQ( m.Get(i)->GhostNodes.size(), m.Get(i)->Points[0].NeighborP.size() );
            int j;
            for(j=0; j<m.Get(i)->GhostNodes.size(); j++){
                Mesh::MeshBase::MeshPoint* g = m.Get(i)->GhostNodes[j];
                EXPECT_EQ(  m.Get(g->ID)->Points[0].ObjectiveCoords, g->ObjectiveCoords );
            }
        }
    }

    TEST(SIMPLEXNODESET, Mapping){
        int dim = 2;

        int PointsPerSide = 5;
        int SubsetsPerSide = 4;

        int i;
        for(i=0; i<Mesh::Simplex::eta( dim, PointsPerSide ); i++){

            std::vector< int > coords = Mesh::Simplex::ind_to_coord( i, dim, PointsPerSide);

            printf( "%d: %d\n", i, Mesh::SimplexNodeSet::WhichSubset( coords, PointsPerSide, SubsetsPerSide ) );
        }
    }
}