#include "gtest/gtest.h"

#include <vector>
#include <set>
#include "mesh.hpp"
#include <cstdlib>

#include <stdio.h>

namespace {

	using namespace std;

	template< typename T >
	void PrintV(std::vector< T > v){
		int i;
		for(i=0; i<v.size(); i++) { cout << v[i] << "\t"; }
		cout << std::endl;
	}

	template< typename T >
	std::vector< T > VecDiff(std::vector< T > v1, std::vector< T > v2){ 
		std::vector< T > res(v1.size());
		int i;
		for(i=0; i<v1.size(); i++) { res[i] = v1[i] - v2[i]; }
		return res;
	}
	template< typename T >
	T VecDot(std::vector< T > v1, std::vector< T > v2){ 
		T res = 0.0;
		int i;
		for(i=0; i<v1.size(); i++) { res += v1[i] * v2[i]; }
		return res;
	}

	//Verify that the coordinate transforms are
	//stable and reversible
	TEST(SimplexTest, CoordinateTransforms){
		int dim = 8;
		int n = 8;

		int i;
		for(i = 0; i < Mesh::Simplex::eta(dim, n); i++){
			std::vector< int > v( Mesh::Simplex::ind_to_coord( i, dim, n) ) ;
			EXPECT_EQ( i,  Mesh::Simplex::coord_to_ind( v , n ) );
		}
	}

	//Verify that the selected neighbors are antiparallel
	TEST(SimplexTest, Neighbors){
		int dim = 8;
		int n = 8;

		//Setup the corners
		srand( 0 );
		std::vector< std::vector < double > > D;
		int i;
		for(i = 0; i < dim; i++){
			std::vector< double > v;
			int j;
			for(j = 0; j < dim; j++){ v.push_back( (double) (rand() % 1000) / 1000.0 ); }
			D.push_back( v );
		}

		//Build a Mesh
		Mesh::Simplex mesh( D, D, D, n);
        mesh.Generate();

		//For each meshpoint
		for(i=0; i<mesh.Points.size(); i++){
			typename Mesh::MeshBase::MeshPoint *p = &mesh.Points[i];
			int j;
			for(j=0; j<p->Neighbors.size(); j+=2){
				//Find the various vector differences
				double dot = VecDot(	VecDiff(p->ObjectiveCoords, mesh.Points[p->Neighbors[j]].ObjectiveCoords), 
										VecDiff(p->ObjectiveCoords, mesh.Points[p->Neighbors[j+1]].ObjectiveCoords) ); 

				double mag1 = VecDot(	VecDiff(p->ObjectiveCoords, mesh.Points[p->Neighbors[j]].ObjectiveCoords), 
										VecDiff(p->ObjectiveCoords, mesh.Points[p->Neighbors[j]].ObjectiveCoords) ) ;
				//Enforce antiparallelism
				EXPECT_NEAR(-1.0, dot / mag1, 1e-6);

				double mag2 = VecDot(	VecDiff(p->ObjectiveCoords, mesh.Points[p->Neighbors[j+1]].ObjectiveCoords), 
										VecDiff(p->ObjectiveCoords, mesh.Points[p->Neighbors[j+1]].ObjectiveCoords) ) ;
				//Enforce equidistance
				EXPECT_NEAR(0.0, mag1 - mag2 , 1e-6);
			}
		}
	}

    //Verify neighborhood access
	TEST(SimplexTest, Neighborhood){
		int dim = 3;
		int n = 3;

		int i;
		for(i = 0; i < Mesh::Simplex::eta(dim, n); i++){
			std::vector< int > v( Mesh::Simplex::ind_to_coord( i, dim, n) ) ;
			PrintV< int >( Mesh::Simplex::getNeighborhood( v , n ) );
		}
    }

    //Verify that the selected neighbors are antiparallel
	TEST(SimplexTest, DISABLED_NeighborAccess){
		int dim = 8;
		int n = 8;

		//Setup the corners
		srand( 0 );
		std::vector< std::vector < double > > D;
		int i;
		for(i = 0; i < dim; i++){
			std::vector< double > v;
			int j;
			for(j = 0; j < dim; j++){ v.push_back( (double) (rand() % 1000) / 1000.0 ); }
			D.push_back( v );
		}

		//Build a Mesh
		Mesh::Simplex mesh( D, D, D, n);
        mesh.Generate();

        //For each meshpoint
		for(i=0; i<mesh.Points.size(); i++){
            //Get the set
            std::set< int > n ( mesh.GetNeighborIDs( i ) );
            //Get the neighbors
            std::vector< int > &vec (mesh.Points[i].Neighbors);

            //Same size?
            EXPECT_EQ( vec.size(), n.size() );

            //Be sure that every element is present
            std::vector< int >::iterator it;
            for(it=vec.begin(); it!=vec.end(); it++) {
                bool elem = ( n.find( *it ) != n.end() );
                EXPECT_EQ( elem, true );
            }
        }
    }

	//Build a 1D mesh (line)
	TEST(SimplexTest, 1D){

		std::vector< std::vector < double > > v, L;
		double vals[3][3] = { {1, 0}, {0, 1} };
		double *a;

		a = &vals[0][0];
		std::vector< double > p;
		p.assign(a, a+2);
		v.push_back(p);

		a = &vals[1][0];
		p.assign(a, a+2);
		v.push_back(p);

		Mesh::Simplex mesh( v, v, v, 5);
        mesh.Generate();
		mesh.Print();
        mesh.WriteOut( "1d.dat" );
	}

	//Build a 2D mesh (trianglar plane)
	TEST(SimplexTest, 2D){

		std::vector< std::vector < double > > D, O, L;
		double Ovals[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 0} };
		double Lvals[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
		double *a;
		double *b;
		int dim = 3;
		
		std::vector< double > p;
		int i;
		for(i=0;i<dim;i++){
			a = &Ovals[i][0];
			p.assign(a, a+dim);
			D.push_back(p);

			p.assign(a, a+dim);
			O.push_back(p);

			a = &Lvals[i][0];
			p.assign(a, a+dim);
			L.push_back(p);
		}

		Mesh::Simplex mesh( D, O, L, 3);
        mesh.Generate();
		mesh.Print();
        mesh.WriteOut( "2d.dat" );
	}

	//Build a 3D mesh (quarter pyramid)
	TEST(SimplexTest, 3D){

		std::vector< std::vector < double > > D, O, L;
		double Ovals[4][4] = { {0, 0, 0, 0}, {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0} };
		double Lvals[4][4] = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
		double *a;
		double *b;
		int dim = 4;
		
		std::vector< double > p;
		int i;
		for(i=0;i<dim;i++){
			a = &Ovals[i][0];
			p.assign(a, a+dim);
			D.push_back(p);

			p.assign(a, a+dim);
			O.push_back(p);

			a = &Lvals[i][0];
			p.assign(a, a+dim);
			L.push_back(p);
		}

		Mesh::Simplex mesh( D, O, L, 3);
        mesh.Generate();
		mesh.Print();
        mesh.WriteOut( "3d.dat" );
	}

    TEST(SimplexNodeSetTest, One){
        int dim = 8;
		int n = 8;

        int i;
        for(i=0; i<Mesh::Simplex::eta(dim, n); i++){
            Mesh::ind_t s = Mesh::SimplexNodeSet::WhichSubset( Mesh::Simplex::ind_to_coord( i, dim, n) , n, 1 );
            EXPECT_EQ( s, 0 );
        }
    }

    TEST(SimplexNodeSetTest, All){
        int dim = 8;
		int n = 8;

        int i;
        for(i=0; i<Mesh::Simplex::eta(dim, n); i++){
            Mesh::ind_t s = Mesh::SimplexNodeSet::WhichSubset( Mesh::Simplex::ind_to_coord( i, dim, n) , n, n );
            EXPECT_EQ( s, i );
        }
    }

    TEST(SimplexSubsetTest, Alive){
        int dim = 3;
        int subsetsperside = 3;
        int pointsperside = 3;

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

        for(i=0; i<(Mesh::SimplexNodeSet::NumberOfSubsets(dim-1, subsetsperside)); i++) {
            Mesh::SimplexSubset mesh( D, D, D, pointsperside, i, subsetsperside );
            mesh.Generate();
            mesh.Print();
        }
    }

        
    
}