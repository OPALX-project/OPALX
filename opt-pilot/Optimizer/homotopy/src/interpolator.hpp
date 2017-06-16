#ifndef interpolator_hpp
#define interpolator_hpp

#include "stdafx.h"
#include "interpolation.h"

#include <vector>
#include <algorithm>
#include <stdio.h>

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>

#include <cmath>
#include <cstdlib>
#include "Eigen/Dense"
#include "Eigen/LU"
using namespace Eigen;

namespace Interpolation {
	class Interpolator { };

	typedef std::vector< std::vector< double > > Data;

	class RBF_AlgLib {
	private:
		int Points;
		int dim;
		
		std::vector< double > mins;
		std::vector< double > ranges;

		alglib::rbfmodel model;
		alglib::rbfreport report;

    public:
		template< typename T >
		class vertical_iterator : public T::const_iterator {
			int column;
		public:
			vertical_iterator( int c ) : column( c ) {}
			typename T::value_type::value_type operator*() const { return (this->T::const_iterator::operator*())[column]; }
			vertical_iterator& operator= (const typename T::const_iterator & rhs)  { this->T::const_iterator::operator=(rhs); return *this; }
		};

		static std::vector< double > FindMins(const Data &d ){
			//Make sure we have some data
			assert ( !d.empty() );
			assert ( !d.front().empty() );

			std::vector< double > result( d[0].size() );
			int i;
			for(i=0; i<result.size(); i++) { 
				vertical_iterator< Data > vi1(i);
				vi1 = d.begin();
				vertical_iterator< Data > vi2(i);
				vi2 = d.end();
				result[i] = *std::min_element( vi1, vi2 ); 
			}
			return result;
		}
		//Maybe we can eliminate dimensions that have no range...
		static std::vector< double > FindRanges(const Data &d ){
			std::vector< double > result( d[0].size() );
			int i;
			for(i=0; i<result.size(); i++) { 
				vertical_iterator< Data > vi1(i);
				vertical_iterator< Data > vi2(i);
				vi1 = d.begin();
				vi2 = d.end();
				result[i] = *std::max_element( vi1, vi2 ) - *std::min_element( vi1, vi2 ); 
				//assert( result[i] != 0 );
			}
			return result;
		}
		inline static double Rescale( double x, double min, double range) { return (range == 0) ? x : ((x-min)/range); }
		inline static double Unscale( double x, double min, double range) { return (range == 0) ? x : ((x*range)+min); }
		double* PrepareData( const Data &data, int N, int dim ) {
			double* d;
			try {
				d = new double[N*dim];
			}
			//catch( std::bad_alloc& ) {
			catch( std::exception& ) {
				printf("Dynamic allocation failed while staging RBF data for copy\n");
				d = NULL;
				return d;
			}

			int i;
			for(i=0; i<N; i++) {
				int j;
				for(j=0; j<dim; j++) { d[i*dim + j] = Rescale( data[i][j], this->mins[j], this->ranges[j] );
                //if( j==dim-1 ) d[i*dim+j] =  1.0 / (1.0 + d[i*dim+j]);
                }
			}

			return d;
		}

	//public:
        enum SCALE {UNSCALED, RESCALED};
        SCALE scale;

		RBF_AlgLib( SCALE s, Data& data ) : mins( FindMins( data ) ), ranges( FindRanges( data ) ), scale(s) {
			try
			{
				//Assume one dependent variable
				this->dim = data.front().size() - 1;
				this->Points = data.size();

				assert( this->dim == 2 || this->dim == 3 );

				double *d = PrepareData( data, this->Points, this->dim + 1);

				alglib::real_2d_array vals;
				vals.setcontent(this->Points, this->dim+1, d);
				delete [] d;

				alglib::rbfcreate(this->dim, 1, this->model);
				alglib::rbfsetpoints(this->model, vals);
				//alglib::rbfsetalgoqnn(this->model, 2.1, 5.0);
                alglib::rbfsetalgomultilayer(this->model, 1.0, 8, 1.0e-2);
                //rbfsetlinterm( this->model );
                //alglib::rbfsetconstterm( this->model );
                //alglib::rbfsetzeroterm( this->model );
				alglib::rbfbuildmodel(this->model, this->report);
				assert( int(this->report.terminationtype) == 1 );
			}
			catch(alglib::ap_error e)
			{
				printf("error msg: %s\n", e.msg.c_str());
			}			
		}

		double evaluate(const std::vector< double > &x){
			assert( x.size() == this->dim );
			if( this->dim == 2 ) 
                switch( this->scale ) {
                    case RESCALED:

                    return alglib::rbfcalc2(this->model, 
				    x[0], 
				    x[1]);
                    break;

                    default: 

                    return Unscale( alglib::rbfcalc2(this->model, 
                    Rescale(x[0], this->mins[0], this->ranges[0] ), 
				    Rescale(x[1], this->mins[1], this->ranges[1] ) ), 
				    this->mins[2], this->ranges[2] ) ;
            }

			if( this->dim == 3 ) 
                switch( this->scale ) {
                    case RESCALED:
                        return alglib::rbfcalc3(this->model, 
				        x[0], 
				        x[1],
                        x[2]);
                        break;

                    default:
				        return Unscale( alglib::rbfcalc3(this->model, 
				        Rescale(x[0], this->mins[0], this->ranges[0] ), 
				        Rescale(x[1], this->mins[1], this->ranges[1] ),
				        Rescale(x[2], this->mins[2], this->ranges[2] ) ),
				        this->mins[3], this->ranges[3] ) ;
                
                }
		}
	};

    class MyRBF {
    private:
        int N;
		int dim;
		
		std::vector< double > mins;
		std::vector< double > ranges;

        friend class RBF_AlgLib;

        VectorXd weights;
        std::vector< VectorXd > points;

        inline static double Rescale( double x, double min, double range) { return (range == 0) ? x : ((x-min)/range); }
		inline static double Unscale( double x, double min, double range) { return (range == 0) ? x : ((x*range)+min); }

        double rho( double r ) { return r; }
        //double rho( double r ) { return exp( -r * r); }
        double dist( VectorXd& l, VectorXd& r ){ 
            return sqrt((l-r).dot(l-r));
        }
        void Setup( Data& data ) {
            VectorXd d(N);

            //store the points & values
            int i;
            for(i=0;i<N;i++){
                VectorXd tmp( dim );
                int j;
                for(j=0;j<dim;j++){
                    tmp[j] = Rescale( data[i][j], this->mins[j], this->ranges[j] );
                }
                this->points.push_back( tmp );
                d[i] = Rescale( data[i][dim], this->mins[dim], this->ranges[dim] );
            }

            //get the A matrix
            MatrixXd A( N, N );
            for(i=0;i<N;i++){
                int j;
                for(j=i;j<N;j++){
                    double r = rho( dist( points[i], points[j] ) );
                    A( i, j ) = r;
                    A( j, i ) = r;
                }
            }

            //Solve for the weights
            //this->weights = A.colPivHouseholderQr().solve( d );
            this->weights = A.fullPivLu().solve( d );
        }

    public:
        enum SCALE {UNSCALED, RESCALED};
        SCALE scale;

        MyRBF( SCALE s , Data& data ) : mins( RBF_AlgLib::FindMins( data ) ), ranges( RBF_AlgLib::FindRanges( data ) ), scale( s ) {
            //Assume one dependent variable
			this->dim = data.front().size() - 1;
			this->N = data.size();

			assert( this->dim == 2 || this->dim == 3 );

            this->Setup( data );
        }

        double evaluate(const std::vector< double > &x){
            Eigen::VectorXd p(this->dim);
            int n;
            for(n=0; n<this->dim; n++) p[n]=x[n];

            
            if( this->scale == UNSCALED ){ 
                int j;
                for(j=0; j<dim; j++) p[j] =  Rescale( p[j], this->mins[j], this->ranges[j] );
            }

            VectorXd r( N );
            int i;
            for(i=0; i<this->N; i++){ r[i] = rho( dist( this->points[i], p ) ); }

            if( this->scale == UNSCALED ) return Unscale( this->weights.dot( r ), this->mins.back(), this->ranges.back() );
            return this->weights.dot( r );
        }
    };

    //Which RBF to use...
    //typedef RBF_AlgLib RBF;
    typedef MyRBF RBF;

    Data GetDataFromFile(const char* filename){
        Data ret;
        std::string line;
        std::ifstream myfile (filename);
        
        if (myfile.is_open())
        {
            while ( myfile.good() )
            {
                std::getline (myfile,line);

                std::vector< double > v;
                double x, y, z;

                sscanf( line.c_str(), "%lf %lf %lf", &x, &y, &z );

                v.push_back(x); v.push_back(y); v.push_back(z);
                ret.push_back(v);
            }
            myfile.close();
        }

        else throw std::runtime_error("Unable to open data file!");

        return ret;
    }
}

#endif