//Simple test to verify interpolator functionality

#include "gtest/gtest.h"

#include <vector>

#include "interpolator.hpp"

namespace {
	TEST(InterpolatorTest, Alive) {
		std::vector< std::vector< double > > data;
		std::vector< double > v(3);
		v[0] = -1.0; v[1] = 0.0; v[2] = 2.0;
		data.push_back( v );
		v[0] = 1.0; v[1] = 0.0; v[2] = 3.0;
		data.push_back( v );

		Interpolation::RBF rbf(Interpolation::RBF::UNSCALED, data);
		std::vector< double > point(2, 0.0);
		double result = rbf.evaluate( point );
		printf("%lf\n", result);
		EXPECT_NEAR( 2.5, double(result), 1e-3 );
	}

    TEST(InterpolatorTest, FromFile) {
		std::vector< std::vector< double > > data = Interpolation::GetDataFromFile("./data/rms_s.dat") ;

		Interpolation::RBF rbf(Interpolation::RBF::UNSCALED, data);

        std::vector< double > point(2, 0.0);
        double sum = 0.0;

        std::vector< std::vector< double > >::iterator it;
        for(it=data.begin(); it != data.end(); it++) {
            point[0] = (*it)[0]; point[1] = (*it)[1];
            double result = rbf.evaluate( point ) - (*it)[2];
		    sum += result*result;
        }

        //Total sum squared error
        EXPECT_NEAR( 0.0, sum, 1e-3 );
	}

    TEST(InterpolatorTest, DISABLED_VerifyOutput) {
		std::vector< std::vector< double > > data = Interpolation::GetDataFromFile("./data/rms_s.dat") ;

		Interpolation::RBF rbf(Interpolation::RBF::RESCALED, data);

        std::vector< double > point(2, 0.0);
        double sum = 0.0;

        int N = 200;
       
        double hix = 1.0;
        double lox = 0.0;

        double hiy = 1.0;
        double loy = 0.0;


        double xstep = (hix - lox)/(double)(N-1);
        double ystep = (hiy - loy)/(double)(N-1);

        double x=lox;
        int i;
        for(i=0; i<N; i++, x+=xstep) {
            double y=loy;
            int j;
            for(j=0; j<N; j++,y+=ystep) {
                point[0] = x; point[1] = y;
                double result = rbf.evaluate( point );
                std::cout << x << "\t" << y << "\t" << result << std::endl;
            }
        }
	}
}