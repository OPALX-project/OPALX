//Simple test to verify AlgLib functionality

#include "gtest/gtest.h"

#include <math.h>
#include <cstdio>
#include <vector>

#include "stdafx.h"
#include "interpolation.h"

namespace {

using namespace alglib;

// Tests that the AlgLib interpolator works
 TEST(AlgLibTest, QNN) {
    // This example illustrates basic concepts of the RBF models: creation, modification,
    // evaluation.
    // 
    // Suppose that we have set of 2-dimensional points with associated
    // scalar function values, and we want to build a RBF model using
    // our data.
    // 
    // NOTE: we can work with 3D models too :)
    // 
    // Typical sequence of steps is given below:
    // 1. we create RBF model object
    // 2. we attach our dataset to the RBF model and tune algorithm settings
    // 3. we rebuild RBF model using QNN algorithm on new data
    // 4. we use RBF model (evaluate, serialize, etc.)
    //
    double v;

    //
    // Step 1: RBF model creation.
    //
    // We have to specify dimensionality of the space (2 or 3) and
    // dimensionality of the function (scalar or vector).
    //
    rbfmodel model;
    rbfcreate(2, 1, model);

    // New model is empty - it can be evaluated,
    // but we just get zero value at any point.
    v = rbfcalc2(model, 0.0, 0.0);
    printf("%.2f\n", double(v)); // EXPECTED: 0.000
	EXPECT_NEAR( 0.0, double(v), 1e-3 );

    //
    // Step 2: we add dataset.
    //
    // XY arrays containt two points - x0=(-1,0) and x1=(+1,0) -
    // and two function values f(x0)=2, f(x1)=3.
    //
    real_2d_array xy = "[[-1,0,2],[+1,0,3]]";
    rbfsetpoints(model, xy);

    // We added points, but model was not rebuild yet.
    // If we call rbfcalc2(), we still will get 0.0 as result.
    v = rbfcalc2(model, 0.0, 0.0);
    printf("%.2f\n", double(v)); // EXPECTED: 0.000
	EXPECT_NEAR( 0.0, double(v), 1e-3 );

    //
    // Step 3: rebuild model
    //
    // After we've configured model, we should rebuild it -
    // it will change coefficients stored internally in the
    // rbfmodel structure.
    //
    // By default, RBF uses QNN algorithm, which works well with
    // relatively uniform datasets (all points are well separated,
    // average distance is approximately same for all points).
    // This default algorithm is perfectly suited for our simple
    // made up data.
    //
    // NOTE: we recommend you to take a look at example of RBF-ML,
    // multilayer RBF algorithm, which sometimes is a better
    // option than QNN.
    //
    rbfreport rep;
    rbfsetalgoqnn(model);
    rbfbuildmodel(model, rep);
    printf("%d\n", int(rep.terminationtype)); // EXPECTED: 1
	EXPECT_NEAR( 1, int(rep.terminationtype), 0 );

    //
    // Step 4: model was built
    //
    // After call of rbfbuildmodel(), rbfcalc2() will return
    // value of the new model.
    //
    v = rbfcalc2(model, 0.0, 0.0);
    printf("%.2f\n", double(v)); // EXPECTED: 2.500
	EXPECT_NEAR( 2.5, double(v), 1e-3 );
	}

 TEST(AlgLibTest, ML) {
	//
    // This example shows how to build models with RBF-ML algorithm. Below
    // we assume that you already know basic concepts shown in the example
    // on RBF-QNN algorithm.
    //
    // RBF-ML is a multilayer RBF algorithm, which fits a sequence of models
    // with decreasing radii. Each model is fitted with fixed number of
    // iterations of linear solver. First layers give only inexact approximation
    // of the target function, because RBF problems with large radii are
    // ill-conditioned. However, as we add more and more layers with smaller
    // and smaller radii, we get better conditioned systems - and more precise models.
    //
    rbfmodel model;
    rbfreport rep;
    double v;

    //
    // We have 2-dimensional space and very simple interpolation problem - all
    // points are distinct and located at straight line. We want to solve it
    // with RBF-ML algorithm. This problem is very simple, and RBF-QNN will
    // solve it too, but we want to evaluate RBF-ML and to start from the simple
    // problem.
    //     X        Y
    //     -2       1
    //     -1       0
    //      0       1
    //     +1      -1
    //     +2       1
    //
    rbfcreate(2, 1, model);
    real_2d_array xy0 = "[[-2,0,1],[-1,0,0],[0,0,1],[+1,0,-1],[+2,0,1]]";
    rbfsetpoints(model, xy0);

    // First, we try to use R=5.0 with single layer (NLayers=1) and moderate amount
    // of regularization.... but results are disappointing: Model(x=0,y=0)=-0.02,
    // and we need 1.0 at (x,y)=(0,0). Why?
    //
    // Because first layer gives very smooth and imprecise approximation of the
    // function. Average distance between points is 1.0, and R=5.0 is too large
    // to give us flexible model. It can give smoothness, but can't give precision.
    // So we need more layers with smaller radii.
    rbfsetalgomultilayer(model, 5.0, 1, 1.0e-3);
    rbfbuildmodel(model, rep);
    
	printf("%d\n", int(rep.terminationtype)); // EXPECTED: 1
	EXPECT_NEAR( 1, int(rep.terminationtype), 0 );

    v = rbfcalc2(model, 0.0, 0.0);
    printf("%.2f\n", double(v)); // EXPECTED: -0.021690
	EXPECT_NEAR( -0.021690, double(v), 1e-6 );

    // Now we know that single layer is not enough. We still want to start with
    // R=5.0 because it has good smoothness properties, but we will add more layers,
    // each with R[i+1]=R[i]/2. We think that 4 layers is enough, because last layer
    // will have R = 5.0/2^3 = 5/8 ~ 0.63, which is smaller than the average distance
    // between points. And it works!
    rbfsetalgomultilayer(model, 5.0, 4, 1.0e-3);
    rbfbuildmodel(model, rep);
    printf("%d\n", int(rep.terminationtype)); // EXPECTED: 1
	EXPECT_NEAR( 1, int(rep.terminationtype), 0 );

    v = rbfcalc2(model, 0.0, 0.0);

    printf("%.2f\n", double(v)); // EXPECTED: 1.000000
	EXPECT_NEAR( 1.0, double(v), 1e-6 );

	}

}  //namespace
