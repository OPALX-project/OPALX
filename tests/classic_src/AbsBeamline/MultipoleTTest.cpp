#include "gtest/gtest.h"
#include "AbsBeamline/MultipoleT.h"
#include<fstream>

using namespace std;

vector< vector<double> > partialsDerivB(const Vector_t &R,const Vector_t B, double stepSize, MultipoleT* dummyField)
{  
    // builds a matrix of all partial derivatives of B -> dx_i B_j 
    vector< vector<double> > allPartials(3, vector<double>(3));
    double t = 0 ;
    Vector_t P, E;
    for(int i = 0; i < 3; i++)
	{ 
	  // B at the previous and next grid points R_prev,  R_next
	  Vector_t R_prev = R, R_next = R;
	  R_prev[i] -= stepSize;
	  R_next[i] += stepSize;
	  Vector_t B_prev, B_next;
	  dummyField->apply(R_prev, P, t, E, B_prev);
	  dummyField->apply(R_next, P, t, E, B_next); 
	  for(int j = 0; j < 3; j++)
	    allPartials[i][j] = (B_next[j] - B_prev[j]) / (2 * stepSize);
	}
     return allPartials;
}

vector< vector<double> > partialsDerivB_5(const Vector_t &R,const Vector_t B, double stepSize, MultipoleT* dummyField)
{  
    // builds a matrix of all partial derivatives of B -> dx_i B_j 
    vector< vector<double> > allPartials(3, vector<double>(3));
    double t = 0 ;
    Vector_t P, E;
    for(int i = 0; i < 3; i++)
	{ 
	  // B at the previous and next grid points R_prev,  R_next
	  Vector_t R_pprev = R, R_prev = R, R_next = R, R_nnext = R;
	  R_pprev(i) -= 2 * stepSize;
	  R_nnext(i) += 2 * stepSize;
	  R_prev(i) -= stepSize;
	  R_next(i) += stepSize;
	  Vector_t B_prev, B_next, B_pprev, B_nnext;
	  dummyField->apply(R_prev, P, t, E, B_prev);
	  dummyField->apply(R_next, P, t, E, B_next);
	  dummyField->apply(R_pprev, P, t, E, B_pprev);
	  dummyField->apply(R_nnext, P, t, E, B_nnext);
	  for(int j = 0; j < 3; j++)
	    allPartials[i][j] = (B_pprev[j] - 8 * B_prev[j] + 8 * B_next[j] - B_nnext[j]) / (12 * stepSize);
	}
     return allPartials;
}
    
double calcDivB(Vector_t &R, Vector_t B, double stepSize, MultipoleT* dummyField )
{
    double div = 0;
    vector< vector<double> > partials (3, vector<double>(3));
    partials = partialsDerivB(R, B, stepSize, dummyField);
    for(int i = 0; i < 3; i++)
        div += partials[i][i];
    return div;
}   

vector<double> calcCurlB(Vector_t &R, Vector_t B, double stepSize, MultipoleT* dummyField)
{
    vector<double> curl(3);
    vector< vector<double> > partials(3, vector<double>(3));
    partials = partialsDerivB(R, B, stepSize, dummyField);
    curl[0] = (partials[1][2] - partials[2][1]);
    curl[1] = (partials[2][0] - partials[0][2]);
    curl[2] = (partials[0][1] - partials[1][0]);
    return curl;
}

TEST(MultipoleT, Field)
{
    MultipoleT* myMagnet = new MultipoleT("Quadrupole");
    double centralField = 5;
    double fringeLength = 0.5;
    // the highest differential in the fringe field
    double max_index = 20;
    myMagnet->setFringeField(centralField, fringeLength, max_index);
    myMagnet->setTransMaxOrder(1);
    myMagnet->setDipoleConstant(0.0);
    myMagnet->setTransProfile(1, 100.0);
    cout<<"test: "<<myMagnet->getTransProfile(1)<<' '<<myMagnet->getTransMaxOrder()<<endl;
    //highest power in the field is z ^ (2 * maxOrder + 1)
    //      !!!  should be less than max_index / 2 !!!
    myMagnet->setMaxOrder(3);
    ofstream fout("3D_quad_");
    Vector_t R(0., 0., 0.), P(3), E(3);
    double t = 0., x, z, s;
    for(x = - 0.02; x <= 0.02 ; x += 0.005) 
      for(z = 0.0; z <= 0.02 ; z += 0.005)
	for(s = -10; s <= 10 ; s += 1) 
	  {
	    R[0] = x;
	    R[1] = z;
	    R[2] = s;
	    Vector_t B(0., 0., 0.);
	    myMagnet->apply(R, P, t, E, B);
	    fout<<x<<' '<<z<<' '<<s<<' '<<B[0]<<' '<<B[1]<<' '<<B[2]<<endl;
	  }
    fout.close();
}

TEST(MultipoleT, Maxwell)
{
    MultipoleT* myMagnet = new MultipoleT("Quadrupole");
    double centralField = 5;
    double fringeLength = 0.5;
    // the highest differential in the fringe field
    double max_index = 20;
    //Set the magnet
    myMagnet->setFringeField(centralField, fringeLength, max_index);
    myMagnet->setTransMaxOrder(1);
    myMagnet->setDipoleConstant(1.0);
    myMagnet->setTransProfile(1, 100.0);
    //highest power in the field is z ^ (2 * maxOrder + 1)
    //      !!!  should be less than max_index / 2 !!!
    myMagnet->setMaxOrder(3);
    //ofstream fout("Quad_CurlB_off");
    Vector_t R(0., 0., 0.), P(3), E(3);
    double t = 0., x, z, s, stepSize= 1e-6;
    for(x = 0.0; x <= 0.0 ; x += 0.001) 
      for(z = 0.0; z <= 0.02 ; z += 0.001)
	for(s = -10; s <= 10 ; s += 0.1) 
	  {
	    R[0] = x;
	    R[1] = z;
	    R[2] = s;
	    Vector_t B(0., 0., 0.);
	    myMagnet->apply(R, P, t, E, B);
	    double div = 0.;
	    div = calcDivB(R, B, stepSize, myMagnet);
	    EXPECT_NEAR(div, 0.0, 0.01);
	    //fout<<x<<' '<<z<<' '<<s<<' '<<B[0]<<' '<<B[1]<<' '<<B[2]<<' '<<div<<endl;
	    vector<double> curl;
	    curl = calcCurlB(R, B, stepSize, myMagnet);
	    EXPECT_NEAR(curl[0], 0.0, 1e-4);
	    EXPECT_NEAR(curl[1], 0.0, 1e-4);
	    EXPECT_NEAR(curl[2], 0.0, 1e-4);
	    //fout<<x<<' '<<z<<' '<<s<<' '<<B[0]<<' '<<B[1]<<' '<<B[2]<<' '<<
	    //sqrt(pow(curl[0], 2)+pow(curl[1],2)+pow(curl[2],2))<<endl;
	  }
    //fout.close();
}

TEST(MultipoleT, Convergence_radius)
{
    MultipoleT* myMagnet = new MultipoleT("Dipole");
    double centralField = 1;
    // the highest differential in the fringe field
    double max_index = 10;
    //Set the magnet
    myMagnet->setTransMaxOrder(0);
    myMagnet->setDipoleConstant(1.0);
    /**for(int order = 1; order <= 10; order ++)
       myMagnet->setTransProfile(order, 1.0);*/
    //highest power in the field is z ^ (2 * maxOrder + 1)
    //      !!!  should be less than max_index / 2 !!!
    myMagnet->setMaxOrder(3);
    ofstream fout("Convergence9");
    Vector_t R(0., 0., 0.), P(3), E(3);
    double t = 0., x, z, s;
    double fringeLength = 0.5;
    myMagnet->setFringeField(centralField, fringeLength, max_index);
    for(x = 0.0; x <= 0.0 ; x += 0.001) 
        for(z = 0.0; z <= 0.7 ; z += 0.01)
	    for(s = -2; s <= 2 ; s += 0.001) 
	        {
		  R[0] = x;
		  R[1] = z;
		  R[2] = s;
		  Vector_t B(0., 0., 0.);
		  myMagnet->apply(R, P, t, E, B);
		  fout<<x<<' '<<z<<' '<<s<<' '<<B[0]<<' '<<B[1]<<' '<<B[2]<<' '<<fringeLength<<endl;
		}
      fout.close();

}

TEST(MultipoleT, Curved_magnet)
{
    MultipoleT* myMagnet = new MultipoleT("Combined function");
    myMagnet->setLength(4.0);
    myMagnet->setBendAngle(1.57);
    myMagnet->setAperture(0.4, 0.4);
    myMagnet->setFringeField(2, 0.5, 0.5); 
    myMagnet->setVarRadius();
    myMagnet->setVarStep(0.1);
    myMagnet->setTransMaxOrder(1);
    myMagnet->setMaxOrder(2);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile(0, 10);
    myMagnet->setTransProfile(1, 1);
    double t=0., x, y, z;
    double stepSize = 1e-7;
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    ofstream fout("Maxwell_test");
    cout<<myMagnet->getBendAngle()<<std::endl;
    for(x = -2.0 ; x <= 2.0 ; x += 0.05) 
        for(y = -6.0; y <= 6.0 ; y += 0.05)
	    for(z = 0.2; z <= 0.21 ; z += 0.1) 
	      {
		R[0] = x;
	        R[1] = z;
		R[2] = y;
		Vector_t B(0., 0., 0.);
		myMagnet->apply(R, P, t, E, B);
		double div = 0.;
	        div = calcDivB(R, B, stepSize, myMagnet);
		vector<double> curl;
	        curl = calcCurlB(R, B, stepSize, myMagnet);
		double abs = 0.0;
		abs += gsl_sf_pow_int(curl[0], 2.0);
		abs += gsl_sf_pow_int(curl[1], 2.0);
		abs += gsl_sf_pow_int(curl[2], 2.0);
		abs = sqrt(abs);
		//if (sqrt(gsl_sf_pow_int(B[0], 2) + gsl_sf_pow_int(B[1], 2) + gsl_sf_pow_int(B[2], 2)) != 0) {
		//  abs /= sqrt(gsl_sf_pow_int(B[0], 2) + gsl_sf_pow_int(B[1], 2) + gsl_sf_pow_int(B[2], 2));
		//}
		if(abs < 10000) {
		  fout<<x<<' '<<y<<' '<<B[0]<<' '<<B[1]<<' '<<B[2]<< ' '<<div<<' '<<abs<<' '<<curl[0]<<' '<<curl[1]<<' '<<curl[2]<<std::endl;
		}
		else {
		  fout<<x<<' '<<y<<' '<<B[0]<<' '<<B[1]<<' '<<B[2]<< ' '<<div<<' '<<0.0<<' '<<curl[0]<<' '<<curl[1]<<' '<<curl[2]<<std::endl;
		}
	      }
		  
    fout.close();
    

}
