#include"Point.hh"
#include"CarbCycField.hh"
#include<iostream>
#include<cmath>
#include<fstream>
#ifndef DEBUG

using namespace std;
const double pi = 3.14159265358979323846;
int main()
{
  const double Rmin=3000;
  const double Rmax=4600;
  const double dR=10.0;
  const double Thetamin=0.0;
  const double Thetamax=60.0;
  const double dTheta=-5.0;
  const int N = 25;// Polygon with 9 vertexes
  const int N2 = 22;// Polygon with 6 vertexes for fringe field
  const string filename="CarbField.dat";

  const double BMax = 40.0; //Gauss
  const double gap = 60; // mm
  // const double a_para = -0.41744;
  const double a_para = 0.0;
  const double b_para = 1.0;
  const double n_para = 0.8;

  Field *myField = new Field(Rmin,Rmax,dR,Thetamin,Thetamax,dTheta);
  
  Point* V =new Point[N+1];

  V[0].x = 3.30 , V[0].y = -11.5 ;     //A

  V[1].x = 3.40 , V[1].y =  -10.5 ;     //AB1
  V[2].x = 3.50 , V[2].y =  -9.7 ;     //AB2
  V[3].x = 3.60 , V[3].y =  -9.0 ;     //AB3
  V[4].x = 3.70 , V[4].y =  -8.2 ;     //AB4
  V[5].x = 3.80 , V[5].y =  -7.5 ;     //AB5
  V[6].x = 3.90 , V[6].y =  -6.9 ;     //AB6
  V[7].x = 4.00 , V[7].y =  -6.3 ;     //AB7
  V[8].x = 4.10 , V[8].y =  -5.7 ;     //AB8
  V[9].x = 4.20 , V[9].y =  -5.0 ;     //AB9 

  V[10].x = 4.36 , V[10].y =  -4.0 ;     //B
  V[11].x = 4.45 , V[11].y =   2.0 ;     //C 
  V[12].x = 4.50 , V[12].y =   8.0 ;     //D
  V[13].x = 4.45 , V[13].y =  14.0 ;     //E
  V[14].x = 4.36 , V[14].y =  20.0 ;     //F
  
  V[15].x = 4.20 , V[15].y =  18.25 ;     //FG1
  V[16].x = 4.10 , V[16].y =  17.15 ;     //FG2
  V[17].x = 4.00 , V[17].y =  15.9 ;     //FG3
  V[18].x = 3.90 , V[18].y =  14.8 ;     //FG4
  V[19].x = 3.80 , V[19].y =  13.7 ;     //FG5
  V[20].x = 3.70 , V[20].y =  12.5 ;     //FG6
  V[21].x = 3.60 , V[21].y =  11.2 ;     //FG7
  V[22].x = 3.50 , V[22].y =  10.1 ;     //FG8
  V[23].x = 3.40 , V[23].y =   8.9 ;     //FG9

  V[24].x = 3.30 , V[24].y =   8.0 ;     //G
  V[25].x = 3.30 , V[25].y = -11.5 ;     //A
    
  for (int ii=0; ii<N+1; ii++){

    double Angle = V[ii].y + 30.0; // rotate 30 degree.
    double cosA = cos(Angle/180.0*pi);
    double sinA = sin(Angle/180.0*pi);

    double temp_x = cosA*V[ii].x;
    double temp_y = sinA*V[ii].x;

    // (r,theat) --> (x,y) in [mm]
    V[ii].x = temp_x * 1000.0;  
    V[ii].y = temp_y * 1000.0;
      
  }

  myField->setField(V, N, BMax);
  
  /*********************************************/ 
  int Pstart = -20;
  int Pend   = 100;

  const double zStart =-2.0;
  const double zEnd   =10.0;

  const double RA  = 3300;  // G

  const double RAB1 = 3400;
  const double RAB2 = 3500;
  const double RAB3 = 3600;
  const double RAB4 = 3700;
  const double RAB5 = 3800;
  const double RAB6 = 3900;
  const double RAB7 = 4000;
  const double RAB8 = 4100;
  const double RAB9 = 4200;

  const double RB  = 4360;  // F

  const double c0=-0.4785;
  const double c1= 2.9215;
  const double c2=-0.8683;
  const double c3= 0.2000;
  const double c4=-0.0300;
  const double c5= 0.0020;
  
  Point* VR1 =new Point[N2+1];
  
  double thetaA  = -11.5;
  
  double thetaAB1 =  -10.5 ;     //AB1
  double thetaAB2 =  -9.7 ;     //AB2 
  double thetaAB3 =  -9.0 ;     //AB3 
  double thetaAB4 =  -8.2 ;     //AB4 
  double thetaAB5 =  -7.5 ;     //AB5 
  double thetaAB6 =  -6.9 ;     //AB6 
  double thetaAB7 =  -6.3 ;     //AB7 
  double thetaAB8 =  -5.7 ;     //AB8 
  double thetaAB9 =  -5.0 ;     //AB9 

  double thetaB  = -4.5;
  
  double length = abs((zEnd-zStart)*gap);

  double dL = abs( length/(Pend-Pstart) );

  double dz = abs (dL/gap);
    
  double dThetaA  = abs( dL/RA/pi*180.0 );

  double dThetaAB1 = abs( dL/RAB1/pi*180.0 );
  double dThetaAB2 = abs( dL/RAB2/pi*180.0 );
  double dThetaAB3 = abs( dL/RAB3/pi*180.0 );
  double dThetaAB4 = abs( dL/RAB4/pi*180.0 );
  double dThetaAB5 = abs( dL/RAB5/pi*180.0 );
  double dThetaAB6 = abs( dL/RAB6/pi*180.0 );
  double dThetaAB7 = abs( dL/RAB7/pi*180.0 );
  double dThetaAB8 = abs( dL/RAB8/pi*180.0 );
  double dThetaAB9 = abs( dL/RAB9/pi*180.0 );

  double dThetaB  = abs( dL /RB/pi*180.0  );

  
  cout <<"---------start add right fringe field------"<<endl;
  
  for (int index = Pstart; index < Pend; index++ )  
  {

    VR1[0].x = RA,  VR1[0].y = thetaA  - dThetaA *(index) ;     //A

    VR1[1].x = RAB1, VR1[1].y = thetaAB1 - dThetaAB1*(index) ;     //AB
    VR1[2].x = RAB2, VR1[2].y = thetaAB2 - dThetaAB2*(index) ;     //AB
    VR1[3].x = RAB3, VR1[3].y = thetaAB3 - dThetaAB3*(index) ;     //AB
    VR1[4].x = RAB4, VR1[4].y = thetaAB4 - dThetaAB4*(index) ;     //AB
    VR1[5].x = RAB5, VR1[5].y = thetaAB5 - dThetaAB5*(index) ;     //AB
    VR1[6].x = RAB6, VR1[6].y = thetaAB6 - dThetaAB6*(index) ;     //AB
    VR1[7].x = RAB7, VR1[7].y = thetaAB7 - dThetaAB7*(index) ;     //AB
    VR1[8].x = RAB8, VR1[8].y = thetaAB8 - dThetaAB8*(index) ;     //AB
    VR1[9].x = RAB9, VR1[9].y = thetaAB9 - dThetaAB9*(index) ;     //AB
                                
    VR1[10].x = RB,   VR1[10].y = thetaB  - dThetaB *(index) ;     //B
    VR1[11].x = RB,   VR1[11].y = thetaB  - dThetaB *(index+1) ;   //C

    VR1[12].x = RAB9, VR1[12].y = thetaAB9 - dThetaAB9*(index+1) ;     //AB
    VR1[13].x = RAB8, VR1[13].y = thetaAB8 - dThetaAB8*(index+1) ;     //AB
    VR1[14].x = RAB7, VR1[14].y = thetaAB7 - dThetaAB7*(index+1) ;     //AB
    VR1[15].x = RAB6, VR1[15].y = thetaAB6 - dThetaAB6*(index+1) ;     //AB
    VR1[16].x = RAB5, VR1[16].y = thetaAB5 - dThetaAB5*(index+1) ;     //AB
    VR1[17].x = RAB4, VR1[17].y = thetaAB4 - dThetaAB4*(index+1) ;     //AB
    VR1[18].x = RAB3, VR1[18].y = thetaAB3 - dThetaAB3*(index+1) ;     //AB
    VR1[19].x = RAB2, VR1[19].y = thetaAB2 - dThetaAB2*(index+1) ;     //AB
    VR1[20].x = RAB1, VR1[20].y = thetaAB1 - dThetaAB1*(index+1) ;     //AB

    VR1[21].x = RA,  VR1[21].y = thetaA  - dThetaA *(index+1) ;   //D
    VR1[22].x = RA,  VR1[22].y = thetaA  - dThetaA *(index) ;     //A

    // const double temp1 = (index+0.5)*dz - a_para ;
    // const double temp2 = 1 - ( temp1/sqrt( temp1*temp1 + b_para*b_para ) );
    // const double h = 0.5*pow(temp2, n_para );
    const double x = (index+0.5)*dz;
    const double s = c0 + c1*x + c2*pow(x,2.0) + c3*pow(x,3.0) + c4*pow(x,4.0) + c5*pow(x,5.0);
    const double h = 1.0/(1+exp(s));
    
    cout<<"Start Fringe field..."<<endl;    
    const double temp3 = VR1[0].y;
    cout<<"x = "<<x<<" [gap], h = " <<h<<endl;
    cout<<"Area index = "<<index<<endl;
    cout<<"Minial angle = "<<temp3<<endl;

    for (int ii=0; ii<N2+1; ii++){

      double Angle = VR1[ii].y + 30.0; // rotate 30 degree.
      double cosA = cos(Angle/180.0*pi);
      double sinA = sin(Angle/180.0*pi);

      double temp_x = cosA*VR1[ii].x;
      double temp_y = sinA*VR1[ii].x;

      // (r,theat) --> (x,y) in [mm]
      VR1[ii].x = temp_x ;  
      VR1[ii].y = temp_y ;
      cout<<VR1[ii].x<<" "<< VR1[ii].y<<endl;  
    }
    cout<<"Bfield [kGs] = " <<h*BMax<<endl;
    cout<<"End Fringe field "<<endl;
    
    myField->setField(VR1, N2, h*BMax);
  }

  /*********************************************/ 
  delete[] VR1;
  VR1 = new Point[N2+1];

  //  Pstart = -20;
  //  Pend   = 100;

  thetaB  = 20.0;

  thetaAB9 =  18.25  ;    //AB1
  thetaAB8 =  17.15 ;     //AB2 
  thetaAB7 =  15.9  ;     //AB3 
  thetaAB6 =  14.8  ;     //AB4 
  thetaAB5 =  13.7  ;     //AB5 
  thetaAB4 =  12.5  ;     //AB6 
  thetaAB3 =  11.2  ;     //AB7 
  thetaAB2 =  10.1  ;     //AB8 
  thetaAB1 =   8.9  ;     //AB9 

  thetaA  =  8.0;


  length = abs((zEnd-zStart)*gap);

  dL = abs( length/(Pend-Pstart) );

  dz = abs (dL/gap);

  dThetaA  = abs( dL/RA/pi*180.0 );

  dThetaAB1 = abs( dL/RAB1/pi*180.0 );
  dThetaAB2 = abs( dL/RAB2/pi*180.0 );
  dThetaAB3 = abs( dL/RAB3/pi*180.0 );
  dThetaAB4 = abs( dL/RAB4/pi*180.0 );
  dThetaAB5 = abs( dL/RAB5/pi*180.0 );
  dThetaAB6 = abs( dL/RAB6/pi*180.0 );
  dThetaAB7 = abs( dL/RAB7/pi*180.0 );
  dThetaAB8 = abs( dL/RAB8/pi*180.0 );
  dThetaAB9 = abs( dL/RAB9/pi*180.0 );

  dThetaB  = abs( dL /RB/pi*180.0  );

  
  cout <<"---------start add left fringe field------"<<endl;
  
  for (int index = Pstart; index < Pend; index++ )  
  {

    VR1[0].x = RA,  VR1[0].y = thetaA  + dThetaA *(index) ;     //A

    VR1[1].x = RAB1, VR1[1].y = thetaAB1 + dThetaAB1*(index) ;     //AB
    VR1[2].x = RAB2, VR1[2].y = thetaAB2 + dThetaAB2*(index) ;     //AB
    VR1[3].x = RAB3, VR1[3].y = thetaAB3 + dThetaAB3*(index) ;     //AB
    VR1[4].x = RAB4, VR1[4].y = thetaAB4 + dThetaAB4*(index) ;     //AB
    VR1[5].x = RAB5, VR1[5].y = thetaAB5 + dThetaAB5*(index) ;     //AB
    VR1[6].x = RAB6, VR1[6].y = thetaAB6 + dThetaAB6*(index) ;     //AB
    VR1[7].x = RAB7, VR1[7].y = thetaAB7 + dThetaAB7*(index) ;     //AB
    VR1[8].x = RAB8, VR1[8].y = thetaAB8 + dThetaAB8*(index) ;     //AB
    VR1[9].x = RAB9, VR1[9].y = thetaAB9 + dThetaAB9*(index) ;     //AB
                                
    VR1[10].x = RB,   VR1[10].y = thetaB  + dThetaB *(index) ;     //B
    VR1[11].x = RB,   VR1[11].y = thetaB  + dThetaB *(index+1) ;   //C

    VR1[12].x = RAB9, VR1[12].y = thetaAB9 + dThetaAB9*(index+1) ;     //AB
    VR1[13].x = RAB8, VR1[13].y = thetaAB8 + dThetaAB8*(index+1) ;     //AB
    VR1[14].x = RAB7, VR1[14].y = thetaAB7 + dThetaAB7*(index+1) ;     //AB
    VR1[15].x = RAB6, VR1[15].y = thetaAB6 + dThetaAB6*(index+1) ;     //AB
    VR1[16].x = RAB5, VR1[16].y = thetaAB5 + dThetaAB5*(index+1) ;     //AB
    VR1[17].x = RAB4, VR1[17].y = thetaAB4 + dThetaAB4*(index+1) ;     //AB
    VR1[18].x = RAB3, VR1[18].y = thetaAB3 + dThetaAB3*(index+1) ;     //AB
    VR1[19].x = RAB2, VR1[19].y = thetaAB2 + dThetaAB2*(index+1) ;     //AB
    VR1[20].x = RAB1, VR1[20].y = thetaAB1 + dThetaAB1*(index+1) ;     //AB

    VR1[21].x = RA,  VR1[21].y = thetaA  + dThetaA *(index+1) ;   //D
    VR1[22].x = RA,  VR1[22].y = thetaA  + dThetaA *(index) ;     //A

    // const double temp1 = (index+0.5)*dz - a_para ;
    // const double temp2 = 1 - ( temp1/sqrt( temp1*temp1 + b_para*b_para ) );
    // const double h = 0.5*pow(temp2, n_para );
    const double x = (index+0.5)*dz;
    const double s = c0 + c1*x + c2*pow(x,2.0) + c3*pow(x,3.0) + c4*pow(x,4.0) + c5*pow(x,5.0);
    const double h = 1.0/(1+exp(s));
    
    cout<<"Start Fringe field..."<<endl;    
    const double temp3 = VR1[0].y;
    cout<<"x = "<<x<<" [gap], h = " <<h<<endl;
    cout<<"Area index = "<<index<<endl;
    cout<<"Minial angle = "<<temp3<<endl;

    for (int ii=0; ii<N2+1; ii++){

      double Angle = VR1[ii].y + 30.0; // rotate 30 degree.
      double cosA = cos(Angle/180.0*pi);
      double sinA = sin(Angle/180.0*pi);

      double temp_x = cosA*VR1[ii].x;
      double temp_y = sinA*VR1[ii].x;

      // (r,theat) --> (x,y) in [mm]
      VR1[ii].x = temp_x ;  
      VR1[ii].y = temp_y ;
      cout<<VR1[ii].x<<" "<< VR1[ii].y<<endl;  
    }
    cout<<"Bfield [kGs] = " <<h*BMax<<endl;
    cout<<"End Fringe field "<<endl;
    
    myField->setField(VR1, N2, h*BMax);
  }

  /*

  thetaA  = -11.5;
  thetaAB = -8.5;
  thetaB  = -4.5;

  length = abs((zEnd-zStart)*gap);

  dL = abs( length/(Pend-Pstart) );

  dz = abs (dL/gap);
    
  dThetaA  = abs( dL/RA/pi*180.0 );
  dThetaAB = abs( dL/RAB/pi*180.0 );
  dThetaB  = abs( dL /RB/pi*180.0  );

    
  for (int index = Pstart; index < Pend; index++ )  
  {
    VR1[0].x = RA,  VR1[0].y = thetaA  - dThetaA *(index) ;     //A
    VR1[1].x = RAB, VR1[1].y = thetaAB - dThetaAB*(index) ;     //AB
    VR1[2].x = RB,  VR1[2].y = thetaB  - dThetaB *(index) ;     //B
    VR1[3].x = RB,  VR1[3].y = thetaB  - dThetaB *(index+1) ;   //C
    VR1[4].x = RAB, VR1[4].y = thetaAB - dThetaAB*(index+1) ;   //CD
    VR1[5].x = RA,  VR1[5].y = thetaA  - dThetaA *(index+1) ;   //D
    VR1[6].x = RA,  VR1[6].y = thetaA  - dThetaA *(index) ;     //A

    // const double temp1 = (index+0.5)*dz - a_para ;
    //const double temp2 = 1 - ( temp1/sqrt( temp1*temp1 + b_para*b_para ) );
    //const double h = 0.5*pow(temp2, n_para );

    const double x = (index+0.5)*dz;
    const double s = c0 + c1*x + c2*pow(x,2.0) + c3*pow(x,3.0) + c4*pow(x,4.0) + c5*pow(x,5.0);
    const double h = 1.0/(1+exp(s));

    cout<<"Start Fringe field..."<<endl;    
    const double temp3 = VR1[0].y;
    cout<<"x = "<<x<<" [gap], h = " <<h<<endl;
    cout<<"Area index = "<<index<<endl;
    cout<<"Minial angle = "<<temp3<<endl;

    for (int ii=0; ii<N2+1; ii++){

      double Angle = VR1[ii].y + 30.0; // rotate 30 degree.
      double cosA = cos(Angle/180.0*pi);
      double sinA = sin(Angle/180.0*pi);

      double temp_x = cosA*VR1[ii].x;
      double temp_y = sinA*VR1[ii].x;

      // (r,theat) --> (x,y) in [mm]
      VR1[ii].x = temp_x ;  
      VR1[ii].y = temp_y ;
      cout<<VR1[ii].x<<" "<< VR1[ii].y<<endl;  
    }
    cout<<"Bfield [kGs] = " <<h*BMax<<endl;
    cout<<"End Fringe field "<<endl;
    
    myField->setField(VR1, N2, h*BMax);
  }
 */

  myField->saveField(filename);

  // cout<<"Field Generation Done!"<<endl;

  ofstream outf;
  outf.setf(ios::scientific, ios::floatfield );
  outf.precision(8);
  outf.open("geom.data");

  Point* V1 =new Point[N+1];
  
  for (int jj = 0 ;jj< 6; jj++) {
    
    double rotateAngle = -60.0*jj;
    double cosA = cos(rotateAngle/180.0*pi);
    double sinA = sin(rotateAngle/180.0*pi);
  
    for (int ii=0; ii<N+1; ii++){
      V1[ii].x = cosA*V[ii].x + sinA*V[ii].y;
      V1[ii].y = -sinA*V[ii].x + cosA*V[ii].y;
    }

    for (int ii=0; ii<N+1; ii++){
      outf<< V1[ii].x<<" "<< V1[ii].y<<endl;
    }
  }

  
  outf.close();
  
  
  return 1;
}

#endif


#ifdef DEBUG
#include"CheckPoint.hh"
// test if a point is in a polygon

int main() 
{
  using namespace checkPoint;
  Point P;
  Point *V;
  
  
  P.x=2.9999999;
  P.y=0.0;

  V = new Point[3];

  V[0].x= 0.0;
  V[0].y= 0.0;

  V[1].x= 3.0;
  V[1].y= 0.0;

  V[2].x= 0.0;
  V[2].y= 4.0;
 
  int result_rn = wn_PnPoly(P, V , 3);
  std::cout<<"Result by Winding is:";
  if (result_rn == 0)
    std::cout<<"Out"<<std::endl;
  else
    std::cout<<"In"<<std::endl;
  
  int result_cn = cn_PnPoly(P, V , 3);
                            
  std::cout<<"Result by Crossing is:";
  if (result_cn == 0)
    std::cout<<"Out"<<std::endl;
  else
    std::cout<<"In"<<std::endl;

  return 1;
}
#endif
