double calculateAngle(const double x,const double y)
{

  double thetaXY;
  if      ((x>0) && (y>=0)) thetaXY=atan(y/x); 
  else if ((x<0) && (y>=0)) thetaXY=pi+atan(y/x);
  else if ((x<0) && (y<=0)) thetaXY=pi+atan(y/x); 
  else if ((x>0) && (y<=0)) thetaXY=2.0*pi+atan(y/x);
  else if ((x==0) && (y> 0)) thetaXY=pi/2.0;
  else if ((x==0) && (y< 0)) thetaXY=3.0/2.0*pi;
  cout<<"angle = "<<thetaXY/pi*180.0<<endl;
  
  return thetaXY;

}
