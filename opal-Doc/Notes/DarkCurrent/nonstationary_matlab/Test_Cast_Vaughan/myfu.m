function Y=myfu(x,vw,vt)
vw_sq=vw*vw;
vt_sq=vt*vt;
x_sq=x*x;
Y=x*vw_sq/vt_sq*exp(-x_sq*vw_sq/2/vt_sq);
% Y=sqrt(2/pi)*vw/vt*exp(-x_sq*vw_sq/2/vt_sq);