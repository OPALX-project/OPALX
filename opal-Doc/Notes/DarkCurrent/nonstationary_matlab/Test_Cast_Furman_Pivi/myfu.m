function Y=myfu(x,vw,vt)
% tmp=1.3271121*1E6/vw;
% syms t;
% Y=subs(dirac(t-tmp),t,x);
% if abs(x-tmp)<=0.001
%     Y=1;
% else
%     Y=0;
% end
vw_sq=vw*vw;
vt_sq=vt*vt;
x_sq=x*x;
Y=x*vw_sq/vt_sq*exp(-x_sq*vw_sq/2/vt_sq);