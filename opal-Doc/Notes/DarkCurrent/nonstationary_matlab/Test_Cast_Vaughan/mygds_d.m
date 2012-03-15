function Y=mygds_d(x,delta_t,vw,vt,lambda,phi_s)
phi_s=mod(phi_s,2*pi);
if (x~=0)
    deri=abs((myuds(x+delta_t,lambda,phi_s)-myuds(x-delta_t,lambda,phi_s))/2/delta_t);
    Y=deri*myfu(myuds(x,lambda,phi_s),vw,vt);
else
    Y=0;
end
% assert(Y>=0,'x=%f ,phi_s=%f,',x,phi_s)
