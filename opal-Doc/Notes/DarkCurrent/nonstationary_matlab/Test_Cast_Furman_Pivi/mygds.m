function Y=mygds(x,delta_t,vw,vt,lambda,phi_s)
deri=abs((myuds(x+delta_t,lambda,phi_s)-myuds(x-delta_t,lambda,phi_s))/2/delta_t);
Y=deri*myfu(myuds(x,lambda,phi_s),vw,vt);

