function Y=mygss(x,delta_t,vw,vt,lambda,phi_s)
deri=abs((myuss(x+delta_t,lambda,phi_s)-myuss(x-delta_t,lambda,phi_s))/2/delta_t);
Y=deri*myfu(myuss(x,lambda,phi_s),vw,vt);
