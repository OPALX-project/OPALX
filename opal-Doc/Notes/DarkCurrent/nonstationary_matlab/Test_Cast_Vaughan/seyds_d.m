function Y=seyds_d(x,vw,lambda,phi_s)
phi_s=mod(phi_s,2*pi);
sey0=0.5;
seym=2.22;
ew0=12.5;%eV
ew1=30.0;%eV
ewm=165.0;%eV
% v_i=myuds(x,lambda,phi_s)*vw;% 
% v_i=(myuds(x,lambda,phi_s)+cos(phi_s)-cos(phi_s+x))*vw;% 
if x~=0
    tmp=(lambda-sin(phi_s)+sin(phi_s+x)-x*cos(phi_s+x))/x*vw;% 
    v_i=tmp/sqrt(1+tmp*tmp/3/1E8/3/1E8);
else
    v_i=0;
end
ewi=(1/sqrt(1-v_i/3.0/1E8*v_i/3.0/1E8)-1)*511000;
if ewi<ew0
    Y=sey0;
    return
else
    u=(ewi-ew0)/(ewm-ew0);
    if u<1
        k=0.56;
    else
        k=0.25;
    end
    if u<=3.6
        Y=seym*power(u*exp(1-u),k);
    else
        Y=seym*1.125/power(u,0.35);
    end
end
        






