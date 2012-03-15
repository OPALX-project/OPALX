function Y=ak(k,l,t,phi,delta_t,vw,vt,lambda)
if (k==1)
    if (l==1)
        Y=mygss_u(t-phi,delta_t,vw,vt,lambda,phi)*seyss_u(t-phi,vw,lambda,phi);
    elseif (l==2)
        Y=mygds_d(t-phi,delta_t,vw,vt,lambda,phi)*seyds_d(t-phi,vw,lambda,phi);
    else
        Y=[];
    end
elseif (k==2)
    if (l==1)
        Y=mygds_u(t-phi,delta_t,vw,vt,lambda,phi)*seyds_u(t-phi,vw,lambda,phi);
    elseif (l==2)
        Y=mygss_d(t-phi,delta_t,vw,vt,lambda,phi)*seyss_d(t-phi,vw,lambda,phi);
    else
        Y=[];
    end
end
y1 = isnan(Y) ; 
if (y1)
    t_phi=t-phi
    K=k
    L=l
    y2=Y
end
        
        