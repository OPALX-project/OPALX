function y=mygauss(segma,a,r)
ybound = (1-exp(-a^2/(2*segma^2)))/a;
if (r==0)
    y=0;
elseif ( r> 0 && r <= abs(a) )
    
    y = (1-exp(-r^2/(2*segma^2)))/r;
    
elseif ( r< 0 && r >= -abs(a) )
    
    y = (1-exp(-r^2/(2*segma^2)))/r;
    
else
    y=ybound/r*a;
end

y=y*10;
