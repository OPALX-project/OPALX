function y = enge_empirical_fun(x)
    c =[-0.04785, 2.9215, -0.8683, 0.2000, -0.0300,0.0020];
    if( x < -2.0)
        y = 1.0;
    elseif( x > 9)
        y = 0.0;
    else  
        y = c(1)*ones(1,length(x));
        for i = 1 : 5
            y = y + c(i+1)*x.^i;
        end
        y = 1.0./(1+exp(y));
    end
end