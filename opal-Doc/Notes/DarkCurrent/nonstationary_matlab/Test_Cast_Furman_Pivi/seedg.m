function Y=seedg(k,x)
if k==2 %down plate;
   
    syms t;
    Y=subs(dirac(t),t,x);
% Y=1;
else
    Y=0;
end