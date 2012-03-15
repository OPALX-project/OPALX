function Y=seedg(k,x)
if k==1 || k==2
    syms t;
    Y=subs(0.5*dirac(t),t,x);
% Y=1;
else
    Y=[];
end