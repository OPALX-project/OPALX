function[Qs]=simpson(a,b,n,func)
% A implementation of the simpson algorithm, used to integrate a function
%
%Arguments:
% a:    integrate from a
% b:    integrate to b
% n:    Number of integration points  
% f:    the function to integrat

h=(b-a)/(2*n);
y=func(a:h:b);
Qs=(y(1)+y(end)+4*sum(y(2:2:end-1))+2*sum(y(3:2:end-2)))*h/3;
end