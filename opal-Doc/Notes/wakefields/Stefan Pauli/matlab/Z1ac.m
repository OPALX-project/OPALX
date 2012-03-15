function y = Z1ac(a,k,c,Z0,sigma,tau)
% function neaded to calculate the wakefield
%
%Arguments:
%   a: radius in [mm]
%   k: wave number
%   c: speed of light
%   Z0: impedance
%   sigma: material constant
%   tau: material constant
y=sqrt(Z0*abs(k)/2) ...
 .* sqrt( sigma./(1-i*c*k*tau)) ...
 .* (i+sign(k));
y = y.*k.^-1;
y = 1./ (y-(i*k*(10^-3*a)/2)); 
y = Z0./(2*pi*10^-3*a)*y;

end