function y = W2(s,c,a,Z0,sigma,tau,psi)
% function neaded to calculate the wakefield
%
%Arguments:
%   s: interpolation point in the wake field
%   c: speed of light
%   a: radius in [mm]
%   Z0: impedance
%   sigma: material constant
%   tau: material constant
%   psi: here not used

f = @(k) real(Z2ac(a,k,c,Z0,sigma,tau)) .* cos(k*s);

q=simpson(-5*10^3,5*10^3,10^6,f);

y = 2*10^-19 * c^2/pi * psi* real(q);

end