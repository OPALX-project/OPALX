function y = W1TransAc(s,c,a,Z0,sigma,tau,psi)
% function neaded to calculate the transversal wakefield
%
%Arguments:
%   s: interpolation point in the wake field
%   c: speed of light
%   a: radius in [mm]
%   Z0: impedance
%   sigma: material constant
%   tau: material constant
%   psi: here not used

f = @(k) real(c./k.*Z1ac(a,k,c,Z0,sigma,tau)) .* cos(k*s);

[q, NPoints] = quadl_SP(f,10^-10,10^6,10^-6);%10^-1); 
y = 10^-12 * 2*c/pi * real(q);

end
