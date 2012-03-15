% This function mkes a plots, which shos how difficult it is to integrate
% the function neaded to calculate the wake

a = 5;
s = 140*10^-6 %   s: interpolation point in the wake field
c = 299792458; % speed of ligth
Z0 = 120*pi;    %impedanz
% used material: copper
  sigma = 6.45337 * 10^7; %   sigma: material constant
   tau = 2.70187 *10^-14; %   tau: material constant

k= linspace(0.01,500000,10000);
y = real(Z1ac(a,k,c,Z0,sigma,tau)) .* cos(k*s);


plot(k,y)
title('Function to integrate for the longitudinal wake');
ylabel('Re(Z(k)) cos(ks)');
xlabel('k');


figure

k= linspace(0.01,1);
y = real(c./k.*Z1ac(a,k,c,Z0,sigma,tau)) .* cos(k*s);

plot(k,y)
title('Function to integrate for the transversal wake');
ylabel('Re(Z(k)) cos(ks)');
xlabel('k');