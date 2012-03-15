function test()
% just used for tests

c = 299792458;
Z0 = 120*pi;
a = 0.75;
   sigma = 6.45337 * 10^7;
   tau = 2.70187 *10^-14;

   s=10^-4;
 w = @W1ac;
 
  w(s,c,a,Z0,sigma,tau)
 
%W1ac(s,c,a,Z0,sigma,tau)

% c = 299792458;
% Z0 = 120*pi;
% a = 0.75;
%    sigma = 6.45337 * 10^7;
%    tau = 2.70187 *10^-14;
%    
%    k=10^5
% real(Z1ac(a,k,c,Z0,sigma,tau))

% f = @(k) real(Z1ac(a,k,c,Z0,sigma,tau)) .* cos(k*s);
% 
% plot(f([0:10:10^6]))

end

