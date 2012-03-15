function profile = generateProfile()
% generate a profile

x=linspace(-10,10,160);
profile = 1/sqrt(2*pi) * exp(-1/2*x.*x);
