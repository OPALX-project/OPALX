function y = compenergy(K, Charge, wake, profile)
% computes the energy sperad
%
%Arguments:
% K:        normalize the length distribution of the particle bunch
% Charge:   Charge of one particle [nC]
% wake:     wakefield
% profile:  length density of the particle bunch

y = 10^-6 * (Charge * 10^3) * K * filter(wake,1,profile);
y = - y / 1000;

