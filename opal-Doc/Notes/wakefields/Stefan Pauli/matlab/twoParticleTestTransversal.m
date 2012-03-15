function [x deltaX] = twoParticleTestTransversal(t0, tend, dt,xStart,a, Direction, Csection, Material, Conductivity, Lbunch, Charge, K, c, wake)
% a fuction to calculate the force resulting form the two particle model
% with a transversal wake
% 
%Arguments:
% t0:           integrate from t = to
% tend:         to t = tend
% dt:           int steps of dt
% xStart:       start at x = xStart
% a:            radius of the beam pipe in [mm]
% Direction:    'longitudinal' or 'transversal'
% Csection:     'circular' (or 'rectangular', not jet implemented)
% Material:     'copper' or 'aluminium'
% Conductivity: 'ac' or 'dc'
% Lbunch:       Length of the computed Wakefield [samples]
% Charge:       Charge of one particle [nC]
% K:            normalize the length distribution of the particle bunch
% c:            speed of light
% wake:         the wake function


%Parameter definition
if nargin<5
    a = 5;
end
if nargin<6
    Direction = 'longitudinal';
end
if nargin<7
    Csection = 'circular';
end
if nargin<8
    Material = 'aluminium';
end
if nargin<9
    Conductivity = 'ac';
end
if nargin<10
    Lbunch = 294;
end
if nargin<11
    Charge = 0.8;
end
if nargin<12
    K = 0.20536314319923724;
end
if nargin<13
    c = 299792458;
end
if nargin<14
    wake = cmpwake(a, Direction, Csection, Material, Conductivity, Lbunch, Charge, K);
end
    % Load the Bunch Profile
    load('bunchProfile.mat');
    energy = compenergy(K, Charge, wake, bunchProfile);
    energy = [energy(:); ones(10,1) * energy(end)];

%Parameter definition
a = 5;
Csection = 'circular';
Material = 'aluminium';
Conductivity = 'ac';
Lbunch = 294;
Charge = 0.8;
K = 0.20536314319923724;
c = 299792458;

[x deltaX] = leapFrogTransversal(t0, tend, dt, xStart, [0; 0; c*.9], @ddx, energy);
%plot(trajektorie);
end
function a = ddx(x,t,energy)
%Parameter definition
a = 5;
Csection = 'circular';
Material = 'aluminium';
Conductivity = 'ac';
Lbunch = 294;
Charge = 0.8;
K = 0.20536314319923724;

%energy = compenergy(K, Charge, wake, bunchProfile);

a = energy(floor(x(3)))*Charge;%/m;
a = [a ; 0; 0];

end

