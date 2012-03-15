function trajektorie = 2ParticleTest()

%Parameter definition
a = 5;
Csection = 'circular';
Material = 'aluminium';
Conductivity = 'ac';
Lbunch = 294;
Charge = 0.8;
K = 0.20536314319923724;

[x trajektorie] = leapFrog(t0, tend, dt, x0, dx0, xBunch, @ddx);
plot(trajekorie);

end
function a = ddx(x,t)
%Parameter definition
a = 5;
Csection = 'circular';
Material = 'aluminium';
Conductivity = 'ac';
Lbunch = 294;
Charge = 0.8;
K = 0.20536314319923724;

% Generate or compute the wakefield
%wake = cmpwake(a, Csection, Material, Conductivity, Lbunch, Charge, K);
load('wakefield_alu_5_matlab.mat');
% Load the Bunch Profile
load('bunchProfile.mat');

energy = compenergy(K, Charge, wake, bunchProfile);

a = energy(floor(x(1)))*Charge;%/m;

end


