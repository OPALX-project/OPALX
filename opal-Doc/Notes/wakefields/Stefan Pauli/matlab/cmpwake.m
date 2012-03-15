function wake = cmpwake(a, Direction, Csection,  Material, Conductivity, Lbunch, Charge, K)
%computes the wakefield with the given parameters
%
%Arguments:
%   a:              radius of the beam pipe in [mm]
%   Direction:      'longitudinal' or 'transversal'
%   Csection:       'circular' (or 'rectangular', not jet implemented)
%   Material:       'copper' or 'aluminium'
%   Conductivity:   'ac' or 'dc'
%   Lbunch:         Length of the computed Wakefield [samples]
%   Charge:         Charge of one particle [nC]
%   K:              normalize the length distribution of the particle bunch
if nargin<1
    a = 5%0.75;
end
if nargin<2
    Direction = 'longitudinal';
end
if nargin<3
    Csection = 'circular';
end
if nargin<4
    Material = 'aluminium'%'copper';
end
if nargin<5
    Conductivity = 'ac';
end
if nargin<6
    Lbunch = 294;
end
if nargin<7
    Charge = 0.8;
end
if nargin<8
    K = 0.20536314319923724;
end

% Define constants
c = 299792458;
Z0 = 120*pi;
psi = 1000;

if strcmp(Material, 'copper')%Material == 'copper111'
   sigma = 6.45337 * 10^7;
   tau = 2.70187 *10^-14;
   r1 = 0;
   r3 = 1;
end


if strcmp(Material, 'aluminium')%Material == 'aluminium'
   sigma = 4.22807 * 10^7;
   tau = 8.00554 *10^-15;
   r1 = 1;
   r3 = 0;
end
if strcmp(Direction, 'longitudinal')
    if Csection == 'circular'
        if Conductivity == 'ac';
            W = @W1ac;
        end
        if Conductivity == 'dc';
            W = @W1;
        end
    end
end
if strcmp(Direction, 'transversal')
    if Csection == 'circular'
        if Conductivity == 'ac';
            W = @W1TransAc;
        end
        if Conductivity == 'dc';
            W = @W1Trans;
        end
    end
end

%start with the code
tic
for j = 0:Lbunch-1
    if(mod(j,10)==0)
        j
    end
    wake(j+1) = real(W(j*10^-6,c,a,Z0,sigma,tau,psi));
end
wake=wake(:);
disp(['Used time in cmpwake.m: ', num2str(toc)]); 

% for j=1:Lbunch-1
%     s = j*10^-6;
%     wake(j) = - 1/(2*pi*a)*sqrt(c/sigma) * 1/(s^(3/2)) * a^2/4;
%     %wake(j) = 4/a^2*exp(-(s/(4*c*tau))) * cos(sqrt(2*sqrt(4*pi*sigma/tau)/a)*s);
% end


