function makePlots()
% make several plots with diffrent parameters

%% definitions
a = 5;
Csection = 'circular';
Material = 'aluminium';
Direction = 'longitudinal';
Conductivity = 'ac';
Lbunch = 294;
Charge = 0.8;
K = 0.20536314319923724;

%% Compare the wakefield form matlab with the wakefield from mathematica
load ('wakefield_alu_5Mathematica.mat')

wake = cmpwake(a, Direction, Csection, Material, Conductivity, Lbunch, Charge, K);
%PrintSignal(wake, Wake_L_a_c_ac_Mat.dat, '//The Wake Filed calculatet in Matlab; Circular, Alu, longitudinal, ac, a=5')
compareData('wakefield from matlab','wakefield from mathematica','s[um]', '[V/pC/m]', wake,wakefield_alu_5Mathematica)

%% Compare the energy spread from matlab with the energy spread from
%mathematica
load ('energy_alu_5_mathematica.mat')
load('bunchProfile.mat')

energy = compenergy(K, Charge, wake, bunchProfile);
figure
compareData('energy sprea from matlab','energy sprea from mathematica','s[um]', '\DeltaE(s) [keV/m]', energy,energy_alu_5_mathematica)

%% Compute the longitudinal wakefield
Direction = 'longitudinal';
figure
%wake = cmpwake(a, Direction, Csection, Material, Conductivity, Lbunch, Charge, K);
plot(wake)
title('Longitudinal wake function')
xlabel('s[\mum]')
ylabel('W_T(s) [V/pC/m^2] Stimmt Das?????????????')

%% Compute the transversal wakefield
Direction = 'transversal';
figure
wakeT = cmpwake(a, Direction, Csection, Material, Conductivity, Lbunch, Charge, K);
plot(wakeT)
title('Transversal wake function')
xlabel('s[\mum]')
ylabel('W_T(s) [V/pC/m^2] Stimmt Das?????????????')
Direction = 'longitudinal';

%% Two particle model transversal
% Doas not work 
% 
% energy = compenergy(K, Charge, wakeT, bunchProfile);
% energy = [energy(:); ones(10,1) * energy(end)]; 
% 
% t =1
% for i = 2:Lbunch
% [x,dx] = twoParticleTestTransversal(0, t, 0.0001,[1 0 i], energy);
% deltaX(i-1) = dx;
% end
% figure
% plot (deltaX)
% title(['Comparison of startposition and position after ' num2str(t) 's'])
% xlabel('(head)                                                               s[\mum]                                                               (tail)')
% ylabel('diffrenc in the relative position [\mum]')
%% Two particle model longitudinal
t=5

for i = 2:Lbunch
[x,dx] = twoParticleTest(0, t, 0.0001,i);
deltaX(i-1) = dx;
end
figure
plot (deltaX)
title(['Comparison of startposition and position after ' num2str(t) 's'])
xlabel('(head)                                                               s[\mum]                                                               (tail)')
ylabel('diffrenc in the relative position [\mum]')


% figure
% plot (deltaX/max(-deltaX))
% xlabel('(head)               s/\mum               (tail)')
% hold on
% plot (energy/max(-energy),'r')
% legend('deltaX','energy')
