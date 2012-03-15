load theory_Vaughan.txt;
load Vaughan_const_part.dat;
load Vaughan_real_part.dat;
plot(theory_Vaughan(:,1),theory_Vaughan(:,2),'-r', Vaughan_const_part(:,1),Vaughan_const_part(:,5)/20000,'-.k', Vaughan_real_part(:,1),Vaughan_real_part(:,4)/2000,'--g');
mylegend = legend('theory', 'OPAL-const-particles', 'OPAL-real-emission');
hold on;
