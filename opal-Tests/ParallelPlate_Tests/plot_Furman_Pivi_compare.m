load theory_Furman_Pivi.txt;
load FurmanPivi_const_part.dat;
load FurmanPivi_real_part.dat;
plot(theory_Furman_Pivi(:,1),theory_Furman_Pivi(:,2),'-r', FurmanPivi_const_part(:,1),FurmanPivi_const_part(:,5)/10000,'-.k', FurmanPivi_real_part(:,1),FurmanPivi_real_part(:,4)/2000,'--g');
mylegend=legend('theory', 'OPAL-const-particles', 'OPAL-real-emission');
hold on;
