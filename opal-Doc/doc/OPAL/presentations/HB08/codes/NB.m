clear
rms=1
segma_m=rms/0.7773
a_m=segma_m*2.145966;
r=-10:0.05:10;

dist=4;
for (k=1:1:length(r))
    y(k)=mygauss(segma_m,a_m,r(k));
end
for (k=1:1:length(r))
    y1(k)=mygauss(segma_m,a_m,r(k)-dist);
end
for (k=1:1:length(r))
    y2(k)=mygauss(segma_m,a_m,r(k)-2*dist);
end
for (k=1:1:length(r))
    y3(k)=mygauss(segma_m,a_m,r(k)-3*dist);
end
for (k=1:1:length(r))
    y4(k)=mygauss(segma_m,a_m,r(k)-4*dist);
end
for (k=1:1:length(r))
    yp1(k)=mygauss(segma_m,a_m,r(k)+dist);
end
for (k=1:1:length(r))
    yp2(k)=mygauss(segma_m,a_m,r(k)+2*dist);
end
for (k=1:1:length(r))
    yp3(k)=mygauss(segma_m,a_m,r(k)+3*dist);
end
for (k=1:1:length(r))
    yp4(k)=mygauss(segma_m,a_m,r(k)+4*dist);
end
yall_3bunch_dist3 = y+y1+yp1;
yall_5bunch_dist3 = y+y1+y2+yp1+yp2;
yall_7bunch_dist3 = y+y1+y2+y3+yp1+yp2+yp3;
yall_9bunch_dist3 = y+y1+y2+y3+yp1+yp2+yp3+y4+yp4;

plot(r,y,'-',r,yall_3bunch_dist3,':',r,yall_5bunch_dist3,'--',r,yall_7bunch_dist3,'-.','LineWidth',3)

legend('N_B = 1','N_B = 3','N_B = 5', 'N_B = 7',2)
AXIS([-3 3 -6 6])
YLABEL('E_{sc}(r) ( A. U.)','FontSize',24)
XLABEL('r (mm)','FontSize',24)
box off
%%%%%%%%%%%%%%%%%%%%%%%%%
pause 
%%%%%%%%%%%%%%%%%%%%%%%%%
dist=6;
for (k=1:1:length(r))
    y(k)=mygauss(segma_m,a_m,r(k));
end
for (k=1:1:length(r))
    y1(k)=mygauss(segma_m,a_m,r(k)-dist);
end
for (k=1:1:length(r))
    y2(k)=mygauss(segma_m,a_m,r(k)-2*dist);
end
for (k=1:1:length(r))
    y3(k)=mygauss(segma_m,a_m,r(k)-3*dist);
end
for (k=1:1:length(r))
    y4(k)=mygauss(segma_m,a_m,r(k)-4*dist);
end
for (k=1:1:length(r))
    yp1(k)=mygauss(segma_m,a_m,r(k)+dist);
end
for (k=1:1:length(r))
    yp2(k)=mygauss(segma_m,a_m,r(k)+2*dist);
end
for (k=1:1:length(r))
    yp3(k)=mygauss(segma_m,a_m,r(k)+3*dist);
end
for (k=1:1:length(r))
    yp4(k)=mygauss(segma_m,a_m,r(k)+4*dist);
end
yall_3bunch_dist5 = y+y1+yp1;
yall_5bunch_dist5 = y+y1+y2+yp1+yp2;
yall_7bunch_dist5 = y+y1+y2+y3+yp1+yp2+yp3;
yall_9bunch_dist5 = y+y1+y2+y3+yp1+yp2+yp3+y4+yp4;
plot(r,y,'-',r,yall_3bunch_dist5,':',r,yall_5bunch_dist5,'--',r,yall_7bunch_dist5,'-.','LineWidth',3)
legend('N_B = 1','N_B = 3','N_B = 5', 'N_B = 7',2)
AXIS([-3 3 -6 6])
YLABEL('E_{sc}(r) ( A. U.)','FontSize',24)
XLABEL('r (mm)','FontSize',24)
box off
%%%%%%%%%%%%%%%%%%%%%%%%%
pause 
%%%%%%%%%%%%%%%%%%%%%%%%%
dist=8;
for (k=1:1:length(r))
    y(k)=mygauss(segma_m,a_m,r(k));
end
for (k=1:1:length(r))
    y1(k)=mygauss(segma_m,a_m,r(k)-dist);
end
for (k=1:1:length(r))
    y2(k)=mygauss(segma_m,a_m,r(k)-2*dist);
end
for (k=1:1:length(r))
    y3(k)=mygauss(segma_m,a_m,r(k)-3*dist);
end
for (k=1:1:length(r))
    y4(k)=mygauss(segma_m,a_m,r(k)-4*dist);
end
for (k=1:1:length(r))
    yp1(k)=mygauss(segma_m,a_m,r(k)+dist);
end
for (k=1:1:length(r))
    yp2(k)=mygauss(segma_m,a_m,r(k)+2*dist);
end
for (k=1:1:length(r))
    yp3(k)=mygauss(segma_m,a_m,r(k)+3*dist);
end
for (k=1:1:length(r))
    yp4(k)=mygauss(segma_m,a_m,r(k)+4*dist);
end
yall_3bunch_dist7 = y+y1+yp1;
yall_5bunch_dist7 = y+y1+y2+yp1+yp2;
yall_7bunch_dist7 = y+y1+y2+y3+yp1+yp2+yp3;
yall_9bunch_dist7 = y+y1+y2+y3+yp1+yp2+yp3+y4+yp4;
plot(r,y,'-',r,yall_3bunch_dist7,':',r,yall_5bunch_dist7,'--',r,yall_7bunch_dist7,'-.','LineWidth',3)
legend('N_B = 1','N_B = 3','N_B = 5', 'N_B = 7',2)
AXIS([-3 3 -6 6])
YLABEL('E_{sc}(r) ( A. U.)','FontSize',24)
XLABEL('r (mm)','FontSize',24)
box off
pause 
plot(r,y,'-',r,yall_5bunch_dist3,':',r,yall_5bunch_dist5,'--',r,yall_5bunch_dist7,'-.', 'LineWidth',3)
legend('N_B = 1','M = 4,N_B = 5','M = 6, N_B = 5', 'M = 8, N_B = 5',2)
AXIS([-3 3 -6 6])
YLABEL('E_{sc}(r) ( A. U.)','FontSize',24)
XLABEL('r (mm)','FontSize',24)
box off