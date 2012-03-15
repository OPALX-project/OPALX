% fid = fopen('fig8_result.txt', 'w', 'n');
% for i=2:n
%     tmp_z=[(i-1)*h/w*1E9,myNp(i,1)];
% 
%     fprintf(fid,'%12.8f %12.8f\n', tmp_z);
% end
% fclose(fid);
load fig8_result.txt;
% plot(fig8_result(:,1),fig8_result(:,2));
% hold on
load fig_simu120v.dat;
n_t=size(fig_simu120v,1);
tmp_ss=zeros(n_t,1);
for j=1:n_t
    tmp_ss(j,1)=(j-1)*1e-11/1e-9;
end
% plot(tmp_s(:,1),fig_simu120v(:,2)/5000);
    
load fig_simu120v1000.dat;
n_t=size(fig_simu120v1000,1);
tmp_s=zeros(n_t,1);
for j=1:n_t
    tmp_s(j,1)=(j-1)*1e-11/1e-9;
end
load const_part_benchmark_Furman.dat;
%plot(fig8_result(:,1),fig8_result(:,2),'-r',tmp_ss(:,1),fig_simu120v(:,2)/5000,'-.g',tmp_s(:,1),fig_simu120v1000(:,2)/1000,'--b',constant_PartNum_10K_Np(:,1),constant_PartNum_10K_Np(:,5)/100000,'--r');
%lg1=legend('theory','OPAL-5000','OPAL-1000','10k_const');
plot(fig8_result(:,1),fig8_result(:,2),'-r',const_part_benchmark_Furman(:,1),const_part_benchmark_Furman(:,5)/50000,'-.k',tmp_s(:,1),fig_simu120v1000(:,2)/1000,'--g');
lg1=legend('theory','OPAL-const-particles','OPAL-real-emission');