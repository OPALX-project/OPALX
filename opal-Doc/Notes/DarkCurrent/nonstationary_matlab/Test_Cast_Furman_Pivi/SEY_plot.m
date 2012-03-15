fid = fopen('SEY_result.txt', 'w', 'n');
for i=2:n
    tmp_z=[(i-1)*h/w*1E9,(myf(1,i)+myf(2,i))/(myI(1,i)+myI(2,i))];
    fprintf(fid,'%12.8f %12.8f\n', tmp_z);
end
fclose(fid);
load SEY_result.txt;
plot(SEY_result(:,1),SEY_result(:,2));
hold on;   
% load impactNum.dat;
% fid = fopen('impact.dat', 'w');
% j=1;
% while (j<=size(impactNum,1))
%     if impactNum(j,3)==0.0
%         impactNum(j,:)=[];
%     else
%         impactNum(j,1)=impactNum(j,2)/1e-9;
%         impactNum(j,2)=impactNum(j,4)/impactNum(j,3);
%         
%         fprintf(fid,'%8.3f %8.3f\n',impactNum(j,1),impactNum(j,2));
%         j=j+1;
%     end
% end
% fclose(fid);
% load impact.dat;
% n_t=size(impact,1);
% tmp_s=zeros(n_t,1);
% for j=1:n_t
%     tmp_s(j,1)=impact(j,1);
% end
% plot(tmp_s(:,1),impact(:,2));