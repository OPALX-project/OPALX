fid = fopen('SEY_result.txt', 'w', 'n');
for i=2:n
    tmp_z=[(i-1)*h/w*1E9,(myf(1,i)+myf(2,i))/(myI(1,i)+myI(2,i))];
    fprintf(fid,'%12.8f %12.8f\n', tmp_z);
end
fclose(fid);
load SEY_result.txt;
plot(SEY_result(:,1),SEY_result(:,2));
hold on;   
