function makeHeader(data,Name)
% write a array to a a c header file
%
%Arguments:
% data:     the data to be written to the c header file
% Name:     name of the header file
%
data = data(:);

nOutLength = length(data);

%ppath='l:\mnt\share\tables';
fidPilotSignal=fopen(Name,'w');
fprintf(fidPilotSignal,'// Profile of the bunch with %d samples \r\n',nOutLength);
fprintf(fidPilotSignal,'// Generated in Matlab');
fprintf(fidPilotSignal,'#define N_PROFILELENGTH %d\r\n',nOutLength);
fprintf(fidPilotSignal,'const double profile[%d] = { \r\n',nOutLength);

%filename=sprintf(Name);
%filename=sprintf('test.h');


for sample=1:nOutLength-1
    fprintf(fidPilotSignal,'%d, \r\n',data(sample));
end
fprintf(fidPilotSignal,'%d \r\n',data(end));
fprintf(fidPilotSignal,'};\r\n');
fclose(fidPilotSignal);
