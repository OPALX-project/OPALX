function PrintSignal(Data, fileName, description)
% write Data in a file
%
%Arguments:
% Data:         the data to be written in to the file
% fileName:     name of file
% description:  a description, written after #
%


if nargin<2
    fileName='Data.dat';
end
if nargin<3
   description='';
end

fidPilotSignal=fopen(fileName,'w');
fprintf(fidPilotSignal,description);

%Print out the Data
for i=1:length(Data)
    fprintf(fidPilotSignal,'%d   %d \r\n',i, Data(i));
end


fclose(fidPilotSignal);
