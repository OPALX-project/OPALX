function compareData(name1,name2,xLabel, yLabel, data1,data2, xStart, xEnd)
% compare 2 diffrent datas and make a plot 
%
%Arguments:
% name1:    Name of the first data 
% name2:    name of the secund data
% xLabel:   xLabel
% yLabel:   yLabel1
% data1:    a arry with the first data
% data2:    a arry with the secund data
% xStart:   start the plot with this x
% xEnd:     end with this x


data1=data1(:);
data2=data2(:);

if nargin<7
   xStart=0;
end
if nargin<8
    xEnd = xStart +length(data1);
end

subplot(2,1,1)
plot(linspace(xStart,xEnd,length(data1)), abs(data1-data2),'LineWidth',2)
title(['Comparison of the ', name1, ' and the ', name2],'FontSize',15)
%legend(['|"',name1,'" - "',name2,'"|'])
xlabel(xLabel,'FontSize',13)
ylabel(['absolut difference in ',yLabel],'FontSize',13)
grid

subplot(2,1,2)
plot(linspace(xStart,xEnd,length(data1)), abs((data1-data2)./data1),'LineWidth',2)
%legend(['|"',name1,'" - "',name2,'" / "',name1,'"|' ])
xlabel(xLabel,'FontSize',13)
ylabel(['relative difference in ',yLabel],'FontSize',13)
grid
end