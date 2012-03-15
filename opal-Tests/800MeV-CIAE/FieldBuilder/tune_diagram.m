clc
load('CYC4.mat');
load('CYC_OPAL.mat');
CYC = CYC4;
plot(CYC(:,1),CYC(:,2),'r');grid on;
xlabel('Energy(MeV)');ylabel('(w0/w - 1)(x1000)');
pause
a = 1;b = length(CYC(:,1));
plot(CYC(a:b,1),CYC(a:b,3),'-+r',CYC(a:b,1),CYC(a:b,4),'-+b',CYC(a:b,1),CYC_OPAL(a:b,1),'-vb',CYC(a:b,1),CYC_OPAL(a:b,2),'-vr');
xlabel('Energy(MeV)');ylabel('Nur or Nuz');
legend('Vr -- CYCLOP','Vz -- CYCLOP','Vr -- OPAL', 'Vz -- OPAL');
pause

plot(CYC(a:b,3),CYC(a:b,4),'-+r',CYC_OPAL(a:b,1),CYC_OPAL(a:b,2),'-xb','LineWidth',2);
xlabel('Vr');ylabel('Vz');
legend('CYCLOP','OPAL');
axis equal;
axis([0.8,2.0,0.6,2.0]);
hold on;

x=[2,2];y=[0.6,2.0];plot(x,y,'k','LineWidth',4);
x=[1,1];y=[0.6,2.0];plot(x,y,'k','LineWidth',4);
x=[0.8,2];y=[1,1];plot(x,y,'k','LineWidth',4);
x=[0.8,2];y=[2,2];plot(x,y,'k','LineWidth',4);


x=[0.8,2]; y = [1.5,1.5];plot(x,y,'k','LineWidth',2);
x=[1.5,1.5];y=[0.6,2.0];plot(x,y,'k','LineWidth',2);

x=[0.8,2]; y = 2-x;plot(x,y,'k','LineWidth',2);
x=[0.8,2]; y = 3-x;plot(x,y,'k','LineWidth',2);

x=[0.8,2]; y = x-1;plot(x,y,'k','LineWidth',2);
x=[0.8,2]; y = x;plot(x,y,'k','LineWidth',2);
x=[0.8,2]; y = x+1;plot(x,y,'k','LineWidth',2);

x=[0.8,2]; y = [2/3,2/3];plot(x,y,'--k','LineWidth',1);
x=[0.8,2]; y = 3-2*x;plot(x,y,'--k','LineWidth',1);

for i = 4:5
    x=[0.8,2];y = [i/3,i/3];plot(x,y,'--k','LineWidth',1);
    x=[i/3,i/3];y = [0.6,2.0];plot(x,y,'--k','LineWidth',1);
    x=[0.8,2]; y = i-2*x;plot(x,y,'--k','LineWidth',1);
    x=[0.8,2]; y = (i-x)/2;plot(x,y,'--k','LineWidth',1);
end
for i = -3:0
    x=[0.8,2];y = i+2*x;plot(x,y,'--k','LineWidth',1);
    x=[0.8,2];y = (x-i)/2;plot(x,y,'--k','LineWidth',1);
end

