function Y=seyds_u(x,vw,lambda,phi_s)
phi_s=mod(phi_s+pi,2*pi);

sePScat_m = 0.02;
sePScatPeak_m = 0.496;
seEScatPeak_m = 0;
seW_m = 60.86;
seP_m = 1.0;

seEOne_m = 0.26;
seETwo_m = 2;

sePRed_m = 0.2;
seERed_m = 0.041;
seR_m = 0.104;

seROne_m = 0.26;
seRTwo_m = 2;

seYPeakTS_m = 1.8848;
seEPeakTS_m = 276.8;
seSTS_m = 1.54;
seTOneTS_m = 0.66;
seTTwoTS_m = 0.8;
seTThreeTS_m = 0.7;
seTFourTS_m = 1.0;
cosTheta=1;
%v_i=myuds(x,lambda,phi_s)*vw;% 
%v_i=(myuds(x,lambda,phi_s)+cos(phi_s)-cos(phi_s+x))*vw;% 
if x~=0
    
    tmp=(lambda-sin(phi_s)+sin(phi_s+x)-x*cos(phi_s+x))/x*vw;% 
    v_i=tmp/sqrt(1+tmp*tmp/3/1E8/3/1E8);
else
    v_i=0;
end
ewi=(1/sqrt(1-v_i/3.0/1E8*v_i/3.0/1E8)-1)*511000;
seypeak = seYPeakTS_m * (1 + seTOneTS_m * (1.0 - power(cosTheta, seTTwoTS_m))); %formula III.E (48a)
seepeak = seEPeakTS_m * (1 + seTThreeTS_m * (1.0 - power(cosTheta, seTFourTS_m))); %formula III.E (48b)
tmpx = ewi / seepeak;
tmpD = seSTS_m * tmpx / (seSTS_m - 1 + power(tmpx, seSTS_m)); %formula III.D (32)
retts = seypeak * tmpD; %formula III.D (31)
tmpdr = power(ewi / seERed_m, seR_m);
retr = sePRed_m * (1.0 - exp(-1 * tmpdr)); %formula III.D (28)
retr = retr * (1.0 + seROne_m * (1.0 - power(cosTheta, seRTwo_m))); %formula III.E (47b)
tmpde = power(abs(ewi - seEScatPeak_m) / seW_m, seP_m) / seP_m;
rete = sePScat_m + (sePScatPeak_m - sePScat_m) * exp(-1 * tmpde); %formula III.D (25)
rete = rete * (1.0 + seEOne_m * (1.0 - power(cosTheta, seETwo_m))); %formula III.E (47a)


Y=retts+retr+rete;

% Y=2.506*exp(-0.0004945*ewi)-1.975*exp(-0.01056*ewi);

        






