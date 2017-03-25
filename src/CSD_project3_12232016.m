close all;
clear;

s=tf('s');
G = (s + 0.1541) / (s*((s^2) + 0.739*s + 0.921));
opentf = 1.151*G;
oriClosedLoop = G/(1+G);

t = 0:0.1:30;
%% a
%steady state error
syms s
a_G = (s + 0.1541) / (s*((s^2) + 0.739*s + 0.921));
a_Kv = limit((s*a_G));
a_ess = 1/a_Kv;

%ramp response
a_y = lsim(oriClosedLoop.Numerator{1},oriClosedLoop.Denominator{1},t,t);
figure(1),plot(t,t,'b',t,a_y,'r')
xlabel('Time(sec)s');
ylabel('Amplitude');
legend('Input','Output');
title('Input: blue, Output: red');


%Refer :http://ctms.engin.umich.edu/CTMS/index.php?aux=Extras_Ess
%Refer :http://digitalcommons.calpoly.edu/cgi/viewcontent.cgi?article=1452&context=theses

%% b
b_ess = 0.1;
b_Kv = 1/b_ess;
b_zeros = roots(G.Numerator{1});
b_poles = roots(G.Denominator{1});
b_K  = (-(10*prod(b_poles(2:3)))/(b_zeros));

%corresponding closed-loop system
b_Closed = b_K*G / (1 + (b_K*G));

%ramp response
b_y = lsim(b_Closed.Numerator{1}, b_Closed.Denominator{1},t,t);
figure(2),plot(t,t,'b',t,b_y,'r')
xlabel('Time(sec)');
ylabel('Amplitude');
legend('Input','Output');
title('Input: blue, Output: red');

%% c
w = 0.1:0.1:100;
num = [1 0.1541];
den = [1 0.739 0.921];
figure(3);
subplot(1,2,1);
bode(b_K*num,den,w);
title('BodePlot: K*G');
[Gm_act,Pm_act,Wcg_act,Wcp_act] = margin(num,den);
phiMax = 50 - Pm_act + 10;
phiRad = (phiMax/180)*pi;
a = (1 + sin(phiRad))/(1 - sin(phiRad));
delta_G = 20 * log10(a);
subplot(1,2,2);
hold on
bode(num,den);
title('Bode Plot of corresponding closed-loop system ramp');
hold off

%constant from the book.
w_max = 10.5;
pc = w_max*sqrt(a);
numc = [a pc];
denc = [1 pc];
hold on;
bode(numc,denc);
hold off;
compNum = conv(num, numc);
compDen = conv(den, denc);
hold on;
bode(compNum, compDen);
hold off;
legend('G','pc','Compensation');

%finding closed-loop transfer function
[closedNum, closedDen] = feedback(G.Numerator{1}, G.Denominator{1},1,1,-1);
[closedNumComp, closedDenComp] = feedback(compNum,compDen,1,1,-1);

%%comments: in step response, the compensation, rising time and settling time
%%are reduced significantly. 
figure(4);
ystep=step(closedNum,closedDen,t);
ystepComp_Lead = step(closedNumComp, closedDenComp,t);
plot(t,ystep,'b',t,ystepComp_Lead,'r');
%axis([0 2 0 0.5])
title('Step Response: uncompensation vs compensation');
legend('uncompensation','compensation');
s=tf('s');
Gcs_lead= (a*s + pc)/(s+pc)

%% d
Wcg_new= 2.1;

% dl=b_K;
% g1 = abs(j*Wcg_new);
% g2 = abs(j*Wcg_new + b_poles(2));
% g3 = abs(j*Wcg_new + b_poles(3));
% dG = dl/(g1*g2*g3)
% mag(dG)
delta_G_new = 30.4;
a_new= 10^(delta_G_new/(-20));
zc_new = Wcg_new/10;
pc_new = a_new*zc_new;
numlag = [1 zc_new];
denlag = [1 pc_new];
newnum= conv(num,numlag);
newden = conv(den,denlag);
Gcs_lag = (a*s +pc_new)/(s+pc_new)
figure(5)
bode(Gcs_lag.Numerator{1},Gcs_lag.Denominator{1})
ystepComp_Lag = step(newnum,newden,t);
figure(6)
plot(t,ystepComp_Lead,'b',t,ystepComp_Lag,'r');
%axis([0 2 0 0.5])
title('Step Response: Lead vs Lag');
legend('Phase Lead','Phase Lag');

[closedLagNum, closedLagDen] = feedback(newnum,newden,1,1,1);

%http://ece.gmu.edu/~gbeale/ece_421/comp_freq_lag.pdf
figure(7)
rlocus(num,den);
hold on
rlocus(closedNumComp,closedDenComp);