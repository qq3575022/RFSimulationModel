clc, clear, close all
% =================================================================== Load Data ======================================================================
% ----------- Position of Four Readers ---------
x1 = [0,    0,    0.865];  
x2 = [2.29, 0,    1.27];   
x3 = [2.29, 2.52, 0.865]; 
x4 = [0,    2.52, 1.27];

% -------------- Time and Coordinates ----------
% time:     IMU Measurement
% coord3:   Ground Truth of 3D Coordinates
% z,z_prev: 3D coordinates for Simulation
[time, coord3, radial, z, z_prev, H1, H2, H3, H4, H1_, H2_, H3_, H4_, r_sim, r_sim1, r_sim2, r_sim3, r_sim4, r_sim1_, r_sim2_, r_sim3_, r_sim4_, r_meas, r_meas1, r_meas2, r_meas3, r_meas4, rphase, rphase1, rphase2, rphase3, rphase4, rdot_sim1, rdot_sim2, rdot_sim3, rdot_sim4, rdot_sim1_, rdot_sim2_, rdot_sim3_, rdot_sim4_, phigt1_1, phigt1_2, phigt1_3, phigt1_4, phigt2_1, phigt2_2, phigt2_3, phigt2_4, phi1_1, phi1_2, phi1_3, phi1_4, phi2_1, phi2_2, phi2_3, phi2_4] = get3Dcoord(x1, x2, x3, x4);

% -------------- Measurement Magnitude ----------
[magD12, magD22, magD32, magD42] = getMeas(time);% Measurement Magnitude of Length 130854

% Start and End Index of 3D Motion
yST  = find(abs(time-107.99)<0.002);   yST = yST(1)-1;
yET  = find(abs(time-111.984)<0.002);  yET = yET(1)+1;

% =================================================================== Simulation ======================================================================

% Parameters of reader
Gt = 0.6*1.462*3/4;    % tag's antenna gain
X  = 0.85;             % polarization mismatch
M  = sqrt(4);          % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

% Parameters of reader
PT = 1;         % reader's transmitted power
GT = sqrt(1.462);     % reader's trasmitter antenna gain 9.5dBi
GR = sqrt(1.462);     % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.0000012; 

offset11 = 0; offset21 = 0; offset31 = 0; offset41 = 0;offset12 = 0; offset22 = 0; offset32 = 0; offset42 = 0;
offset13 = 0; offset23 = 0; offset33 = 0; offset43 = 0;offset14 = 0; offset24 = 0; offset34 = 0; offset44 = 0;

for k = yST-1:1:yET  

[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1),r_meas1(k+1),rphase1(k+1),rdot_sim1(k+1),rdot_sim1_(k+1),phigt1_1(k+1),phigt2_1(k+1),phi1_1(k+1),phi2_1(k+1),offset11,offset21,offset31,offset41] = noisysimNoMultiStatic(x1,f1,Gt,M,X,PT,0.7*GT,GR,R,sigma,k,z,z_prev,phigt1_1(k),phi1_1(k),time(k+1)-time(k),magD12(k),offset11,offset21,offset31,offset41);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1),r_meas2(k+1),rphase2(k+1),rdot_sim2(k+1),rdot_sim2_(k+1),phigt1_2(k+1),phigt2_2(k+1),phi1_2(k+1),phi2_2(k+1),offset12,offset22,offset32,offset42] = noisysimNoMultiStatic(x2,f2,Gt,M,X,PT,  7*GT,GR,R,sigma,k,z,z_prev,phigt1_2(k),phi1_2(k),time(k+1)-time(k),magD22(k),offset12,offset22,offset32,offset42);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1),r_meas3(k+1),rphase3(k+1),rdot_sim3(k+1),rdot_sim3_(k+1),phigt1_3(k+1),phigt2_3(k+1),phi1_3(k+1),phi2_3(k+1),offset13,offset23,offset33,offset43] = noisysimNoMultiStatic(x3,f3,Gt,M,X,PT,    GT,GR,R,sigma,k,z,z_prev,phigt1_3(k),phi1_3(k),time(k+1)-time(k),magD32(k),offset13,offset23,offset33,offset43);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1),r_meas4(k+1),rphase4(k+1),rdot_sim4(k+1),rdot_sim4_(k+1),phigt1_4(k+1),phigt2_4(k+1),phi1_4(k+1),phi2_4(k+1),offset14,offset24,offset34,offset44] = noisysimNoMultiStatic(x4,f4,Gt,M,X,PT,0.5*GT,GR,R,sigma,k,z,z_prev,phigt1_4(k),phi1_4(k),time(k+1)-time(k),magD42(k),offset14,offset24,offset34,offset44);  

r_sim(:,k+1)  = lsqnonlin(@(xx)getMulPath(r_sim1(k+1), r_sim2(k+1), r_sim3(k+1),  r_sim4(k+1), x1, x2, x3, x4, xx), [0,0,0]);
rphase(:,k+1) = lsqnonlin(@(xx)getMulPath(rphase1(k+1),rphase2(k+1),rphase3(k+1),rphase4(k+1), x1, x2, x3, x4, xx), [0,0,0]);
r_meas(:,k+1) = lsqnonlin(@(xx)getMulPath(r_meas1(k+1),r_meas2(k+1),r_meas3(k+1),r_meas4(k+1), x1, x2, x3, x4, xx), [0,0,0]);

end
%%
figure
subplot(411),plot(H1_(yST:yET)); hold on; plot(magD12(yST:yET)+0.0002284);
subplot(412),plot(H2_(yST:yET)); hold on; plot(magD22(yST:yET)+0.0006892);
subplot(413),plot(H3_(yST:yET)); hold on; plot(magD32(yST:yET));
subplot(414),plot(H4_(yST:yET)); hold on; plot(magD42(yST:yET)+0.00010078);
%% ==================================================================== Results ======================================================================
% -------------- Radial Distances RMS: Reader1 Reader2 Reader3 Reader4 ----------
rmeas  = [rms(r_meas1(yST:yET) - radial(1,yST:yET)), rms(r_meas2(yST:yET) - radial(2,yST:yET)), rms(r_meas3(yST:yET) - radial(3,yST:yET)), rms(r_meas4(yST:yET)  - radial(4,yST:yET))]
rmag   = [rms(r_sim1(yST:yET) - radial(1,yST:yET)),  rms(r_sim2(yST:yET) - radial(2,yST:yET)),  rms(r_sim3(yST:yET) - radial(3,yST:yET)),  rms(r_sim4(yST:yET)  - radial(4,yST:yET))]
rphase = [rms(rphase1(yST:yET) - radial(1,yST:yET)), rms(rphase2(yST:yET) - radial(2,yST:yET)), rms(rphase3(yST:yET) - radial(3,yST:yET)), rms(rphase4(yST:yET) - radial(4,yST:yET))]

meas  = sqrt(rmeas(1)^2  + rmeas(2)^2  + rmeas(3)^2)
mag   = sqrt(rmag(1)^2   + rmag(2)^2   + rmag(3)^2)
phase = sqrt(rphase(1)^2 + rphase(2)^2 + rphase(3)^2)
%%
figure
subplot(411), plot(time, 3/(4*pi)*(phi1_1-phi2_1))
subplot(412), plot(time, 3/(4*pi)*(phi1_2-phi2_2))
subplot(413), plot(time, 3/(4*pi)*(phi1_3-phi2_3))
subplot(414), plot(time, 3/(4*pi)*(phi1_4-phi2_4))

% -------------- Measurement Radial Distances: Reader1 Reader2 Reader3 Reader4 ----------

figure
subplot(4,1,1),plot(time(yST:yET), r_meas1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), r_meas2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), r_meas3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), r_meas4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

% -------------- Magnitude Radial Distances: Reader1 Reader2 Reader3 Reader4 ----------

figure
subplot(4,1,1),plot(time(yST:yET), r_sim1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), r_sim2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), r_sim3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), r_sim4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

% -------------- Two Frequencies Radial Distances: Reader1 Reader2 Reader3 Reader4 ----------

figure
subplot(4,1,1),plot(time(yST:yET), rphase1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), rphase2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), rphase3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), rphase4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
