clc, clear, close all
[time, coord3, z, z_prev, H1, H2, H3, H4, H1_, H2_, H3_, H4_, r_sim1, r_sim2, r_sim3, r_sim4, r_sim1_, r_sim2_, r_sim3_, r_sim4_, rdot_sim1, rdot_sim2, rdot_sim3, rdot_sim4, rdot_sim1_, rdot_sim2_, rdot_sim3_, rdot_sim4_, phi1, phi2, phi3, phi4, phi21, phi22, phi23, phi24, phi31, phi32, phi33, phi34, phi41, phi42, phi43, phi44] = get3Dcoord();

[magD12, magD22, magD32, magD42] = getMeas(time);


yST  = find(abs(time-107.99)<0.002);   yST = yST(1)-1;
yET  = find(abs(time-111.984)<0.002);  yET = yET(1)+1;

x1 = [0,    0,    0.865];  
x2 = [2.29, 0,    1.27];   
x3 = [2.29, 2.52, 0.865]; 
x4 = [0,    2.52, 1.27];

radial = NaN(4, length(coord3));
for i = 1:1:length(coord3)
radial(1,i) = sqrt((coord3(1,i) - x1(1))^2+(coord3(2,i) - x1(2))^2+(coord3(3,i) - x1(3))^2);
radial(2,i) = sqrt((coord3(1,i) - x2(1))^2+(coord3(2,i) - x2(2))^2+(coord3(3,i) - x2(3))^2);
radial(3,i) = sqrt((coord3(1,i) - x3(1))^2+(coord3(2,i) - x3(2))^2+(coord3(3,i) - x3(3))^2);
radial(4,i) = sqrt((coord3(1,i) - x4(1))^2+(coord3(2,i) - x4(2))^2+(coord3(3,i) - x4(3))^2);
end

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

r_sim = NaN(3,length(coord3)); rdot_sim = NaN(3,length(coord3)); measurement = NaN(3,length(coord3));
for k = yST-1:1:yET
[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1), r_sim1__(k+1),rdot_sim1(k+1), rdot_sim1_(k+1),phi1(k+1), phi21(k+1), phi31(k+1), phi41(k+1), offset11, offset21, offset31, offset41] = noisysim(x1,f1,Gt,M,X,PT,0.7*GT, GR,R,sigma,k,z,z_prev,phi1(k),phi31(k),time(k+1)-time(k),magD12(k), offset11, offset21, offset31, offset41);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1), r_sim2__(k+1),rdot_sim2(k+1), rdot_sim2_(k+1),phi2(k+1), phi22(k+1), phi32(k+1), phi42(k+1), offset12, offset22, offset32, offset42] = noisysim(x2,f2,Gt,M,X,PT,  7*GT, GR,R,sigma,k,z,z_prev,phi2(k),phi32(k),time(k+1)-time(k),magD22(k), offset12, offset22, offset32, offset42);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1), r_sim3__(k+1),rdot_sim3(k+1), rdot_sim3_(k+1),phi3(k+1), phi23(k+1), phi33(k+1), phi43(k+1), offset13, offset23, offset33, offset43] = noisysim(x3,f3,Gt,M,X,PT,    GT, GR,R,sigma,k,z,z_prev,phi3(k),phi33(k),time(k+1)-time(k),magD32(k), offset13, offset23, offset33, offset43);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1), r_sim4__(k+1),rdot_sim4(k+1), rdot_sim4_(k+1),phi4(k+1), phi24(k+1), phi34(k+1), phi44(k+1), offset14, offset24, offset34, offset44] = noisysim(x4,f4,Gt,M,X,PT,0.5*GT, GR,R,sigma,k,z,z_prev,phi4(k),phi34(k),time(k+1)-time(k),magD42(k), offset14, offset24, offset34, offset44);  
% [r_sim1(k+1), r_sim2(k+1), r_sim3(k+1), r_sim4(k+1)];
% [rdot_sim1(k+1), rdot_sim2(k+1), rdot_sim3(k+1), rdot_sim4(k+1)];
r_sim(:,k+1) = lsqnonlin(@(xx)getMulPath(r_sim1(k+1), r_sim2(k+1), r_sim3(k+1), r_sim4(k+1), x1, x2, x3, x4, xx), [0,0,0]);
rdot_sim(:,k+1) = lsqnonlin(@(xx)getMulPath(rdot_sim1(k+1), rdot_sim2(k+1), rdot_sim3(k+1), rdot_sim4(k+1), x1, x2, x3, x4, xx), [0,0,0]);
measurement(:,k+1) = lsqnonlin(@(xx)getMulPath(r_sim1__(k+1), r_sim2__(k+1), r_sim3__(k+1), r_sim4__(k+1), x1, x2, x3, x4, xx), [0,0,0]);
end
%%
meas = [rms(measurement(1, yST:yET) - coord3(1,yST:yET)),    rms(measurement(2, yST:yET) - coord3(2,yST:yET)),    rms(measurement(3, yST:yET) - coord3(3,yST:yET))]
rsim = [rms(r_sim(1, yST:yET) - coord3(1,yST:yET)),    rms(r_sim(2, yST:yET) - coord3(2,yST:yET)),    rms(r_sim(3, yST:yET) - coord3(3,yST:yET))]
rdot = [rms(rdot_sim(1, yST:yET)+0.4 - coord3(1,yST:yET)), rms(rdot_sim(2, yST:yET)+0.5 - coord3(2,yST:yET)), rms(rdot_sim(3, yST:yET)+0.2 - coord3(3,yST:yET))]

sqrt(meas(1)^2 + meas(2)^2 + meas(3)^2)
sqrt(rsim(1)^2 + rsim(2)^2 + rsim(3)^2)
sqrt(rdot(1)^2 + rdot(2)^2 + rdot(3)^2)
%%
figure
subplot(3,1,1),plot(time(yST:yET), rdot_sim(1, yST:yET)+0.4,'LineWidth',2);hold on; plot(time(yST:yET), coord3(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(3,1,2),plot(time(yST:yET), rdot_sim(2, yST:yET)+0.4,'LineWidth',2);hold on; plot(time(yST:yET), coord3(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(3,1,3),plot(time(yST:yET), rdot_sim(3, yST:yET)+0.2,'LineWidth',2);hold on; plot(time(yST:yET), coord3(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

%%
figure
subplot(411), plot(time(2:end), 1*10^2/(4*pi)*(phi1(1:end-1)-phi21(2:end))), hold on; plot(time(2:end), radial(1,2:end),'LineWidth',2);
subplot(412), plot(time(2:end), 1*10^2/(4*pi)*(phi2(1:end-1)-phi22(2:end))), hold on; plot(time(2:end), radial(2,2:end),'LineWidth',2);
subplot(413), plot(time(2:end), 1*10^2/(4*pi)*(phi3(1:end-1)-phi23(2:end))), hold on; plot(time(2:end), radial(3,2:end),'LineWidth',2);
subplot(414), plot(time(2:end), 1*10^2/(4*pi)*(phi4(1:end-1)-phi24(2:end))), hold on; plot(time(2:end), radial(4,2:end),'LineWidth',2);
%%
meas = [rms(r_sim1__(yST:yET)+1.2 - radial(1,yST:yET)),    rms(r_sim2__(yST:yET)+1.2 - radial(2,yST:yET)),    rms(r_sim3__(yST:yET) - radial(3,yST:yET)),    rms(r_sim4__(yST:yET)+0.6 - radial(4,yST:yET))]
mag = [rms(r_sim1(yST:yET) - radial(1,yST:yET)),    rms(r_sim2(yST:yET) - radial(2,yST:yET)),    rms(r_sim3(yST:yET) - radial(3,yST:yET)),    rms(r_sim4(yST:yET) - radial(4,yST:yET))]
tof = [rms(rdot_sim1(yST:yET) - radial(1,yST:yET)), rms(rdot_sim2(yST:yET) - radial(2,yST:yET)), rms(rdot_sim3(yST:yET) - radial(3,yST:yET)), rms(rdot_sim4(yST:yET) - radial(4,yST:yET))]
%%
figure
subplot(4,1,1),plot(time(yST:yET), r_sim1__(yST:yET)+1.2,'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), r_sim2__(yST:yET)+1.2,'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), r_sim4__(yST:yET)+0.6,'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
%%
figure
subplot(4,1,1),plot(time(yST:yET), r_sim1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), r_sim2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), r_sim3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), r_sim4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
%%
figure
subplot(4,1,1),plot(time(yST:yET), rdot_sim1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), rdot_sim2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), rdot_sim3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), rdot_sim4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

