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

for k = 1:1:length(time)-1  
[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1), r_sim1__(k+1),rdot_sim1(k+1), rdot_sim1_(k+1),phi1(k+1), phi21(k+1), phi31(k+1), phi41(k+1), offset11, offset21, offset31, offset41] = noisysim(x1,f1,Gt,M,X,PT,0.7*GT, GR,R,sigma,k,z,z_prev,phi1(k),phi31(k),time(k+1)-time(k),magD12(k), offset11, offset21, offset31, offset41);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1), r_sim2__(k+1),rdot_sim2(k+1), rdot_sim2_(k+1),phi2(k+1), phi22(k+1), phi32(k+1), phi42(k+1), offset12, offset22, offset32, offset42] = noisysim(x2,f2,Gt,M,X,PT,  7*GT, GR,R,sigma,k,z,z_prev,phi2(k),phi32(k),time(k+1)-time(k),magD22(k), offset12, offset22, offset32, offset42);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1), r_sim3__(k+1),rdot_sim3(k+1), rdot_sim3_(k+1),phi3(k+1), phi23(k+1), phi33(k+1), phi43(k+1), offset13, offset23, offset33, offset43] = noisysim(x3,f3,Gt,M,X,PT,    GT, GR,R,sigma,k,z,z_prev,phi3(k),phi33(k),time(k+1)-time(k),magD32(k), offset13, offset23, offset33, offset43);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1), r_sim4__(k+1),rdot_sim4(k+1), rdot_sim4_(k+1),phi4(k+1), phi24(k+1), phi34(k+1), phi44(k+1), offset14, offset24, offset34, offset44] = noisysim(x4,f4,Gt,M,X,PT,0.5*GT, GR,R,sigma,k,z,z_prev,phi4(k),phi34(k),time(k+1)-time(k),magD42(k), offset14, offset24, offset34, offset44);  
end
%%
figure
subplot(411), plot(time, radial(1,:))
subplot(412), plot(time, radial(2,:))
subplot(413), plot(time, radial(3,:))
subplot(414), plot(time, radial(4,:))

%%
figure
subplot(411), plot(time, 3/(4*pi)*(phi31-phi41))
subplot(412), plot(time, 3/(4*pi)*(phi32-phi42))
subplot(413), plot(time, 3/(4*pi)*(phi33-phi43))
subplot(414), plot(time, 3/(4*pi)*(phi34-phi44))
%%
[rms(r_sim1(yST:yET) - radial(1,yST:yET)), rms(r_sim2(yST:yET) - radial(2,yST:yET)), rms(r_sim3(yST:yET) - radial(3,yST:yET))]
[rms(rdot_sim1(yST:yET) - radial(1,yST:yET)), rms(rdot_sim2(yST:yET) - radial(2,yST:yET)), rms(rdot_sim3(yST:yET) - radial(3,yST:yET))]
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

%%
figure
subplot(4,1,1),plot(time, r_sim1,'LineWidth',2);hold on; plot(time, r_sim1_,'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time, r_sim2,'LineWidth',2);hold on; plot(time, r_sim2_,'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time, r_sim3,'LineWidth',2);hold on; plot(time, r_sim3_,'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time, r_sim4,'LineWidth',2);hold on; plot(time, r_sim4_,'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
%
figure
subplot(4,1,1),plot(time, rdot_sim1_,'LineWidth',3);%hold on; plot(time, radial(1,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,2),plot(time, rdot_sim2_,'LineWidth',3);%hold on; plot(time, radial(2,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,3),plot(time, rdot_sim3_,'LineWidth',3);%hold on; plot(time, radial(3,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,4),plot(time, rdot_sim4_,'LineWidth',3);%hold on; plot(time, radial(4,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
%%
figure
subplot(4,1,1),plot(time, rdot_sim1,'LineWidth',3);%hold on; plot(time, radial(1,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,2),plot(time, rdot_sim2,'LineWidth',3);%hold on; plot(time, radial(2,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,3),plot(time, rdot_sim3,'LineWidth',3);%hold on; plot(time, radial(3,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,4),plot(time, rdot_sim4,'LineWidth',3);%hold on; plot(time, radial(4,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')

%%
figure
subplot(4,1,1),plot(time, r_sim1_,'LineWidth',3);hold on; plot(time, radial(1,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,2),plot(time, r_sim2_,'LineWidth',3);hold on; plot(time, radial(2,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,3),plot(time, r_sim3_,'LineWidth',3);hold on; plot(time, radial(3,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,4),plot(time, r_sim4_,'LineWidth',3);hold on; plot(time, radial(4,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')

%

% 
% [-mean(H1(2:124840) - magD12(2:124840)), rms(H1(2:124840) - magD12(2:124840)), sqrt(rms(H1(2:124840) - magD12(2:124840)))]
% [-mean(H2(2:124840) - magD22(2:124840)), rms(H2(2:124840) - magD22(2:124840)), sqrt(rms(H2(2:124840) - magD22(2:124840)))]
% [-mean(H3(2:124840) - magD32(2:124840)), rms(H3(2:124840) - magD32(2:124840)), sqrt(rms(H3(2:124840) - magD32(2:124840)))]
% [-mean(H4(2:124840) - magD42(2:124840)), rms(H4(2:124840) - magD42(2:124840)), sqrt(rms(H4(2:124840) - magD42(2:124840)))]
% 
% 
% [-mean(H1_(2:124840) - magD12(2:124840)), rms(H1_(2:124840) - magD12(2:124840)), sqrt(rms(H1_(2:124840) - magD12(2:124840)))]
% [-mean(H2_(2:124840) - magD22(2:124840)), rms(H2_(2:124840) - magD22(2:124840)), sqrt(rms(H2_(2:124840) - magD22(2:124840)))]
% [-mean(H3_(2:124840) - magD32(2:124840)), rms(H3_(2:124840) - magD32(2:124840)), sqrt(rms(H3_(2:124840) - magD32(2:124840)))]
% [-mean(H4_(2:124840) - magD42(2:124840)), rms(H4_(2:124840) - magD42(2:124840)), sqrt(rms(H4_(2:124840) - magD42(2:124840)))]
% 
% [-mean(H1_(2:124840) - H1(2:124840)), rms(H1_(2:124840) - H1(2:124840)), sqrt(rms(H1_(2:124840) - H1(2:124840)))]
% [-mean(H2_(2:124840) - H2(2:124840)), rms(H2_(2:124840) - H2(2:124840)), sqrt(rms(H2_(2:124840) - H2(2:124840)))]
% [-mean(H3_(2:124840) - H3(2:124840)), rms(H3_(2:124840) - H3(2:124840)), sqrt(rms(H3_(2:124840) - H3(2:124840)))]
% [-mean(H4_(2:124840) - H4(2:124840)), rms(H4_(2:124840) - H4(2:124840)), sqrt(rms(H4_(2:124840) - H4(2:124840)))]

% --------------------------------------------------------- End ==========================================================================

% Parameters of reader
% PT = 1;         % reader's transmitted power
% GT = sqrt(1.462);     % reader's trasmitter antenna gain 9.5dBi
% GR = sqrt(1.462);     % reader's receiver   antenna gain 9.5dBi
% R = 15;
% 
% Gt = 0.2193;    % tag's antenna gain
% X  = 0.35;          % polarization mismatch
% M  = 4;             % load modulation factor of the tag
% f1 = 5.8*10^9;
% f2 = 5.83*10^9;
% f3 = 5.82*10^9;
% f4 = 5.85*10^9;

% 
% figure
% subplot(411), plot(time1, phaseD1);
% subplot(412), plot(time2, phaseD2);
% subplot(413), plot(time3, phaseD3);
% subplot(414), plot(time4, phaseD4);
% 
% figure
% subplot(411), plot(time1, magD1)
% subplot(412), plot(time2, magD2)
% subplot(413), plot(time3, magD3)
% subplot(414), plot(time4, magD4)
%
% figure
% subplot(411), plot(time1+13.3287, phaseD1,'LineWidth',2), xlim([14, 145]), hold on, plot(time, rdot_sim1+1.57,'LineWidth',1), xlim([14, 145]), legend('Measurement','Simulation','location','BestOutside'), title('Phase from Reader $\#1$ in 3D','interpreter','latex'), ylabel('Phase [rad]'), xlabel('time [s]')
% subplot(412), plot(time2+13.3287, phaseD2,'LineWidth',2), xlim([14, 145]), hold on, plot(time, rdot_sim2+1.57,'LineWidth',1), xlim([14, 145]), legend('Measurement','Simulation','location','BestOutside'), title('Phase from Reader $\#2$ in 3D','interpreter','latex'), ylabel('Phase [rad]'), xlabel('time [s]')
% subplot(413), plot(time3+13.6187, phaseD3,'LineWidth',2), xlim([14, 145]), hold on, plot(time, rdot_sim3+1.56,'LineWidth',1), xlim([14, 145]), legend('Measurement','Simulation','location','BestOutside'), title('Phase from Reader $\#3$ in 3D','interpreter','latex'), ylabel('Phase [rad]'), xlabel('time [s]')
% subplot(414), plot(time4+13.3287, phaseD4,'LineWidth',2), xlim([14, 145]), hold on, plot(time, rdot_sim4+1.55,'LineWidth',1), xlim([14, 145]), legend('Measurement','Simulation','location','BestOutside'), title('Phase from Reader $\#4$ in 3D','interpreter','latex'), ylabel('Phase [rad]'), xlabel('time [s]')

% mag1   = magD1(145388:315021);     mag2 = magD2(145388:315021);     mag3 = magD3(295791:640912);     mag4 = magD4(145388:315021);
% phase1 = phaseD1(145388:315021); phase2 = phaseD2(145388:315021); phase3 = phaseD3(295791:640912); phase4 = phaseD4(145388:315021);
% t1 = time1(145388:315021);           t2 = time2(145388:315021);       t3 = time3(295791:640912);      t4 = time4(145388:315021); 
% 
% %
% magD12 = NaN(1,length(time)); magD22 = NaN(1,length(time)); magD32 = NaN(1,length(time)); magD42 = NaN(1,length(time));
% magD12(1) = magD1(1); magD22(1) = magD1(1);  magD32(1) = magD1(1);  magD42(1) = magD1(1); 
% indexT = 2;
% 
% for indexMag = 1:1:length(time1)
% %     timeT = time(indexT)
% %     time1T = time1(indexMag)+13.3287
%     if abs(time1(indexMag)+13.3287 - time(indexT)) < 0.0002 || ((time1(indexMag)+13.3287) > time(indexT))
%         %indexT
%         if indexT  <= length(time)
%             magD12(indexT) = magD1(indexMag);
%             magD22(indexT) = magD2(indexMag);
%             
%             if indexMag < length(time4)
%                 magD42(indexT) = magD4(indexMag);
%             end
%         end
%         
%         indexT = indexT + 1;
%     end
% end
% 
% indexT = 1;
% 
% for indexMag = 1:1:length(time3)
%     if abs(time3(indexMag)+13.3287 - time(indexT)) < 0.0002 || ((time3(indexMag)+13.3287) > time(indexT))
%         
%         if indexT  <= length(time)
%             magD32(indexT) = magD3(indexMag);
%         end
%         
%         indexT = indexT + 1;
%     end
% end