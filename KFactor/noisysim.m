function [v, v2, r, r2, rdot, rdot2] = noisysim(x,f,Gt,M,X,PT,GT,GR,R,sigma,index,k,z,z_prev,T,time1,multi)
% ============================================= Magnitude =================================================

lambda = 3*10^8/f; lambda2 = 3*10^8/(f+0.1);
xcoord = z(1,k)-x(1);  ycoord = z(2,k)-x(2);   zcoord = z(3,k)-x(3);

% ============================================= Magnitude =================================================
direct = sqrt((xcoord)^2+(ycoord)^2+(zcoord)^2);
mul    = 1/direct*exp(-i*2*pi*direct/lambda);
mul2   = 1/direct*exp(-i*2*pi*direct/lambda2);

direct_prev =  sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2+(z_prev(3,k)-x(3))^2);
mul_prev    = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda);

drefl1 = sqrt((x(1) + 2)^2 + ((x(2) - z(2,k))/2)^2 + ((x(3) - z(3,k))/2)^2) + sqrt((z(1,k) + 2)^2 + ((x(2) - z(2,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef1  = 1/drefl1*exp(-i*2*pi*drefl1/lambda);
mulRef12 = 1/drefl1*exp(-i*2*pi*drefl1/lambda2);

drefl2 = sqrt((3- x(1))^2 + ((x(2) - z(2,k))/2)^2 + ((x(3) - z(3,k))/2)^2) + sqrt((3 - z(1,k))^2 + ((x(2) - z(2,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef2 = 1/drefl2*exp(-i*2*pi*drefl2/lambda);

drefl3 = sqrt((x(2) + 2)^2 + ((x(1) - z(1,k))/2)^2 + ((x(3) - z(3,k))/2)^2) + sqrt((z(2,k) + 2)^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef3 = 1/drefl3*exp(-i*2*pi*drefl3/lambda);

drefl4 = sqrt((3 - x(2))^2 + ((x(1) - z(1,k))/2)^2 + ((x(3) - z(3,k))/2)^2) + sqrt((3 - z(2,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef4 = 1/drefl4*exp(-i*2*pi*drefl4/lambda);

drefl5 = sqrt((x(3) + 1)^2 + ((x(2) - z(2,k))/2)^2 + ((x(1) - z(1,k))/2)^2) + sqrt((z(3,k) + 1)^2 + ((x(2) - z(2,k))/2)^2+ ((x(1) - z(1,k))/2)^2);%direct + lambda*50;
mulRef5 = 1/drefl5*exp(-i*2*pi*drefl5/lambda);

drefl6 = sqrt((5 - x(3))^2 + ((x(1) - z(1,k))/2)^2 + ((x(2) - z(2,k))/2)^2) + sqrt((5 - z(3,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(2) - z(2,k))/2)^2);%direct + lambda*50;
mulRef6 = 1/drefl6*exp(-i*2*pi*drefl6/lambda);


H  = abs(sqrt(2*R*PT*GT*GR*Gt^2*lambda^2*X^2*M)/(4*pi)*(mul+ 0.1*mulRef1 + 0.1*mulRef2 + 0.1*mulRef3 + 0.1*mulRef4 + 0.1*mulRef5 +0.1*mulRef6));%  + 0.6*mulSum3 %mul+ 0.25* mulSum + 0.25*mulSum2
H2 = abs(sqrt(2*R*PT*GT*GR*Gt^2*lambda^2*X^2*M)/(4*pi)*(mul));%mul+ 0.25* mulSum + 0.25*mulSum2

% ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++

%pd = makedist('Rician','s',sqrt(H),'sigma',0.01*sqrt(H));
 
v = H;%random(pd).^2;
v2 = H2;

unirand = rand;

r  = sqrt(2*PT*GT*GR*Gt^2*lambda^2*X^2*M*R)/(4*pi*H);
r2 = sqrt(2*PT*GT*GR*Gt^2*lambda^2*X^2*M*R)/(4*pi*H2);
%%

% % ============================================== Phase ====================================================

phi_prev_mod = angle(sqrt(2*R*PT*GT*GR*Gt^2*lambda^2*X^2*M)/(4*pi)*(mul_prev));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod      = angle(sqrt(2*R*PT*GT*GR*Gt^2*lambda^2*X^2*M)/(4*pi)*(mul + 0.2*mulRef1));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
phi_mod2     = angle(sqrt(2*R*PT*GT*GR*Gt^2*lambda2^2*X^2*M)/(4*pi)*(mul2 + 0.2*mulRef12));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

phi_mod_gt      = angle(sqrt(2*R*PT*GT*GR*Gt^2*lambda^2*X^2*M)/(4*pi)*(mul))%  + 0.2*mulRef1));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
phi_mod2_gt     = angle(sqrt(2*R*PT*GT*GR*Gt^2*lambda2^2*X^2*M)/(4*pi)*(mul2))%+ 0.2*mulRef12));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);


diff = phi_mod - phi_prev_mod;

delta_phi = diff;

% +++++++++++++++++++++++++++++++++++++++++++++ Noise of Phase +++++++++++++++++++++++++++++++++++++++++++++++

new_sigma = sigma/phi_mod + sigma/phi_prev_mod;
phi_noise = 800000*new_sigma.*rand;

delta_phi2 = exp(phi_noise)*delta_phi;


% rdot2 = lambda/(4*pi)*delta_phi2*1/T;

rdot  = 2*3*10^8/(4*pi)*(phi_mod -phi_mod2)/(0.1);%lambda/(4*pi)*phi_mod*1/T;
rdot2 = lambda/(4*pi)*diff*1/T;

if abs(rdot2) > 2
    rdot2 = 0;
end

end


% global l; global l1; global l2;global l3;
% 
% if index == 0
%     if (phi_prev_mod - phi_mod ) > 4.5
%         l = l + 1;
%         phi_prev_mu = phi_prev_mod + (l-1)*2*pi;
%         phi_mu = phi_mod + l*2*pi;
%     elseif phi_mod > 5 && phi_prev_mod < 1 
%         l = l - 1;
%         phi_prev_mu = phi_prev_mod + l*2*pi;
%         phi_mu = phi_mod + (l-1)*2*pi;
%     else
%         phi_prev_mu = phi_prev_mod + l*2*pi;
%         phi_mu = phi_mod + l*2*pi;
%     end
%     
% elseif index == 1
%     if (phi_prev_mod - phi_mod)> 1.8 %5 && phi_mod < 1 
%         l1 = l1 + 1;
%         phi_prev_mu = phi_prev_mod + (l1-1)*2*pi;
%         phi_mu = phi_mod + l1*2*pi;
% 
%     elseif (phi_mod - phi_prev_mod)> 2.8 %phi_prev_mod < 1 && phi_mod > 5
%         l1 = l1 - 1;
%         phi_prev_mu = phi_prev_mod + l1*2*pi;
%         phi_mu = phi_mod + (l1-1)*2*pi;
% 
%     else
%         phi_prev_mu = phi_prev_mod + l1*2*pi;
%         phi_mu = phi_mod + l1*2*pi;
%     end
%     
% elseif index == 2
%     if k < 0.47*length(z)
%         if (phi_prev_mod - phi_mod)> 2.9 %5 && phi_mod < 1 
%             l2 = l2 + 1;
%             phi_prev_mu = phi_prev_mod + (l2-1)*2*pi;
%             phi_mu = phi_mod + l2*2*pi;
% 
%         else
%             phi_prev_mu = phi_prev_mod + l2*2*pi;
%             phi_mu = phi_mod + l2*2*pi;
%         end
%     else
%         
%         if (phi_mod - phi_prev_mod)> 3.0 %phi_prev_mod < 1 && phi_mod > 5
%             l2 = l2 - 1;
%             phi_prev_mu = phi_prev_mod + l2*2*pi;
%             phi_mu = phi_mod + (l2-1)*2*pi;
% 
%         else
%             phi_prev_mu = phi_prev_mod + l2*2*pi;
%             phi_mu = phi_mod + l2*2*pi;
%         end
%     end
%     
% else 
%     if k < 0.5*length(z)
%         if (phi_mod - phi_prev_mod)> 2.3 %phi_prev_mod < 1 && phi_mod > 5
%             l3 = l3 - 1;
%             phi_prev_mu = phi_prev_mod + l3*2*pi;
%             phi_mu = phi_mod + (l3-1)*2*pi;
%             
%         elseif (phi_prev_mod - phi_mod)> 5.0 %5 && phi_mod < 1 
%             l3 = l3 + 1;
%             phi_prev_mu = phi_prev_mod + (l3-1)*2*pi;
%             phi_mu = phi_mod + l3*2*pi;
%         else
%             phi_prev_mu = phi_prev_mod + l3*2*pi;
%             phi_mu = phi_mod + l3*2*pi;
%         end
%         
%     else
%         if (phi_prev_mod - phi_mod)> 3.0%5 && phi_mod < 1 
%             l3 = l3 + 1;
%             phi_prev_mu = phi_prev_mod + (l3-1)*2*pi;
%             phi_mu = phi_mod + l3*2*pi;
%    
%         else
%             phi_prev_mu = phi_prev_mod + l3*2*pi;
%             phi_mu = phi_mod + l3*2*pi;
%         end
%                
%     end
% end

%phi_prev_conc = phi_prev_mu;  %icdf('Normal',uniphase,phi_prev_mu,sigma);
%phi_conc      = phi_mu;       %icdf('Normal',uniphase2,phi_mu,sigma);