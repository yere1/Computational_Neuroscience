clear all
clc

dt=0.001;
tmax=4;
t=dt:dt:tmax;
delta_t=1e-4;
% tau_D=0.25;
 
g_syn=zeros(size(t));
delta_g=1e-9;
tau=100e-3;

% g_sec=zeros(size(t));
% g_max=2e-9;
% tau_sec=100e-3;

C=20e-12;
E_l=-65e-3;
E_syn=0;
g_l=2e-9;
V=zeros(size(t));

r=repelem([20 100 10 50], length(t)/4);
p=r*delta_t;
p_0=0.5;
% D=ones(size(t));
% F=ones(size(t));
% f_fac=0.25;
% tau_fac=0.25;
% F_max=1/p_0;

spikes=rand(size(t))<p;

for ind=2:length(t)
     
   
     D_gsyn=-g_syn(ind-1)/tau;
     g_syn(ind)=g_syn(ind-1)+D_gsyn*dt;
     
%      D_Depr=(1-D(ind-1))/tau_D;
%      D(ind)=D(ind-1)+D_Depr*dt;
%      
%      D_gsec=-g_sec(ind-1)/tau_sec;
%      g_sec(ind)=g_sec(ind-1)+D_gsec*dt;
%     
%      D_F=(1-F(ind-1))/tau_fac;
%      F(ind)=F(ind-1)+D_F*dt;
%      
     D_V=((g_l)*(E_l-V(ind-1))+(g_syn(ind-1))*(E_syn-V(ind-1)))/C;
     V(ind)=V(ind-1)+D_V*dt;
     
    if spikes(ind)==1
        g_syn(ind)=g_syn(ind)+delta_g;
%        D(ind)=D(ind-1)-D(ind-1)*p_0;
%        F(ind)=F(ind-1)+f_fac*(F_max-F(ind-1));
%        g_sec(ind)=g_sec(ind)+g_max*p_0*D(ind);
    end
    
    if V(ind)>(-50e-3)
        V(ind)=-80e-3;
    end
end
subplot(2,1,2)
plot(t,g_syn)
ylabel('\fontsize{14}g_{syn}, Siemens')
xlabel('\fontsize{14}Time, seconds')
title('\fontsize{15}g_{syn}(t)')
grid on

subplot(2,1,1)
plot(t,V)
ylabel('\fontsize{14}V, Volts')
xlabel('\fontsize{14}Time, seconds')
title('\fontsize{15}V(t)')
grid on
