%% Written by Vardges on 04/02/2017.

clear all

%% Defining parameters
r_max=100;
theta_e=-5;
theta_i=0;
alpha_e=0.05;
alpha_i=1;

W_ee=2;
W_ei=2.5;
W_ie=-2.5;
W_ii=-2;

W_x=1.75;

tmax=3;
dt=1e-5;
t=0:dt:tmax;

tau_e=5e-3;
tau_i=5e-3;

Ie_base=25;
Ii_base=20;
Ie_stim=20;

%% Preallocating dynamic variables and defining currents
r_e=zeros(length(t),2);
r_i=zeros(length(t),2);

I_e=zeros(length(t),2);
I_i=zeros(length(t),2);

sig=0.5;
I_eapp=Ie_base*ones(length(t),2);
I_iapp=Ii_base*ones(length(t),2);
start=round(1/dt);
finish=round(1.1/dt);
I_eapp(start:finish,1)=(Ie_stim+Ii_base)*ones(length(start:finish),1);
start=round(2/dt);
finish=round(2.1/dt);
I_eapp(start:finish,2)=(Ie_stim+Ii_base)*ones(length(start:finish),1);
I_eapp=I_eapp.*(1+sig*randn(size(I_eapp)));

for k=2:length(t)
    %% First neuron
    I_e(k,1)=W_ee*r_e(k-1,1)+W_ie*r_i(k-1,1)+I_eapp(k-1,1);
    I_i(k,1)=W_ei*r_e(k-1,1)+W_ii*r_i(k-1,1)+I_iapp(k-1,1)+W_x*r_e(k-1,2);

    D_re=(1/tau_e)*(-r_e(k-1,1) + alpha_e*sign(I_e(k-1,1)-theta_e)*((I_e(k-1,1)-theta_e)^2));
    r_e(k,1)=r_e(k-1,1)+D_re*dt;
    
    D_ri=(1/tau_i)*(-r_i(k-1,1) + alpha_i*(I_i(k-1,1)-theta_i));
    r_i(k,1)=r_i(k-1,1)+D_ri*dt;
    
    r_e(k,1)=min(r_e(k,1), r_max);
    r_i(k,1)=min(r_i(k,1), r_max);
    
    r_e(k,1)=max(r_e(k,1), 0);
    r_i(k,1)=max(r_i(k,1), 0);
    
    
    %% Second neuron
    I_e(k,2)=W_ee*r_e(k-1,2)+W_ie*r_i(k-1,2)+I_eapp(k-1,2);
    I_i(k,2)=W_ei*r_e(k-1,2)+W_ii*r_i(k-1,2)+I_iapp(k-1,2)+W_x*r_e(k-1,1);
    
    D_re=(1/tau_e)*(-r_e(k-1,2) + alpha_e*sign(I_e(k-1,2)-theta_e)*((I_e(k-1,2)-theta_e)^2));
    r_e(k,2)=r_e(k-1,2)+D_re*dt;
    
    D_ri=(1/tau_i)*(-r_i(k-1,2) + alpha_i*(I_i(k-1,2)-theta_i));
    r_i(k,2)=r_i(k-1,2)+D_ri*dt;
    
    r_e(k,2)=min(r_e(k,2), r_max);
    r_i(k,2)=min(r_i(k,2), r_max);
    
    r_e(k,2)=max(r_e(k,2), 0);
    r_i(k,2)=max(r_i(k,2), 0);
    
        
end

close gcf
subplot(2,1,1)
plot(t, r_i(:,1), 'r')
hold on
plot(t, r_e(:,1), 'b')
set(gcf, 'WindowStyle', 'docked')
xlabel('\fontsize{14} Time, seconds')
ylabel('\fontsize{14} Firing rate, Hz')
title(['\fontsize{14} Firing rates, I_E^{(base)}=' num2str(Ie_base) ', I_I^{(base)}=' num2str(Ii_base)])
legend('Inhibitory units', 'Excitatory units')

subplot(2,1,2)
plot(t, r_i(:,2), 'r')
hold on
plot(t, r_e(:,2), 'b')
set(gcf, 'WindowStyle', 'docked')
xlabel('\fontsize{14} Time, seconds')
ylabel('\fontsize{14} Firing rate, Hz')
title(['\fontsize{14} Firing rates, I_E^{(base)}=' num2str(Ie_base) ', I_I^{(base)}=' num2str(Ii_base)])
legend('Inhibitory units', 'Excitatory units')
