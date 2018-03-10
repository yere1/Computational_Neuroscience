% Written by Vardges on February 14th, 2017
% This script is a simulation of a Hodgkin-Huxley neuron
% with different input currents.

clear all
disp('Input "a" though "e". Each letter corresponds to a different external current to the cell.');
section=input('Enter a letter. \n', 's');

if section=='c'
    param=input('Enter the delay parameter in seconds. \n');
end

global t dt;
dt=1e-7;
tmax=0.3;
t=0:dt:tmax;
I_0=zeros(1,length(t)); 
m(1)=0; h(1)=0; n(1)=0;

switch section 

    case 'a'
        
        V=HH(I_0,m(1),h(1),n(1));
        
    case 'b'
        
        step_start=(50e-3)/dt;  
        step_finish=(150e-3)/dt;
        I_0(step_start:step_finish)= floor(t(1:(step_finish-step_start+1))/0.01)*(0.018e-9);
        I_0(step_finish:end)=(0.18e-9)*ones(1, length(I_0)-step_finish+1);
        V=HH(I_0,m(1),h(1),n(1)); 
        
      case 'c'
          
        cycle=[zeros(1,round(param/dt)) ones(1, round((5e-3)/dt))];
        full_cycle=repmat(cycle, 1, 10);
        I_0(1:length(full_cycle))=(0.18e-9)*full_cycle;        
        V=HH(I_0,m(1),h(1),n(1));
        
     case 'd'
         
        m(1)=0.05; h(1)=0.5; n(1)=0.35; 
        cycle=[zeros(1,round((20e-3)/dt)) ones(1, round((5e-3)/dt))];
        full_cycle=repmat(cycle,1,10);        
        I_0=(0.65e-9)*ones(1,length(t));
        I_0(1:length(full_cycle))=I_0(1:length(full_cycle))-(0.65e-9)*full_cycle;
        V=HH(I_0,m(1),h(1),n(1));
         
         
     case 'e'
        
        pulse_start=round((50e-3)/dt);
        pulse_end=pulse_start+round((5e-3)/dt);
        m(1)=0.05; h(1)=0.5; n(1)=0.35;
        I_0=(0.70e-9)*ones(1,length(t));
        I_0(pulse_start:pulse_end)=(1e-9)*ones(1,(5e-3)/dt+1);        
        V=HH(I_0,m(1),h(1),n(1));
                 
     case 'f'
         
        pulse_start=round((50e-3)/dt);
        pulse_end=pulse_start+round((5e-3)/dt);
        I_0=(0.70e-9)*ones(1,length(t));
        I_0(pulse_start:pulse_end)=(1e-9)*ones(1,50001);        
        V=HH(I_0,m(1),h(1),n(1));
        
    otherwise
        
        error('Section does not exist')        
          
end
        
subplot(2,1,1)
plot(t,V, 'r')
xlabel('\fontsize{12} Time, seconds')
ylabel('\fontsize{12} Voltage, Volts')
title(['\fontsize{17} Voltage vs. time, Section ', section])
xlim([0 tmax])
grid on

subplot(2,1,2)
plot(t,I_0, 'b')
xlabel('\fontsize{12} Time, seconds')
ylabel('\fontsize{12} Current, Amperes')
title(['\fontsize{17} Current vs. time, Section ', section])
ylim([0 1.1e-9])
xlim([0 tmax])
grid on


function V=HH(I_0, M, H, N)

%---Defining parameters---
G_l=30e-9; %Siemens
G_na=12e-6; %Siemens
G_k=4e-6; %Siemens
E_na=50e-3; %Volts
E_k=-75e-3; %Volts
E_l=-65e-3; %Volts
C_m=100e-12; %Farad
global t dt;

%---Defining dynamical variables---
m=zeros(size(t)); m(1)=M;
h=zeros(size(t)); h(1)=H;
n=zeros(size(t)); n(1)=N;
V=zeros(size(t)); V(1)=E_l; 


%---Equation solver---
for k=2:length(t)
    
    %---Section taken from Latte---
    Vm = V(k-1);          % membane potential for calculations
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta.
    
    % First, sodium activation rate
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else                    % potassium activation rate
        alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    
    beta_n = 125*exp((-Vm-0.075)/0.08);     % potassium deactivation rate
    
    
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    %---End of section taken from Latte---
    
%---Dynamical equations---

% Sodium activation gating 
D_m=alpha_m*(1-m(k-1))-beta_m*m(k-1);
m(k)=m(k-1)+D_m*dt;

% Sodium inactivation gating
D_h=alpha_h*(1-h(k-1))-beta_h*h(k-1);
h(k)=h(k-1)+D_h*dt;

% Potassium gating
D_n=alpha_n*(1-n(k-1))-beta_n*n(k-1);
n(k)=n(k-1)+D_n*dt;

% Membrane voltage
D_V=(1/C_m)*(G_l*(E_l-Vm) + G_na*m(k)^3*h(k)*(E_na-Vm) + G_k*n(k)^4*(E_k-Vm) + I_0(k));
V(k)=Vm+D_V*dt;

end

end
