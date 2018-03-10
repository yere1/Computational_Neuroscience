% Written by Vardges on the 25th of April, 2017
% This script executes a simple neural net with  five inputs and
% two outputs. It takes the five inputs to deciede, if the 
% weather will be "rainy" or "sunny". The script outputs plots, which
% show the evolution of synaptic weights through time. The script
% offers choice of 5 learning rules.
%%
clear all

rule=input('Which rule?  ');

P=[0.05, 0.25, 0.5, 0.75, 0.95; 0.95, 0.75, 0.5, 0.25, 0.05]';
tau=20e-3;
r_max=100;
I_th=50; 
I_sigma=5;
tmax=0.5;
dt=1e-3;
t=0:dt:tmax;
revisor=0;

start=round(0.1/dt);
I_app=zeros(5, length(t));

W_inputs=0.2*ones(2,5);
W_decision=[0.5 -0.5; -0.5 0.5];
W=[W_inputs W_decision];

I_input=zeros(2,length(t));
r_input=zeros(5,length(t));
r_decision=zeros(2, length(t));

%%
for trial=1:800
    
      cue=1+double(rand>0.5);
      is_there_current=double(rand(5,1)<P(:, cue));
      I_app(1:5, start:end)=50*repmat(is_there_current,1, length(t(start:end)));
          
for k=2:length(t)
    
    % Euler method calculations for the input units
    
    D_rinput=(1/tau)*(-r_input(:, k-1) + r_max./(1+exp((I_th-I_app(:, k))./I_sigma)));
    r_input(:, k)= r_input(:, k-1)+D_rinput*dt;
    
    % Euler method calculations for decision-making units
    rate_mat=repmat(([r_input(:,k); r_decision(:,k-1)]'),2,1);
    I_rnd=(1/sqrt(dt))*randn(2,1);
    I_input(:,k)=sum(W.*rate_mat,2) + I_rnd;
    D_rdecision=(1/tau)*(-r_decision(:, k-1) + r_max./(1+exp((I_th-I_input(:, k))./I_sigma)));
    r_decision(:, k)= r_decision(:, k-1)+D_rdecision*dt;
    
        if exist('active')==0
            
            if r_decision(1,k)>=40
                active=1;
                inactive=2;
            end
            
            if r_decision(2,k)>=40
                active=2;
                inactive=1;
            end
                
        end
    
end

    if exist('active')==0
       
        if r_decision(1,end)>r_decision(2,end)
            active=1;
            inactive=2;
        else
            active=2;
            inactive=1;   
            
        end
    end

    switch rule
        
        case 1
        
            input=[r_input(:, end)>30]';
            E= double(active==cue)-0.5;
            d_W(active,:)= 0.04.*E.*input;
            d_W(inactive,:)= -0.04.*E.*input;
            W(:,1:5)=W(:,1:5)+d_W;
            Z(:,:, trial)=W;
            
            revisor=revisor+double(cue==active);
                   
        case 2
            
            R(trial)=double(cue==active);
            if trial<11
                Ravg=0.5;
            else
                Ravg=mean(R(trial-10 : trial));
            end
            
            input=[r_input(:, end)>30]';
            E=R(trial)-Ravg;
            d_W(active,:)= 0.04.*E.*input;
            d_W(inactive,:)= -0.04.*E.*input;
            W(:,1:5)=W(:,1:5)+d_W;
            Z(:,:, trial)=W;
            
            revisor=revisor+double(cue==active);
            
        case 3
            
            input=[r_input(:, end)>30]';
            E= double(active==cue)-0.5;
                        
            if E>0
                d_W(active,:)= 0.04.*E.*input.*(0.4-W(active,1:5));
                d_W(inactive,:)= -0.04.*E.*input.*W(inactive,1:5);
            end
            
            if E<0
                d_W(active,:)= 0.04.*E.*input.*W(active,1:5);
                d_W(inactive,:)= -0.04.*E.*input.*(0.4-W(inactive,1:5));
            end
            
            W(:,1:5)=W(:,1:5)+d_W;
            Z(:,:, trial)=W;
            
            revisor=revisor+double(cue==active);
            
        case 4
            
            R(trial)=double(cue==active);
            if trial<11
                Ravg=0.5;
            else
                Ravg=mean(R(trial-10 : trial));
            end
            
            input=[r_input(:, end)>30]';
            E=R(trial)-Ravg;
            d_W(active,:)= 0.04.*E.*input;
            W(:,1:5)=W(:,1:5)+d_W;
            Z(:,:, trial)=W;
            
            revisor=revisor+double(cue==active);
            
        case 5
            
            R(trial)=double(cue==active);
            if trial<11
                Ravg=0.5;
            else
                Ravg=mean(R(trial-10 : trial));
            end
            
            input=[r_input(:, end)>30]';
            E=R(trial)-Ravg;
            d_W(active,:)= 0.04.*E.*input+0.04.*E.*(input-1);
            W(:,1:5)=W(:,1:5)+d_W;
            Z(:,:, trial)=W;
            
            revisor=revisor+double(cue==active);
                       
    end
    
    clear active inactive
    
end
 
figure
subplot(2,1,1)
plot(squeeze(Z(1,1,:)))
hold on
plot(squeeze(Z(1,2,:)))
plot(squeeze(Z(1,3,:)))
plot(squeeze(Z(1,4,:)))
plot(squeeze(Z(1,5,:)))
hold off
title(sprintf('Synapse evolution, Rule #%d\nConnections to the 1st unit\nCorrect responses: %d', rule, revisor), 'FontSize', 14);
ylabel('Synaptic strengths', 'FontSize', 14)
xlabel('Trial number', 'FontSize', 14)
legend('A', 'B', 'C', 'D', 'E', 'Location', 'southwest', 'Orientation', 'horizontal')
grid on

subplot(2,1,2)
plot(squeeze(Z(2,1,:)))
hold on
plot(squeeze(Z(2,2,:)))
plot(squeeze(Z(2,3,:)))
plot(squeeze(Z(2,4,:)))
plot(squeeze(Z(2,5,:)))
hold off
title(sprintf('Synapse evolution, Rule #%d\nConnections to the 2nd unit\nCorrect responses: %d', rule, revisor), 'FontSize', 14);
ylabel('Synaptic strengths', 'FontSize', 14)
xlabel('Trial number', 'FontSize', 14)
legend('A', 'B', 'C', 'D', 'E', 'Location', 'southwest', 'Orientation', 'horizontal')
grid on

figure
x=log(P(:,1)./P(:,2));
y=diff(W);
plot(x, y(1:5), 'LineWidth', 2)
title(sprintf('Rule #%d', rule), 'FontSize', 14)
ylabel('Final weight diff.', 'FontSize', 14)
xlabel('Log prob. ratio', 'FontSize', 14)
grid on
