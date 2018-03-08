% Written by Vardges on the 29th of April, 2017.
% This script achieves 2D pattern recognition of noise-degraded stimuli.

%% Defining parameters
clear all

tau_r=10e-3;
r_max=50;
I_th=10;
delta=1;
r_T=25;
e_plus=0.1/289;
e_minus=e_plus/4;
Wmax=8/289;
Wmin=-Wmax;


tmax=1;
dt=10e-3;
t=0:dt:tmax;

%% Creating the four patterns 

% Making the tools
background=zeros(17,17);
line=ones(1,17);

% Creating letter Z
pattern_z=flip(diag(line), 2);
pattern_z(1, 1:end)=line;
pattern_z(end, 1:end)=line;

% Creating letter O
pattern_o=background;
pattern_o(1, 1:end)=line;
pattern_o(end, 1:end)=line;
pattern_o(1:end, 1)=line;
pattern_o(1:end, end)=line;

% Creating letter I
pattern_i=background;
pattern_i(1:end, 8)=line;

% Creating letter E
pattern_e=background;
pattern_e(1:end, 1)=line;
pattern_e(1, 1:8)=line(1:8);
pattern_e(end, 1:8)=line(1:8);
pattern_e(8, 1:8)=line(1:8);

% Stacking the patterns
stack=cat(3, pattern_z, pattern_o, pattern_i, pattern_e);


%% Learning the four patterns
for num=1:4
    
    sample=stack(1:end, 1:end, num);
    r_vec=train(sample);
   
    subplot(2,2,num)
    imagesc(reshape(r_vec, 17,17))    
    title('Firing rate matrix at the end of training', 'FontSize', 12)
    
end


function r_vec=train(sample)


tau_r=10e-3;
r_max=50;
I_th=10;
delta=1;
r_T=25;
e_plus=0.2*0.1/289;
e_minus=e_plus/4;
Wmax=8/289;
Wmin=-Wmax;


tmax=1;
dt=10e-3;
t=0:dt:tmax;
W=(-0.3/289)*ones(289);
plus=0;
minus=0;

for l=1:400
    
    noisy_input=sample+(rand(17)>0.5);
    noisy_input= (noisy_input==1);
   
    r_vec=noisy_input(:); 
    I_app=zeros(17, 17, length(t));
    I_app(:,:, 1:round(0.5/dt))=50*repmat(noisy_input, 1, 1, length(1:round(0.5/dt)));
    
   for k=2:length(t)
        
        % Calculating firing rates
        I_curr=I_app(:,:,k);
        I_curr=I_curr(:);
        I=(W*r_vec) + I_curr;
        D_r=(1/tau_r)*( -r_vec + (r_max./(1+exp(-(I-I_th)/delta))) );
        r_vec=r_vec+D_r*dt;
        
        % Calculating connectivity integrals
        [r_j, r_i]=meshgrid(r_vec);
        plus=plus+double(r_i>r_T).*double(r_j>r_T).*dt;
        minus=minus+double(r_T>r_i).*double(r_j>r_T).*dt;
            
    end
    
    delta_W= e_plus*plus-e_minus*minus;
    W=W+delta_W;
    W=W.*double(W<=Wmax)+Wmax*double(W>=Wmax);
    W=W.*double(W>=Wmin)+Wmin*double(W<=Wmin);    
    
end

end
