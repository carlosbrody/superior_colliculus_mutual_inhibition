clear
clc

%%% try changing the actual numbers and put it in julia



%WEIGHTS
% wv=41; % vertical inhibitory weights
% wh=4;    % horizfffontal inhibitory weights

wv=9; % vertical inhibitory weights
wh=0.25;    % horizfffontal inhibitory weights


%INPUTS
bias=0.2135;  %constant pro bias input
delay=0.25; %delay input
light=3; %light input



opto=0.9; %opto reduction



noise_intracell_on=0; %apply noise to intracellular potentials
noisefr_on=0; %apply noise to firing rates

if(noise_intracell_on)
    noise_intracell=.1; %noise applied to "intracellular potential"
else
    noise_intracell=0; %noise applied to "intracellular potential"
end

if(noisefr_on)
    noisefr=.03; %noise applied to firing rates
else
    noisefr=0; %noise applied to firing rates
end




%%% weight matrix
w=[  0  -wv  -wh     0 ; ... %PRO_Contra_goL (1)
    -wv   0    0   -wh ; ... %ANTI_Ipsi_goR (2)
    -wh   0    0   -wv ; ... %PRO_Contra_goR (3)
    0  -wh  -wv     0 ] ;   %ANTI_Ipsi_goL (4)

%Parameters
dt=10; %integration timestep
totalTime=200+200+200; %simulation time (ms)
tau=17.6; %time constant
totalSteps=totalTime/dt; %total time steps
totalSteps

types=[0 1 1 1;... %pro, no opto
    0 opto opto opto;... %pro, opto all
    0 opto 1 1;... %pro, opto cue
    0 1 opto 1;... %pro, opto delay
    0 1 1 opto;... %pro, opto choice
    1 1 1 1;...    %anti, no opto
    1 opto opto opto;... %anti, opto all
    1 opto 1 1;... %anti, opto cue
    1 1 opto 1;... %anti, opto delay
    1 1 1 opto];   %anti, opto choice

labels={'pro, no opto';...
    'pro, opto all';...
    'pro, opto cue';...
    'pro, opto delay';...
    'pro, opto choice';...
    'anti, no opto';...
    'anti, opto all';...
    'anti, opto cue';...
    'anti, opto delay';...
    'anti, opto choice'};

%neurons, timepoints, trialtypes, repeats
yall=nan(4,totalSteps,10);


%loop over trial types
for iii=1:10
    
    
    anti=types(iii,1);
    optoC=types(iii,2);
    optoD=types(iii,3);
    optoCh=types(iii,4);
    
    
    %initialize matrix of voltages
    x=zeros(4,totalSteps);
    %initialize matrix of activations
    y=zeros(4,totalSteps);
    
    %set initial states
    x(:,1) = [-25 -25 -25 -25];
    y(:,1) = (tanh( (x(:,1)-5)/50 )+1)/2;
    
    %set inputs to target proper nodes
    bias_in=bias*[1 0 1 0];
    if anti
        delay_in=delay*[0 1 0 1];
    else
        delay_in=delay*[1 0 1 0];
    end
    light_in=light*[1 1 0 0];
    
    
    
    %------------------------ Cue period
    
    nume=totalSteps/3;
    for i=2:nume
        dx = -x(:,i-1)/tau + w*y(:,i-1) + delay_in' + bias_in';
        x(:,i) = x(:,i-1) + dx*dt;
        x(:,i) = x(:,i) + noise_intracell*randn(4,1);
        y(:,i) = (tanh( (x(:,i)-5)/50 )+1)/2;
        y(:,i) = optoC*y(:,i);
        y(:,i) = y(:,i) + noisefr*randn(4,1);
        
    end
    
    
    %------------------------ Delay period
    
    for i=nume+1:nume*2
        dx = -x(:,i-1)/tau + w*y(:,i-1) + delay_in' + bias_in';
        x(:,i) = x(:,i-1) + dx*dt;
        x(:,i) = x(:,i) + noise_intracell*randn(4,1);
        y(:,i) = (tanh( (x(:,i)-5)/50 )+1)/2;
        y(:,i) = optoD*y(:,i);
        y(:,i) = y(:,i) + noisefr*randn(4,1);
        
    end
    
    
    %------------------------ Choice period
    
    for i=nume*2+1:totalSteps
        dx = -x(:,i-1)/tau + w*y(:,i-1) + light_in' + bias_in' ;
        x(:,i) = x(:,i-1) + dx*dt;
        x(:,i) = x(:,i) + noise_intracell*randn(4,1);
        y(:,i) = (tanh( (x(:,i)-5)/50 )+1)/2;
        y(:,i) = optoCh*y(:,i);
        y(:,i) = y(:,i) + noisefr*randn(4,1);
        
    end
    
    
    yall(:,:,iii)=y;
end

colo='brcm';
figure
for i=1:10
    subplot(2,5,i)
    hold on
    for j=1:4
        a=squeeze(yall(j,:,i));
        plot(1:size(yall,2),a,colo(j))
    end
    title(labels{i})
end


