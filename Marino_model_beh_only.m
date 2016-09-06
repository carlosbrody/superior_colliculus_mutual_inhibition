clear
clc


%WEIGHTS
wv=36; % vertical inhibitory weights
wh=1;    % horizontal inhibitory weights


%INPUTS
bias=0.854;  %constant pro bias input
delay=1; %delay input
light=12; %light input



opto=0.8; %opto reduction



noise_intracell_on=0; %apply noise to intracellular potentials
noisefr_on=0; %apply noise to firing rates

if(noise_intracell_on)
    noise_intracell=.1; %noise applied to "intracellular potential"
else
    noise_intracell=0; %noise applied to "intracellular potential"
end

if(noisefr_on)
    noisefr=.01; %noise applied to firing rates
else
    noisefr=0; %noise applied to firing rates
end




%%% weight matrix
w=[  0  -wv  -wh     0 ; ... %PRO_Contra_goL (1)
    -wv   0    0   -wh ; ... %ANTI_Ipsi_goR (2)
    -wh   0    0   -wv ; ... %PRO_Contra_goR (3)
     0  -wh  -wv     0 ] ;   %ANTI_Ipsi_goL (4)

%Parameters
dt=.05; %integration timestep
totalTime=200+200+200; %simulation time (ms)
tau=4.4; %time constant
totalSteps=totalTime/dt; %total time steps


types=[0 1 1 1;... %pro, no opto
    0 opto 1 1;... %pro, opto cue
    0 1 opto 1;... %pro, opto delay
    0 1 1 opto;... %pro, opto choice
    1 1 1 1;...    %anti, no opto
    1 opto 1 1;... %anti, opto cue
    1 1 opto 1;... %anti, opto delay
    1 1 1 opto];   %anti, opto choice

labels={'pro, no opto';... 
    'pro, opto cue';... 
    'pro, opto delay';...
    'pro, opto choice';...
    'anti, no opto';...
    'anti, opto cue';...
    'anti, opto delay';...
    'anti, opto choice'};


yall=nan(4,totalSteps,8);
for iii=1:8
    
    
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
    
    for i=2:4000
        dx = -x(:,i-1)/tau + w*y(:,i-1) + delay_in' + bias_in';
        x(:,i) = x(:,i-1) + dx*dt;
        x(:,i) = x(:,i) + noise_intracell*randn(4,1);
        y(:,i) = (tanh( (x(:,i)-5)/50 )+1)/2;
        y(:,i) = optoC*y(:,i);
        y(:,i) = y(:,i) + noisefr*randn(4,1);
        
    end
        
    
    %------------------------ Delay period
    
    for i=4001:8000
        dx = -x(:,i-1)/tau + w*y(:,i-1) + delay_in' + bias_in';
        x(:,i) = x(:,i-1) + dx*dt;
        x(:,i) = x(:,i) + noise_intracell*randn(4,1);
        y(:,i) = (tanh( (x(:,i)-5)/50 )+1)/2;
        y(:,i) = optoD*y(:,i);
        y(:,i) = y(:,i) + noisefr*randn(4,1);
        
    end
    
    
    %------------------------ Choice period
    
    for i=8001:12000
        dx = -x(:,i-1)/tau + w*y(:,i-1) + light_in' + bias_in' ;
        x(:,i) = x(:,i-1) + dx*dt;
        x(:,i) = x(:,i) + noise_intracell*randn(4,1);
        y(:,i) = (tanh( (x(:,i)-5)/50 )+1)/2;
        y(:,i) = optoCh*y(:,i);
        y(:,i) = y(:,i) + noisefr*randn(4,1);
        
    end
    
    
    yall(:,:,iii)=y;
end



%%% PLOTS! %%%

ylb=-0.01;
yub=41;
yscale = 45;
linew = 2;
T=dt:dt:totalTime;

figure
for iii=1:8
    
    subplot(2,4,iii);
    
    
    hold on
    plot([-.7 -.7],[0 41],'k','linewidth',.2)
    hold on
    plot([0 0],[0 41],'k','linewidth',.2)
    
    hold on
    plot((T-400)/(200/.7),yall(1,:,iii)*yscale,'Color',[0, 0, 1],'linewidth',linew)
    plot((T-400)/(200/.7),yall(3,:,iii)*yscale,'Color',[1, 0.7, 0.5],'linewidth',linew)
    plot((T-400)/(200/.7),yall(2,:,iii)*yscale,'Color',[1, 0, 0],'linewidth',linew)
    plot((T-400)/(200/.7),yall(4,:,iii)*yscale,'Color',[0, 1, 1],'linewidth',linew)
    ylim([ylb yub])
    xlim([-1.4 0.7])
    
    ylabel('Firing rate (Hz)')
    xlabel('Time (s)')
    set(gca,'fontsize',14)
    set(gca,'linewidth',1)
    title(labels{iii});
    
    set(gca,'TickDir','out')
    set(gca,'FontSize',23)
    box off
    
end







