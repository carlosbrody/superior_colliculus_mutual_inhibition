function run_solution(farm_id, examples, p )
dl          = 1.1;          %x/y limits for each PC space
delay       = 0.1;         % delay in seconds between plotting each time step
numtrials   = 50;           % number of trials in data set for each solution
nump        = numtrials/2;
numsteps    = 61;           % total number of timesteps
numrdsteps  = 41;           % number of timesteps in rule and delay period
plot_hits   = p.plot_hits;  
plot_traj   = p.plot_traj;

% get current farm, split into trial epochs
x  = squeeze(examples(farm_id, 1,:,:,:,:));
rd = squeeze(examples(farm_id, 1,:,:,1:numrdsteps,:));
t  = squeeze(examples(farm_id, 1,:,:,(numrdsteps+1):end,:));

% split into pro/anti
p = squeeze(x(1,:,:,:));
a = squeeze(x(2,:,:,:));

% determine hits and misses
phits = [squeeze(p(1,end,:)) > squeeze(p(4,end,:))];
phits((nump+1):end) = ~phits((nump+1):end);
ahits = [squeeze(a(4,end,:)) > squeeze(a(1,end,:))];
ahits((nump+1):end) = ~ahits((nump+1):end);

% determine which trials to plot
if plot_hits;
    pplot = phits;
    aplot = ahits;
else;
    pplot = ~phits;
    aplot = ~ahits;  
end


% do pca on rule and delay period
temp = reshape(rd,2,4,numrdsteps*numtrials);
d = squeeze([squeeze(temp(1,:,:)) squeeze(temp(2,:,:))]);
% get mean here, then need to reshape to apply to p/a
rd_mean = mean(d,2);
c = cov(d');
[dvecs, dvals] = eig(c);

% do pca on target period
temp = reshape(t,2,4,(numsteps-numrdsteps)*numtrials);
d = squeeze([squeeze(temp(1,:,:)) squeeze(temp(2,:,:))]);
% get mean here, then need to reshape to apply to p/a
t_mean = mean(d,2);
c = cov(d');
[tvecs, tvals] = eig(c);

% subtract off mean for projecting data
% have to subtract different means for first and second half of each trial
p(:,1:numrdsteps,:) = p(:,1:numrdsteps,:) - repmat(rd_mean, 1, length(1:numrdsteps), size(p,3));
p(:,(numrdsteps+1):end,:) = p(:,(numrdsteps+1):end,:) - repmat(t_mean, 1, length((numrdsteps+1):size(p,2)), size(p,3));
a(:,1:numrdsteps,:) = a(:,1:numrdsteps,:) - repmat(rd_mean, 1, length(1:numrdsteps), size(a,3));
a(:,(numrdsteps+1):end,:) = a(:,(numrdsteps+1):end,:) - repmat(t_mean, 1, length((numrdsteps+1):size(a,2)), size(a,3));

% plot solution
close all; figure; 
subplot(1,2,1); hold on;        set(gca, 'fontsize',20);
ylabel('Delay Period PC 2');    xlabel('Delay Period PC 1')
ylim([-dl,dl]); xlim([-dl,dl]); axis square
title('Pro')
subplot(1,2,2); hold on;        set(gca, 'fontsize',20);
ylabel('Target Period PC 2');   xlabel('Target Period PC 1')
ylim([-dl,dl]); xlim([-dl,dl]); axis square
title('timestep -')

% plot PRO TRIALS
if plot_traj
for i=1:numsteps
    subplot(1,2,1);
    title('Pro')
    plot(dvecs(:,4)'*squeeze(p(:,i,pplot)), dvecs(:,3)'*squeeze(p(:,i,pplot)),'ko')

    subplot(1,2,2);
    title(['timestep - ' num2str(i)])
    plot(tvecs(:,4)'*squeeze(p(:,i,pplot)), tvecs(:,3)'*squeeze(p(:,i,pplot)),'ko')

    if i >= (numrdsteps+1)
    subplot(1,2,1);
    plot(dvecs(:,4)'*squeeze(p(:,i,logical([pplot(1:nump); zeros(nump,1)]))), dvecs(:,3)'*squeeze(p(:,i,logical([pplot(1:nump); zeros(nump,1)]))),'ro')
    plot(dvecs(:,4)'*squeeze(p(:,i,logical([zeros(nump,1); pplot((nump+1):end)]))), dvecs(:,3)'*squeeze(p(:,i,logical([zeros(nump,1);pplot((nump+1):end)]))),'mo')

    subplot(1,2,2);
    plot(tvecs(:,4)'*squeeze(p(:,i,logical([pplot(1:nump); zeros(nump,1)]))), tvecs(:,3)'*squeeze(p(:,i,logical([pplot(1:nump); zeros(nump,1)]) )),'ro')
    plot(tvecs(:,4)'*squeeze(p(:,i,logical([zeros(nump,1); pplot((nump+1):end)]) )), tvecs(:,3)'*squeeze(p(:,i,logical([zeros(nump,1); pplot((nump+1):end)]))),'mo')
    end

    pause(delay)

end
end
% plot PRO ENDING POINTS
subplot(1,2,1);
plot(dvecs(:,4)'*squeeze(p(:,end,logical([pplot(1:nump); zeros(nump,1)]))), dvecs(:,3)'*squeeze(p(:,end,logical([pplot(1:nump); zeros(nump,1)]))),'ro','markerfacecolor','r')
plot(dvecs(:,4)'*squeeze(p(:,end,logical([zeros(nump,1); pplot((nump+1):end)]))), dvecs(:,3)'*squeeze(p(:,end,logical([zeros(nump,1); pplot((nump+1):end)]))),'mo','markerfacecolor','m')
                                   
subplot(1,2,2);
title('Final')
plot(tvecs(:,4)'*squeeze(p(:,end,logical([pplot(1:nump); zeros(nump,1)]))), tvecs(:,3)'*squeeze(p(:,end,logical([pplot(1:nump); zeros(nump,1)]))),'ro','markerfacecolor','r')
plot(tvecs(:,4)'*squeeze(p(:,end,logical([zeros(nump,1); pplot((nump+1):end)]))), tvecs(:,3)'*squeeze(p(:,end,logical([zeros(nump,1); pplot((nump+1):end)]))),'mo','markerfacecolor','m')

pause(2)
% plot ANTI TRIALS
if plot_traj
for i=1:numsteps
    subplot(1,2,1);
    title('Anti')
    plot(dvecs(:,4)'*squeeze(a(:,i,aplot)), dvecs(:,3)'*squeeze(a(:,i,aplot)),'bo')

    subplot(1,2,2);
    title(['timestep - ' num2str(i)])
    plot(tvecs(:,4)'*squeeze(a(:,i,aplot)), tvecs(:,3)'*squeeze(a(:,i,aplot)),'bo')
    if i >= (numrdsteps+1)
    subplot(1,2,1);
    plot(dvecs(:,4)'*squeeze(a(:,i,logical([aplot(1:nump); zeros(nump,1)]))), dvecs(:,3)'*squeeze(a(:,i,logical([aplot(1:nump); zeros(nump,1)]))),'go')
    plot(dvecs(:,4)'*squeeze(a(:,i,logical([zeros(nump,1); aplot((nump+1):end)]))), dvecs(:,3)'*squeeze(a(:,i,logical([zeros(nump,1); aplot((nump+1):end)]))),'co')
    subplot(1,2,2);
    plot(tvecs(:,4)'*squeeze(a(:,i,logical([aplot(1:nump); zeros(nump,1)]))), tvecs(:,3)'*squeeze(a(:,i,logical([aplot(1:nump); zeros(nump,1)]))),'go')
    plot(tvecs(:,4)'*squeeze(a(:,i,logical([zeros(nump,1); aplot((nump+1):end)]))), tvecs(:,3)'*squeeze(a(:,i,logical([zeros(nump,1); aplot((nump+1):end)]))),'co')
    end
    pause(delay)

end
end
% PLOT ANTI ENDING POINTS
subplot(1,2,1);
title('')
plot(dvecs(:,4)'*squeeze(a(:,end,logical([aplot(1:nump); zeros(nump,1)]))), dvecs(:,3)'*squeeze(a(:,end,logical([aplot(1:nump); zeros(nump,1)]))),'go','markerfacecolor','g')
plot(dvecs(:,4)'*squeeze(a(:,end,logical([zeros(nump,1); aplot((nump+1):end)]) )), dvecs(:,3)'*squeeze(a(:,end,logical([zeros(nump,1); aplot((nump+1):end)]))),'co','markerfacecolor','c')
                                   
subplot(1,2,2);
title('final')
plot(tvecs(:,4)'*squeeze(a(:,end,logical([aplot(1:nump); zeros(nump,1)]))), tvecs(:,3)'*squeeze(a(:,end,logical([aplot(1:nump); zeros(nump,1)]))),'go','markerfacecolor','g')
plot(tvecs(:,4)'*squeeze(a(:,end,logical([zeros(nump,1); aplot((nump+1):end)]))), tvecs(:,3)'*squeeze(a(:,end,logical([zeros(nump,1); aplot((nump+1):end)]))),'co','markerfacecolor','c')


% Plot the decision plane
% Not sure how to do this correctly
% method 1 , take vector orthogonal to decision hyperplane, project into space, then find perpendicular line
%ortho = [-1;0;0;1];
%ortho = ortho./norm(ortho);
%portho = [tvecs(:,4)'*ortho; tvecs(:,3)'*ortho];
%pdec = [-portho(2); portho(1)];
%pdec = pdec./norm(pdec);
%subplot(1,2,2);
%plot([-pdec(1)*2 pdec(1)*2], [-pdec(2)*2 pdec(2)*2], 'b--')

% method 2, projection decision line into space. Not sure why these are different???
dec = [1;0;0;1];
dec = dec./norm(dec);
ppdec= [tvecs(:,4)'*dec; tvecs(:,3)'*dec];
ppdec = ppdec./norm(ppdec);
plot([-ppdec(1)*2 ppdec(1)*2], [-ppdec(2)*2 ppdec(2)*2], 'r--')



% Data has been mean centered, which means decision line should be offset. Need to figure this out. 
% 
%pdisp = [tvecs(:,4)'*t_mean; tvecs(:,3)'*t_mean];
%new_vec = ppdec + pdisp;
%plot([-new_vec(1)*2 new_vec(1)*2], [-new_vec(2)*2 new_vec(2)*2], 'r--')
%new_vec2 = pdec + pdisp;
%plot([-new_vec2(1)*2 new_vec2(1)*2], [-new_vec2(2)*2 new_vec2(2)*2], 'c--')





