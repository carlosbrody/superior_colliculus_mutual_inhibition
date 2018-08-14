% Script for plotting Extended figure 7, which shows the average response for each cluster, as well as an example farm
close all
clear all


% loads a 6 d matrix, 350 farms x 3 opto x pro/anti x 4 nodes, x 76 time steps, x 10 examples
load MiniC30_C30_examples.mat

% set up time vector
ntim    =76;
tvec    = 1:ntim;
dt      = 0.024;
tvec    = tvec.*dt;
target_on   = 1.2+dt;
target_off  = target_on + 0.3+dt;

% need to filter examples by cost
costs   = [results.cost{:}];
cdex    = costs < -.00025;
good_examples = examples(cdex, :,:,:,:,:);

% make list of clusters
ids = sort(unique(cluster_ids));
ids(isnan(ids))= [];

% cluster colors
colors{1}='\color{blue}';
colors{2}='\color{green}';
colors{3}='\color{red}';
colors{4}='\color{cyan}';
colors{5}='\color{magenta}';
colors{6}='\color{yellow}';

% node colors
colo2{1}=[0 0 1];
colo2{2}=[1 0 0];
colo2{3}=[1 .5 .5];
colo2{4}=[0 1 1];

for i=1:length(ids)
    figure;
%% plot cluster average dynamics
    subplot(4,4,1); hold on;    set(gca, 'fontsize',20)
    cluster= squeeze(mean(mean(good_examples(cluster_ids == ids(i),1,1,:,:,:),6)));
    plot([target_on target_on], [0 1], 'k--')
    plot([target_off target_off], [0 1], 'k--')
    plot(tvec,cluster(1,:),'Color', colo2{1}, 'linewidth',2)
    plot(tvec,cluster(2,:),'Color', colo2{2}, 'linewidth',2)
    plot(tvec,cluster(3,:),'Color', colo2{3}, 'linewidth',2)
    plot(tvec,cluster(4,:),'Color', colo2{4}, 'linewidth',2)
    ylabel({'Cluster Average';'Activation (au)'},'fontsize',20)
    str=['Cluster ' num2str(i) '; Pro Trial'];
    title([colors{i} str],'fontsize',20)    
    cluster_dex(1,:) = cluster(1,:) - cluster(4,:);
    rule_dex(1,:) = .5*(cluster(1,:) + cluster(4,:)) - .5*(cluster(2,:) + cluster(3,:));
    xlim([0 tvec(end)]); ylim([0 1])

    subplot(4,4,2); hold on;    set(gca, 'fontsize',20)
    cluster = squeeze(mean(mean(good_examples(cluster_ids == ids(i),1,2,:,:,:),6)));
    plot([target_on target_on], [0 1], 'k--')
    plot([target_off target_off], [0 1], 'k--')
    plot(tvec,cluster(1,:),'Color', colo2{1}, 'linewidth',2)
    plot(tvec,cluster(2,:),'Color', colo2{2}, 'linewidth',2)
    plot(tvec,cluster(3,:),'Color', colo2{3}, 'linewidth',2)
    plot(tvec,cluster(4,:),'Color', colo2{4}, 'linewidth',2)
    str=['Anti Trial'];
    title([colors{i} str],'fontsize',20)        
    cluster_dex(2,:) = cluster(1,:) - cluster(4,:);
    rule_dex(2,:) = .5*(cluster(1,:) + cluster(4,:)) - .5*(cluster(2,:) + cluster(3,:));
    xlim([0 tvec(end)]); ylim([0 1])

    subplot(4,4,3); hold on;    set(gca, 'fontsize',20)
    cluster = squeeze(mean(mean(good_examples(cluster_ids == ids(i),2,2,:,:,:),6)));
    plot([target_on target_on], [0 1], 'k--')
    plot([target_off target_off], [0 1], 'k--')
    plot(tvec,cluster(1,:),'Color', colo2{1}, 'linewidth',2)
    plot(tvec,cluster(2,:),'Color', colo2{2}, 'linewidth',2)
    plot(tvec,cluster(3,:),'Color', colo2{3}, 'linewidth',2)
    plot(tvec,cluster(4,:),'Color', colo2{4}, 'linewidth',2)
    str=['Anti - delay opto'];
    title([colors{i} str],'fontsize',20)        
    cluster_dex(3,:) = cluster(1,:) - cluster(4,:);
    rule_dex(3,:) = .5*(cluster(1,:) + cluster(4,:)) - .5*(cluster(2,:) + cluster(3,:));
    xlim([0 tvec(end)]); ylim([0 1])

    subplot(4,4,4); hold on;    set(gca, 'fontsize',20)
    cluster = squeeze(mean(mean(good_examples(cluster_ids == ids(i),3,2,:,:,:),6)));
    plot([target_on target_on], [0 1], 'k--')
    plot([target_off target_off], [0 1], 'k--')
    plot(tvec,cluster(1,:),'Color', colo2{1}, 'linewidth',2)
    plot(tvec,cluster(2,:),'Color', colo2{2}, 'linewidth',2)
    plot(tvec,cluster(3,:),'Color', colo2{3}, 'linewidth',2)
    plot(tvec,cluster(4,:),'Color', colo2{4}, 'linewidth',2)
    str=['Anti - choice opto'];
    title([colors{i} str],'fontsize',20)        
    cluster_dex(4,:) = cluster(1,:) - cluster(4,:);
    rule_dex(4,:) = .5*(cluster(1,:) + cluster(4,:)) - .5*(cluster(2,:) + cluster(3,:));
    xlim([0 tvec(end)]); ylim([0 1])

% set up axis
    subplot(4,4,9); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    plot([tvec(1) tvec(end)], [0 0], 'k--')
    ylabel({'Choice Index';'(Pro-Anti)'}, 'fontsize',20)
    xlim([0 tvec(end)])
    ylim([-1 1])

    subplot(4,4,10); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    plot([tvec(1) tvec(end)], [0 0], 'k--')
    xlim([0 tvec(end)])
    ylim([-1 1])

    subplot(4,4,11); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    plot([tvec(1) tvec(end)], [0 0], 'k--')
    xlim([0 tvec(end)])
    ylim([-1 1])

    subplot(4,4,12); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    plot([tvec(1) tvec(end)], [0 0], 'k--')
    xlim([0 tvec(end)])
    ylim([-1 1])

    subplot(4,4,13); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    plot([tvec(1) tvec(end)], [0 0], 'k--')
    ylabel({'Rule Index';'(Pro-Anti)'}, 'fontsize',20)
    xlim([0 tvec(end)])
    ylim([-1 1])
    xlabel('time(s)')

    subplot(4,4,14); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    plot([tvec(1) tvec(end)], [0 0], 'k--')
    xlim([0 tvec(end)])
    ylim([-1 1])
    xlabel('time(s)')

    subplot(4,4,15); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    plot([tvec(1) tvec(end)], [0 0], 'k--')
    xlim([0 tvec(end)])
    ylim([-1 1])
    xlabel('time(s)')

    subplot(4,4,16); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    plot([tvec(1) tvec(end)], [0 0], 'k--')
    xlim([0 tvec(end)])
    ylim([-1 1])
    xlabel('time(s)')

%%%%% plot example solution
    subplot(4,4,5); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    dex = find(cluster_ids == ids(i),1);
    trajs = squeeze(good_examples(dex, 1, 1, :, :, :));
    ylabel({'Example Model';'Activation (au)'}, 'fontsize',20)
    xlim([0 tvec(end)]); ylim([0 1])
    for j=1:10
        subplot(4,4,5); hold on; set(gca, 'fontsize',20);
        plot(tvec,trajs(1,:,j),'Color', colo2{1}, 'linewidth',1)
        plot(tvec,trajs(2,:,j),'Color', colo2{2}, 'linewidth',1)
        plot(tvec,trajs(3,:,j),'Color', colo2{3}, 'linewidth',1)
        plot(tvec,trajs(4,:,j),'Color', colo2{4}, 'linewidth',1)

        subplot(4,4,9); hold on; set(gca, 'fontsize',20);    
        plot(tvec, trajs(1,:,j) - trajs(4,:,j), 'color',[.7 .7 .7], 'linewidth',1)

        subplot(4,4,13); hold on; set(gca, 'fontsize',20);    
        plot(tvec, .5*(trajs(1,:,j)+trajs(4,:,j))-.5*(trajs(2,:,j)+trajs(3,:,j)), 'color',[.7 .7 .7], 'linewidth',1)

    end


    subplot(4,4,6); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    dex = find(cluster_ids == ids(i),1);
    trajs = squeeze(good_examples(dex, 1, 2, :, :, :));
    xlim([0 tvec(end)]); ylim([0 1])
    for j=1:10
        subplot(4,4,6); hold on; set(gca, 'fontsize',20);    
        plot(tvec,trajs(1,:,j),'Color', colo2{1}, 'linewidth',1)
        plot(tvec,trajs(2,:,j),'Color', colo2{2}, 'linewidth',1)
        plot(tvec,trajs(3,:,j),'Color', colo2{3}, 'linewidth',1)
        plot(tvec,trajs(4,:,j),'Color', colo2{4}, 'linewidth',1)

        subplot(4,4,10); hold on; set(gca, 'fontsize',20);    
        plot(tvec, trajs(1,:,j) - trajs(4,:,j), 'color',[.7 .7 .7], 'linewidth',1)

        subplot(4,4,14); hold on; set(gca, 'fontsize',20);    
        plot(tvec, .5*(trajs(1,:,j)+trajs(4,:,j))-.5*(trajs(2,:,j)+trajs(3,:,j)), 'color',[.7 .7 .7], 'linewidth',1)
    end

    subplot(4,4,7); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    dex = find(cluster_ids == ids(i),1);
    trajs = squeeze(good_examples(dex, 2, 2, :, :, :));
    xlim([0 tvec(end)]); ylim([0 1])
    for j=1:10
        subplot(4,4,7); hold on; set(gca, 'fontsize',20);    
        plot(tvec,trajs(1,:,j),'Color', colo2{1}, 'linewidth',1)
        plot(tvec,trajs(2,:,j),'Color', colo2{2}, 'linewidth',1)
        plot(tvec,trajs(3,:,j),'Color', colo2{3}, 'linewidth',1)
        plot(tvec,trajs(4,:,j),'Color', colo2{4}, 'linewidth',1)

        subplot(4,4,11); hold on; set(gca, 'fontsize',20);    
        plot(tvec, trajs(1,:,j) - trajs(4,:,j), 'color',[.7 .7 .7], 'linewidth',1)

        subplot(4,4,15); hold on; set(gca, 'fontsize',20);    
        plot(tvec, .5*(trajs(1,:,j)+trajs(4,:,j))-.5*(trajs(2,:,j)+trajs(3,:,j)), 'color',[.7 .7 .7], 'linewidth',1)
    end


    subplot(4,4,8); hold on; set(gca, 'fontsize',20);
    plot([target_on target_on], [-1 1], 'k--')
    plot([target_off target_off], [-1 1], 'k--')
    dex = find(cluster_ids == ids(i),1);
    trajs = squeeze(good_examples(dex, 3, 2, :, :, :));
    xlim([0 tvec(end)]); ylim([0 1])
    for j=1:10
        subplot(4,4,8); hold on; set(gca, 'fontsize',20);    
        plot(tvec,trajs(1,:,j),'Color', colo2{1}, 'linewidth',1)
        plot(tvec,trajs(2,:,j),'Color', colo2{2}, 'linewidth',1)
        plot(tvec,trajs(3,:,j),'Color', colo2{3}, 'linewidth',1)
        plot(tvec,trajs(4,:,j),'Color', colo2{4}, 'linewidth',1)

        subplot(4,4,12); hold on; set(gca, 'fontsize',20);    
        plot(tvec, trajs(1,:,j) - trajs(4,:,j), 'color',[.7 .7 .7], 'linewidth',1)

        subplot(4,4,16); hold on; set(gca, 'fontsize',20);    
        plot(tvec, .5*(trajs(1,:,j)+trajs(4,:,j))-.5*(trajs(2,:,j)+trajs(3,:,j)), 'color',[.7 .7 .7], 'linewidth',1)
    end


%%%%% plot choice 
    subplot(4,4,9); hold on; set(gca, 'fontsize',20);
    plot(tvec, cluster_dex(1,:), 'k', 'linewidth',2)
    ylabel({'Choice Index';'(Pro-Anti)'}, 'fontsize',20)
    xlim([0 tvec(end)])
    ylim([-1 1])

    subplot(4,4,10); hold on; set(gca, 'fontsize',20);
    plot(tvec, cluster_dex(2,:), 'k', 'linewidth',2)
    xlim([0 tvec(end)])
    ylim([-1 1])

    subplot(4,4,11); hold on; set(gca, 'fontsize',20);
    plot(tvec, cluster_dex(3,:), 'k', 'linewidth',2)
    xlim([0 tvec(end)])
    ylim([-1 1])

    subplot(4,4,12); hold on; set(gca, 'fontsize',20);
    plot(tvec, cluster_dex(4,:), 'k', 'linewidth',2)
    xlim([0 tvec(end)])
    ylim([-1 1])

%%%%% plot rule index
    subplot(4,4,13); hold on; set(gca, 'fontsize',20);
    plot(tvec, rule_dex(1,:), 'k', 'linewidth',2)
    ylabel({'Rule Index';'(Pro-Anti)'}, 'fontsize',20)
    xlim([0 tvec(end)])
    ylim([-1 1])
    xlabel('time(s)')

    subplot(4,4,14); hold on; set(gca, 'fontsize',20);
    plot(tvec, rule_dex(2,:), 'k', 'linewidth',2)
    xlim([0 tvec(end)])
    ylim([-1 1])
    xlabel('time(s)')

    subplot(4,4,15); hold on; set(gca, 'fontsize',20);
    plot(tvec, rule_dex(3,:), 'k', 'linewidth',2)
    xlim([0 tvec(end)])
    ylim([-1 1])
    xlabel('time(s)')

    subplot(4,4,16); hold on; set(gca, 'fontsize',20);
    plot(tvec, rule_dex(4,:), 'k', 'linewidth',2)
    xlim([0 tvec(end)])
    ylim([-1 1])
    xlabel('time(s)')


    fig = gcf;
    set(fig, 'PaperUnits', 'inches')
    set(fig, 'PaperPosition', [0 0 17 12])
    set(fig, 'PaperPositionMode', 'Manual')
    set(fig, 'PaperSize', [17 12])
    set(gca, 'fontsize',20)
    print(['cfig_example_' num2str(i)] ,'-dpdf')
end







close all;
