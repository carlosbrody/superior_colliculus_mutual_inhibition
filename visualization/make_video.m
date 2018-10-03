close all;
clear all;

% load cluster_ids, examples, results
load MiniC32_C32_examples_50.mat

% get data with L/R synthetic trials
load MiniC32_C32_examples_50_LR.mat

% filter out bad solutions
dex = [results.cost{:}] <= -0.0001;
ex = examples(dex, :,:,:,:,:);


% some settings
p.plot_hits = true;     % if true, plot hits only. If false, plot errors only
p.plot_traj = true;    %if false, only plot the final values
farm_id = 1;        % which solution to plot
run_solution(farm_id, ex,p)

