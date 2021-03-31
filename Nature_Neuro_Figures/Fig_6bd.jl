# This script makes Figure 6 B,D
# This script runs in Julia 0.6.4

# load relevant files
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")
include("unilateral_analysis.jl")
include("cluster_farms.jl")
include("cluster_example.jl")

# load results of entire farm run
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)
full_output = load("MiniC32_C32_full_trial_inactivation_2.jld","output");

# Making psychometric with relative error increase
## Delay period
targets = [90 70; 85 50; 90 70];
diff_targets = [5 20; 0 0];
mean_diff_pro_delay = mean(full_output[:,4,1] - full_output[:,3,1]).*-100;
mean_diff_anti_delay = mean(full_output[:,4,2] - full_output[:,3,2]).*-100;
std_diff_pro_delay = std(full_output[:,4,1] - full_output[:,3,1]).*-100;
std_diff_anti_delay = std(full_output[:,4,2] - full_output[:,3,2]).*-100;
figure(figsize=(1.5,2))
plot(vec([1 1.5]), vec([diff_targets[1,1] diff_targets[1,1]]),"k")
plot(vec([2 2.5]), vec([diff_targets[1,2] diff_targets[1,2]]),"k")
plot(vec([1.25]), mean_diff_pro_delay, "o",color=(0, 128/255, 0))
plot(vec([2.25]), mean_diff_anti_delay, "o", color=(253/255,137/255,57/255))
plot(vec([1.25 1.25]), vec([mean_diff_pro_delay + 1.96*[-std_diff_pro_delay, std_diff_pro_delay]][1]),color=(0, 128/255, 0))
plot(vec([2.25 2.25]), vec([mean_diff_anti_delay + 1.96*[-std_diff_anti_delay, std_diff_anti_delay]][1]), color=(253/255,137/255,57/255))
plot(vec([0 3]), vec([0 0]), color=(179/255,179/255,179/255))
ylim(-5, 25)
xlim(0.5, 3)
ylabel("% Error Increase")
xlabel("Delay")
xticks([1.25, 2.25], ["Pro", "Anti"])

## Choice period
mean_diff_pro_choice = mean(full_output[:,5,1] - full_output[:,3,1]).*-100;
mean_diff_anti_choice = mean(full_output[:,5,2] - full_output[:,3,2]).*-100;
std_diff_pro_choice = std(full_output[:,5,1] - full_output[:,3,1]).*-100;
std_diff_anti_choice = std(full_output[:,5,2] - full_output[:,3,2]).*-100;
figure(figsize=(1.5,2))
plot(vec([1 1.5]), vec([diff_targets[2,1] diff_targets[2,1]]),"k")
plot(vec([2 2.5]), vec([diff_targets[2,2] diff_targets[2,2]]),"k")
plot(vec([1.25]), mean_diff_pro_choice, "o",color=(0, 128/255, 0))
plot(vec([2.25]), mean_diff_anti_choice, "o", color=(253/255,137/255,57/255))
plot(vec([1.25 1.25]), vec([mean_diff_pro_choice + 1.96*[-std_diff_pro_choice, std_diff_pro_choice]][1]),color=(0, 128/255, 0))
plot(vec([2.25 2.25]), vec([mean_diff_anti_choice + 1.96*[-std_diff_anti_choice, std_diff_anti_choice]][1]), color=(253/255,137/255,57/255))
plot(vec([0 3]), vec([0 0]), color=(179/255,179/255,179/255))
ylim(-5, 25)
xlim(0.5, 3)
ylabel("% Error Increase")
xlabel("Choice")
xticks([1.25, 2.25], ["Pro", "Anti"])

## Full Trial
if true
    diff_pro_full  = (full_output[:,1,1]    - .9).*-100;
    diff_anti_full = (full_output[:,1,2]   - .7).*-100;
    plb = mean(diff_pro_full) - 3*std(diff_pro_full);
    pub = mean(diff_pro_full) + 3*std(diff_pro_full);
    pdex = (diff_pro_full .> plb) .& (diff_pro_full .< pub);
    alb = mean(diff_anti_full) - 3*std(diff_anti_full);
    aub = mean(diff_anti_full) + 3*std(diff_anti_full);
    adex = (diff_anti_full .> alb) .& (diff_anti_full .< aub);
    dex = adex .* pdex;
    mean_diff_pro_full = mean(full_output[dex,1,1]    - .9).*-100;
    mean_diff_anti_full = mean(full_output[dex,1,2]   - .7).*-100;
    std_diff_pro_full = std(full_output[dex,1,1]      - .9).*-100;
    std_diff_anti_full = std(full_output[dex,1,2]     - .7).*-100;
else
    mean_diff_pro_full = mean(full_output[:,1,1]    - .9).*-100;
    mean_diff_anti_full = mean(full_output[:,1,2]   - .7).*-100;
    std_diff_pro_full = std(full_output[:,1,1]      - .9).*-100;
    std_diff_anti_full = std(full_output[:,1,2]     - .7).*-100;
end
figure(figsize=(1.5,2))
plot(vec([1.25]), mean_diff_pro_full, "o",color=(0, 128/255, 0))
plot(vec([2.25]), mean_diff_anti_full, "o", color=(253/255,137/255,57/255))
plot(vec([1.25 1.25]), vec([mean_diff_pro_full  + 1.96*[-std_diff_pro_full,  std_diff_pro_full]][1]),color=(0, 128/255, 0))
plot(vec([2.25 2.25]), vec([mean_diff_anti_full + 1.96*[-std_diff_anti_full, std_diff_anti_full]][1]), color=(253/255,137/255,57/255))
plot(vec([0 3]), vec([0 0]), color=(179/255,179/255,179/255))
ylim(-15, 35)
xlim(0.5, 3)
ylabel("% Error Increase")
xlabel("Full trial")
xticks([1.25, 2.25], ["Pro", "Anti"])

## Rule period
if true
    diff_pro_rule  = (full_output[:,2,1]    - .9).*-100;
    diff_anti_rule = (full_output[:,2,2]   - .7).*-100;
    plb = mean(diff_pro_rule) - 3*std(diff_pro_rule);
    pub = mean(diff_pro_rule) + 3*std(diff_pro_rule);
    pdex = (diff_pro_rule .> plb) .& (diff_pro_rule .< pub);
    alb = mean(diff_anti_rule) - 3*std(diff_anti_rule);
    aub = mean(diff_anti_rule) + 3*std(diff_anti_rule);
    adex = (diff_anti_rule .> alb) .& (diff_anti_rule .< aub);
    dex = adex .* pdex;
    mean_diff_pro_rule = mean(full_output[dex,2,1]    - .9).*-100;
    mean_diff_anti_rule = mean(full_output[dex,2,2]   - .7).*-100;
    std_diff_pro_rule = std(full_output[dex,2,1]      - .9).*-100;
    std_diff_anti_rule = std(full_output[dex,2,2]     - .7).*-100;
else
    mean_diff_pro_rule = mean(full_output[:,2,1]    - output[:,1,1]).*-100;
    mean_diff_anti_rule = mean(full_output[:,2,2]   - output[:,1,2]).*-100;
    std_diff_pro_rule = std(full_output[:,2,1]      - output[:,1,1]).*-100;
    std_diff_anti_rule = std(full_output[:,2,2]     - output[:,1,2]).*-100;
end

figure(figsize=(1.5,2))
plot(vec([1.25]), mean_diff_pro_rule, "o",color=(0, 128/255, 0))
plot(vec([2.25]), mean_diff_anti_rule, "o", color=(253/255,137/255,57/255))
plot(vec([1.25 1.25]), vec([mean_diff_pro_rule  + 1.96*[-std_diff_pro_rule,  std_diff_pro_rule]][1]),color=(0, 128/255, 0))
plot(vec([2.25 2.25]), vec([mean_diff_anti_rule + 1.96*[-std_diff_anti_rule, std_diff_anti_rule]][1]), color=(253/255,137/255,57/255))
plot(vec([0 3]), vec([0 0]), color=(179/255,179/255,179/255))
ylim(-15, 35)
xlim(0.5, 3)
ylabel("% Error Increase")
xlabel("Task cue")
xticks([1.25, 2.25], ["Pro", "Anti"])


