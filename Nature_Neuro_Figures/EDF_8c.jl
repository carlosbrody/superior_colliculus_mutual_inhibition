# This script makes EDF panel 8c
# The variance explained within each solution.

# load relevant files
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")
include("unilateral_analysis.jl")
include("cluster_farms.jl")
include("cluster_example.jl")

# Make examples if you haven't already
make_examples()

# Do analysis
dims = plot_dimension_analysis(;return_all=true,limited_delay=true);
figure(figsize=(3,3))
means = mean(dims,1)[1:12].*100;
stds  = std(dims,1)[1:12].*100;
plot(vec([0.5 3.5]), vec([86 86]), color=(179/255,179/255,179/255));
plot(vec([0.5 3.5]), vec([97 97]), color=(179/255,179/255,179/255));
plot(vec([1 2 3]), vec(means[1:4:end]),"ko-")
plot(vec([1 2 3]), vec(means[2:4:end]),"ro-")
plot(vec([1 2 3]), vec(means[3:4:end]),"bo-")
plot(vec([1,1]), vec([means[1]+stds[1], means[1]-stds[1]]), "k");
plot(vec([1,1]), vec([means[2]+stds[2], means[2]-stds[2]]), "r");
plot(vec([1,1]), vec([means[3]+stds[3], means[3]-stds[3]]), "b");
plot(vec([2,2]), vec([means[5]+stds[5], means[5]-stds[5]]), "k");
plot(vec([2,2]), vec([means[6]+stds[6], means[6]-stds[6]]), "r");
plot(vec([2,2]), vec([means[7]+stds[7], means[7]-stds[7]]), "b");
plot(vec([3,3]), vec([means[9]+stds[9], means[9]-stds[9]]), "k");
plot(vec([3,3]), vec([means[10]+stds[10], means[10]-stds[10]]), "r");
plot(vec([3,3]), vec([means[11]+stds[11], means[11]-stds[11]]), "b");

xlim(0.5,3.5); 
ylim(0,100);
ylabel("Variance Explained %",fontsize=16)
yticks(fontsize=12)
xticks([1,2,3],["Full","Delay","Target"],fontsize=12)
plt[:tight_layout]()

