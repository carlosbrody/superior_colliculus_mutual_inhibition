# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.




# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


using JLD

nruns = 0;
for f in filter(x -> startswith(x, "farm_C10_"), readdir("../NewFarms")); nruns += 1; end
@printf("%d runs found\n", nruns)

outs = zeros(nruns,3)
for i=1:nruns;
    myname = @sprintf("%d", i)
    while length(myname)<4; myname = "0" * myname; end
    print()
    outs[i,1] = load("../NewFarms/farm_C6_" * myname * ".jld", "traj1")[2,end]
    traj3, qu_out = load("../NewFarms/farm_C10_" * myname * ".jld", "traj3", "qu_out")
    outs[i,2] = traj3[2,end]
    outs[i,3] = qu_out[1]
    if rem(i,10)==0; @printf("  did run %d\n", i); end
end


# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


using StatsBase
using PyPlot
pygui(true)

x = -0.001:0.000001:maximum(outs[:,1:2])

h1 = fit(Histogram, outs[:,1], x, closed=:right)
h2 = fit(Histogram, outs[:,2], x, closed=:right)

figure(2); clf();
subplot(2,1,1)
plot(x, [0 ; cumsum(h1.weights)]/sum(h1.weights), "b-")
plot(x, [0 ; cumsum(h2.weights)]/sum(h2.weights), "r-")

legend(["standard minimization", "with pre-search"])
title(@sprintf("%d runs for each method", sum(h1.weights)))

subplot(2,1,2)
plot(x, [0 ; cumsum(h1.weights)]/sum(h1.weights), "b-")
plot(x, [0 ; cumsum(h2.weights)]/sum(h2.weights), "r-")
xlim([-0.001, -0.0007])
vlines(-0.00095, ylim()[1], ylim()[2], color="g", linestyle="--")

xlabel("Training cost (left of green dashed is a good run)")
ylabel("fraction of runs")


pre_fail = outs[outs[:,3].==0,2]
pre_succ = outs[.!(outs[:,3].==0),2]

text(-0.000925, 0.85, 
    @sprintf("Success rate after pre-search fail is %.2f%%\n", 100*length(find(pre_fail.<-0.0009))/length(pre_fail)))
text(-0.000925, 0.7, 
    @sprintf("Success rate after pre-search success is %.2f%%\n", 100*length(find(pre_succ.<-0.0009))/length(pre_succ)))




