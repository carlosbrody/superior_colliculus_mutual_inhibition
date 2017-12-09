# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.


if !isdefined(:plot_PA)
    @printf("Including pro_anti.jl\n")
    include("pro_anti.jl")
end

if !isdefined(:farmload)
    @printf("Including results_analysis.jl\n")
    include("results_analysis.jl")
end

if !isdefined(:reso)
    @printf("Loading MiniOptimized farms\n")
    reso = farmload("C17"; farmdir="MiniOptimized")
end

if !isdefined(:resn) || true   # forcing it to run right now so we always get latest updated files
    @printf("Loading MiniOptimized_Redux farms\n")
    resn = farmload("C17"; farmdir="MiniOptimized_Redux")
end

# The next two lines define how we threshold the optimized runs, to choose which ones go into the PCA plot
cost_choice = "cost"     # "cost" means test cost, "tcost" means training cost
threshold   = -0.0001    # Any run with a cost below this goes in

# Next we compute PCA space.
paramso = reso["params"][find(res[cost_choice].<=threshold),:]
# Whiten (mean zero, std 1):
paramso = paramso - ones(size(paramso,1),1)*mean(paramso,1)
paramso = paramso ./ (ones(size(paramso,1),1)*std(paramso,1))
# compute and diagonalize covariance matrix:
C = (paramso' * paramso)/size(paramso,1)
D, V = eig(C)
# Project old and new parameters onto eigenvectors, reverse so top PCA is first
oVparams = (inv(V)*reso["params"]')'[:,end:-1:1]
nVparams = (inv(V)*resn["params"]')'[:,end:-1:1]

figure(100); clf(); 

# Plot a red dot for all the old parameters:
u = find(res[cost_choice].<=threshold)
plot(oVparams[u,1], oVparams[u,2], "r.", markersize=11)

hg = hr = 0  # define here so vars are available outside for loop
# Now plot a green dot and a connecting line for all the new ones:
for f in readdir("MiniOptimized_Redux/")
    if startswith(f, "farm_C17") 
        iold = find(reso["files"] .== ("MiniOptimized/" * f))
        inew = find(resn["files"] .== ("MiniOptimized_Redux/" * f))
        hg = plot([oVparams[iold,1], nVparams[inew,1]], [oVparams[iold,2], nVparams[inew,2]], 
            "g.-", markersize=11)
        # Add another red dot on top of old values to keep it red
        hr = plot(oVparams[iold,1], oVparams[iold,2], "r.-", markersize=11)
    end
end

xlabel("PCA Dim 1")
ylabel("PCA Dim 2")
legend([hg[1], hr[1]], ["after re-optimizing with new test trials", "training trials"])


