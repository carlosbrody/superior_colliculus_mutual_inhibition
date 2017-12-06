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

# 
cost_choice = "cost"
threshold   = -0.0001

paramso = reso["params"]
paramso = paramso - ones(size(paramso,1),1)*mean(paramso,1)
paramso = paramso ./ (ones(size(paramso,1),1)*std(paramso,1))

C = (paramso' * paramso)/size(paramso,1)
D, V = eig(C)

oVparams = (inv(V)*reso["params"]')'[:,end:-1:1]
nVparams = (inv(V)*resn["params"]')'[:,end:-1:1]

figure(1); clf(); iold = inew = 0

u = find(res[cost_choice].<=threshold)
plot(oVparams[u,1], oVparams[u,2], "r.", markersize=11)

nloop = hg = hr = 1
for f in readdir("MiniOptimized_Redux/")
    if startswith(f, "farm_C17") # && nloop > 10
        iold = find(reso["files"] .== ("MiniOptimized/" * f))
        inew = find(resn["files"] .== ("MiniOptimized_Redux/" * f))
        hg = plot([oVparams[iold,1], nVparams[inew,1]], [oVparams[iold,2], nVparams[inew,2]], 
            "g.-", markersize=11)
        hr = plot(oVparams[iold,1], oVparams[iold,2], "r.-", markersize=11)
    end
    nloop += 1
end

xlabel("PCA Dim 1")
ylabel("PCA Dim 2")
legend([hg[1], hr[1]], ["after re-optimizing with new test trials", "training trials"])


