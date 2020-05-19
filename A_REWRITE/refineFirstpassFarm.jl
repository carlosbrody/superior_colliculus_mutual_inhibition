include("startup.jl")
include("commonSetup.jl")

extra_pars[:nPro]  = extra_pars[:many_trials]
extra_pars[:nAnti] = extra_pars[:many_trials]

func =  x -> JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=false,
make_dict(args, x, merge(mypars, extra_pars))...)[1]

extra_pars[:seedrand] = Int64(my_run_number*round(time()*1000))
Random.seed!(extra_pars[:seedrand])
args, seed, bounder = argsSeedBounder();
G = readCostsFile(maxCost=extra_pars[:first_pass_cost_threshold])
seed = Array{Float64}(G[my_run_number+1,3:end])

bfunc, old2new, new2old = boundedVersion(func, bounder)

g = x -> ForwardDiff.gradient(bfunc, x)
h = x -> ForwardDiff.hessian(bfunc, x)

ResultsStash = ["seedrand"             "cost"      "minimizer" ;
extra_pars[:seedrand] func(seed)  [seed]]

println("\n\n**************************\n\n")
println("\nStarting seed on for ResultsStash:")
display(ResultsStash)
println("\nInitial cost is ", func(seed))
println()



result = optimize(
    bfunc, g, h, # seems overall faster without?
    old2new(seed), NewtonTrustRegion(), # seems overall faster without?
    Optim.Options(store_trace=true, show_trace=true,
    iterations=40);
    inplace=false);

truecost = func(new2old(Optim.minimizer(result)))
hostname = chomp(read(`hostname`, String))
fp = open("negRefinedCosts_$hostname.csv", "a")
println(fp, extra_pars[:seedrand], ", ", truecost, ", ")
writedlm(fp, new2old(Optim.minimizer(result))[:]', ',')
close(fp)

println("With seedrand ", extra_pars[:seedrand], " true cost was ", truecost)

ResultsStash = vcat(ResultsStash,
[extra_pars[:seedrand] Optim.minimum(result) [new2old(Optim.minimizer(result))]])

println("\nResultsStash:")
display(ResultsStash)
println()

println("\nJJ:")
display(JJ(extra_pars[:nPro], extra_pars[:nAnti]; asDict=true, verbose=false,
    make_dict(args, new2old(Optim.minimizer(result)), merge(mypars, extra_pars))...))

println()
