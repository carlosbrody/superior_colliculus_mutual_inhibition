include("startup.jl")
include("commonSetup.jl")

# extra_pars[:seedrand] = Int64(my_run_number*round(time()*1000))
#Random.seed!(extra_pars[:seedrand])
# args, seed, bounder = argsSeedBounder();

extra_pars[:nPro]  = extra_pars[:few_trials]
extra_pars[:nAnti] = extra_pars[:few_trials]

func =  x -> JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=false,
    make_dict(args, x, merge(mypars, extra_pars))...)[1]

for looper=1:400

    extra_pars[:seedrand] = Int64(my_run_number*round(time()*1000))
    Random.seed!(extra_pars[:seedrand])
    args, seed, bounder = argsSeedBounder();

    bfunc, old2new, new2old = boundedVersion(func, bounder)

    g = x -> ForwardDiff.gradient(bfunc, x)
    h = x -> ForwardDiff.hessian(bfunc, x)

    ResultsStash = ["seedrand"             "cost"      "minimizer" ;
                     extra_pars[:seedrand] func(seed)  [seed]]

    println("\n\n**************************\n\n")
    println("\nStarting seed on loop $looper for ResultsStash:")
    display(ResultsStash)
    println("\nInitial cost is ", func(seed))
    println()



    result = optimize(
        bfunc, g, h, # seems overall faster without?
        old2new(seed), NewtonTrustRegion(), # seems overall faster without?
        Optim.Options(store_trace=true, show_trace=true,
            iterations=500);
        inplace=false);

    truecost = func(new2old(Optim.minimizer(result)))
    if truecost < extra_pars[:first_pass_cost_threshold]
        println("------> BELOW THRESHOLD COST <-------")
        hostname = chomp(read(`hostname`, String))
        fp = open("negCosts_$hostname.csv", "a")
        println(fp, extra_pars[:seedrand], ", ", truecost, ", ")
        writedlm(fp, new2old(Optim.minimizer(result))[:]', ',')
        close(fp)
    end

    end
    println("With seedrand ", extra_pars[:seedrand], " true cost was ", truecost)

    ResultsStash = vcat(ResultsStash,
        [extra_pars[:seedrand] Optim.minimum(result) [new2old(Optim.minimizer(result))]])

    println("\nResultsStash:")
    display(ResultsStash)
    println()

    println("\nJJ:")
    display(JJ(extra_pars[:nPro], extra_pars[:nAnti]; asDict=true, verbose=false,
        make_dict(G["args"], new2old(Optim.minimizer(result)), merge(mypars, extra_pars))...))

    println()
end
