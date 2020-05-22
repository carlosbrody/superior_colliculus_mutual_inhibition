include("startup.jl")
include("commonSetup.jl")

# extra_pars[:seedrand] = Int64(my_run_number*round(time()*1000))
#Random.seed!(extra_pars[:seedrand])

args, seed, bounder = argsSeedBounder();
extra_pars[:nPro]  = extra_pars[:few_trials]
extra_pars[:nAnti] = extra_pars[:few_trials]

func =  x -> JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=false,
    make_dict(args, x, merge(mypars, extra_pars))...)[1]


function firstPassCallback(x)
    if typeof(x) <: Array
        x = x[end]
    end
    return x.value <= extra_pars[:first_pass_cost_threshold]
end

function secondPassCallback(x)
    if typeof(x) <: Array
        x = x[end]
    end
    return x.value <= extra_pars[:stoppingCostThreshold]
end


for looper=1  # :400

    extra_pars[:nPro]  = extra_pars[:few_trials]
    extra_pars[:nAnti] = extra_pars[:few_trials]
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
        bfunc, g, h,
        old2new(seed), NewtonTrustRegion(),
        Optim.Options(store_trace=true, show_trace=true,
        callback  = firstPassCallback,
        iterations=extra_pars[:firstPassNIter]);
        inplace=false);

    truecost = func(new2old(Optim.minimizer(result)))
    if truecost < extra_pars[:first_pass_cost_threshold]
        println("------> BELOW THRESHOLD FIRST PASS COST <-------")
        hostname = chomp(read(`hostname`, String))
        fp = open("neg50Costs_$hostname.csv", "a")
        print(fp, Int64(extra_pars[:seedrand]), ", ", truecost, ", ")
        writedlm(fp, new2old(Optim.minimizer(result))[:]', ',')
        close(fp)

        extra_pars[:nPro]  = extra_pars[:many_trials]
        extra_pars[:nAnti] = extra_pars[:many_trials]
        seed = new2old(Optim.minimizer(result))

        result = optimize(
            bfunc, g, h, # seems overall faster without?
            old2new(seed), NewtonTrustRegion(), # seems overall faster without?
            Optim.Options(store_trace=true, show_trace=true,
            callback = secondPassCallback,
            iterations=extra_pars[:secondPassNIter]);
            inplace=false);

        truecost = func(new2old(Optim.minimizer(result)))
        hostname = chomp(read(`hostname`, String))
        costsfile = "neg1600Costs_$hostname.csv"
        fp = open(costsfile, "a")
        print(fp, Int64(extra_pars[:seedrand]), ", ", truecost, ", ")
        writedlm(fp, new2old(Optim.minimizer(result))[:]', ',')
        close(fp)

    end
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
end
