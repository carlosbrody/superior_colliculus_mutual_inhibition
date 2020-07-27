include("startup.jl")
include("commonSetup.jl")

using Dates

# extra_pars[:seedrand] = Int64(my_run_number*round(time()*1000))
#Random.seed!(extra_pars[:seedrand])

args, seed, bounder = argsSeedBounder();
extra_pars[:nPro]  = extra_pars[:pass1NTrials]
extra_pars[:nAnti] = extra_pars[:pass1NTrials]

global latest_minimizer
global latest_minimum = Inf
##
func = x -> begin
    answer = JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=false,
        make_dict(args, x, merge(mypars, extra_pars))...)[1]
    if answer < latest_minimum
        global latest_minimum = answer
        global latest_minimizer = x
    end
    return answer
end

##

function passNCallback(x, n)
    if typeof(x) <: Array
        x = x[end]
    end
    println("PASS $n.$(x.iteration) "*
        "$(Dates.format(now(), "e, dd u yyyy HH:MM:SS")):"*
        " pars=", get_value(latest_minimizer))
    return x.value <= extra_pars[Symbol("pass$(n)CostThreshold")]
end

extra_pars[:seedrand] = Int64(my_run_number*round(time()*1000))
Random.seed!(extra_pars[:seedrand])
global seed
args, seed, bounder = argsSeedBounder();

bfunc, old2new, new2old = boundedVersion(func, bounder)

g = x -> ForwardDiff.gradient(bfunc, x)
h = x -> ForwardDiff.hessian(bfunc, x)

ResultsStash = ["seedrand"             "cost"      "minimizer" ;
    extra_pars[:seedrand] func(seed)  [seed]]

"""
    passNOptimization(seed, n)

    Given a seed and a pass number, runs the optimization, and
    returns the result, the final cost, and a Boolean, stopFlag,
    which is true if the pass's cost threshold was below the final cost.
"""
function passNOptimization(seed, n)
    mypars[:rule_and_delay_periods] = extra_pars[Symbol("pass$(n)RnD")]

    extra_pars[:nPro]  = extra_pars[Symbol("pass$(n)NTrials")]
    extra_pars[:nAnti] = extra_pars[Symbol("pass$(n)NTrials")]

    println("\nWill now do $(extra_pars[Symbol("pass$(n)NIter")]) iterations "*
        "with $(extra_pars[:nPro]) trials, up to a threshold cost " *
        "$(extra_pars[Symbol("pass$(n)CostThreshold")])." )
    println("rule and delay period range is $(extra_pars[Symbol("pass$(n)RnD")]).\n")

    result = optimize(
        bfunc, g, h,
        old2new(seed), NewtonTrustRegion(),
        Optim.Options(store_trace=true, show_trace=true,
            callback  = x -> passNCallback(x, n),
            iterations=extra_pars[Symbol("pass$(n)NIter")]);
        inplace=false);

    stopFlag = !evaluateModel(extra_pars[:testruns],
        args, new2old(Optim.minimizer(result)), mypars, extra_pars,
        reportFile=stdout)

    if !stopFlag
        println("------> BELOW REQUIREMENTS for PASS $n <-------")
        hostname = chomp(read(`hostname`, String))
        fp = open("negPass$(n)Costs_$hostname.csv", "a")
        print(fp, Int64(extra_pars[:seedrand]), ", ", truecost, ", ")
        writedlm(fp, new2old(Optim.minimizer(result))[:]', ',')
        close(fp)
    end

    return result, truecost, stopFlag
end

##

global truecost, seed, result

println("\n\n**************************\n\n")
println("\n$(Dates.format(now(), "e, dd u yyyy HH:MM:SS")) : "*
    "Starting seedrand $(extra_pars[:seedrand]) on main loop for ResultsStash:")
display(ResultsStash)


for passnum=1:extra_pars[:nPasses]
    global truecost, seed, result

    result, truecost, stopFlag = passNOptimization(seed, passnum)
    if stopFlag; break; end

    seed = new2old(Optim.minimizer(result))
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
