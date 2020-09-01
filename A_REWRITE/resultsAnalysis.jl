include("startup.jl")

using DelimitedFiles


"""

answer = runSolution(id, ntrials)

Loads Solutions/solutions6.mat, and from there picks out a solution and runs
ntrials of it.  Returns a dictionary, here referred to as "answers".

Unit activities will be in answers["proValls"] (all the Pro trials) and
answers["antiValls"] (all the Anti trials). Each one of those will have size
nOptoConditions by nRunConditions where the nRunConditions will be the
product of the number of rule_and_delay_periods and the target_periods.
(values of these two for each run are found in answers["rule_and_delay_periods"]
and in answers["target_periods"]).

So, for example, answers["proValls"][1,3] will be unit activities for the first
opto condition (control), and for answers["rule_and_delay_period"][3] and
answers["target_period"][3]. Opto conditions are 1=control, 2=delay_period,
3=target_period.
    answers["hBP"] will also be nOptoConditions by nRunConditions, and holds
        fraction of hits for Pro trials
    answers["hBA"] will also be nOptoConditions by nRunConditions, and holds
        fraction of hits for Anti trials

answers["proValls"][i,j] will itself have size nunits-by-ntimebins-by-ntrials,
where nunits=6.

NodeIDs are indicated by answers["ProNodeID"], answers["AntiNodeID"],
answers["LeftNodeID"] and answers["RightNodeID"]. answers[""]

answers["paramNames"]      contains parameter names.
answers["paramvals"]       contains parameter values for this solution,
answers["weightMatrixMap"] contains strings showing weights
    Weights are symmetrical across the two sides of the brain, such that
    node 2 (the first Anti node) will have the same weight pattern as node 4
    (the first Anti node on the other side of the brain).


"""
function runSolution(id::Int64, ntrials::Int64)

    G = matread("Solutions/solutions6.mat")
    pvals  = G["paramvals"]
    args   = G["argnames"]
    srands = G["randseeds"]

    if id < 0  || id > length(srands)
        error("id must be between 1 and $(length(srands)), the number of solutions.")
    end
    id = size(pvals,1)-(id-1)  # We choose a number between 1 and size(pvals,1) to identify
    # the solution in pvals that we'll run with

    params = pvals[id,:]
    include("sixNodeSetup_C32.jl")
    extra_pars[:seedrand] = srands[id]
    extra_pars[:few_trials]                = 50
    extra_pars[:many_trials]               = 1600

    extra_pars[:nPro] = ntrials
    extra_pars[:nAnti] = ntrials

    # This will now run the same way as when the cost during training was calculated,
    # meaning that costs[id] and answers["cost"] will be the same if we run with
    # extra_pars[:many_trials] trials:
    answers = JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=true, asDict=true,
        model_details=true, merge(make_dict(args, params, mypars), extra_pars)...)

    answers["AntiNodeID"]    =          mypars[:AntiNodeID];
    answers["ProNodeID"]     =          mypars[:ProNodeID];
    answers["RightNodeID"]   =          mypars[:RightNodeID];
    answers["LeftNodeID"]    =          mypars[:LeftNodeID];

    answers["rule_and_delay_periods"] = mypars[:rule_and_delay_periods][:][[1,1,2,2]]
    answers["target_periods"]         = mypars[:target_periods][:][[1,2,1,2]]

    answers["weightMatrixMap"] = G["weightMatrixMap"]
    answers["paramvals"]       = G["paramvals"]
    answers["paramNames"]      = G["argnames"]

    return answers
end



"""
findSolutionsFromReportsFiles()

Looks in ../../Reports/ for files that start with "r6_spockRun3", and parses
those to find those that went into further trainibg (phase 2), and from those,
the latest cost and parameters; stashes any that have cost below cost_threshold
in a .mat file Solutions/solution6.mat
"""
function findSolutionsFromReportsFiles(; cost_threshold = -0.0001)

    include("sixNodeSetup_C32.jl")
    args = argsSeedBounder()[1]   # get the string list of model params

    pullReports=false; if pullReports   # if there are new runs on the VMs, get them
        remotes = "spock"

        if remotes == "googleVMs"
            machines = ["proanti006", "proanti007", "proanti008",
                        "proanti009", "proanti010", "proanti011"]

            for i in machines
                run(`./pull_reports.sh $i`)
            end
        elseif remotes == "spock"
            run(`./pull_reports.sh`)
        else
            error("Don't know remotes $remotes")
        end

        # At the of this block, report files will have been put in ../../Reports
    end

    global costs  = fill(0.0, 0)   # the cost for each solution
    global fnames = fill("",  0)   # the filename for each solution
    global pvals  = fill(0.0, 0, length(args))  # model params for each solution
    global srands = fill(0, 0)     # random seeds for each solution

    # Now we'll parse each Report file to find where it ended up at.
    u = readdir("../../Reports")
    for i in filter(x -> startswith(x, "r6_spockRun3"), u)
        # Read the file into string s:
        fp = open("../../Reports/$i", "r")
        s = read(fp, String)
        close(fp)
        println(i)
        # Find the last cost report:
        z = findlast("OVERALL", s)
        if z != nothing  && findlast("random seed", s)!=nothing
            z2 = findfirst("cost=", s[z[end]:end])
            z3 = findfirst(",", s[z[end]:end])
            mycost = parse(Float64, s[z[end]+z2[end]:z[end]+z3[1]-2])
            # println("$i: mycost=$mycost")
            global costs  = vcat(costs, mycost)
            global fnames = vcat(fnames, i)

            z = findlast("ps=[", s)[end].+1
            n = findfirst("]", s[z+1:end])[1] + z-1
            global pvals = vcat(pvals, readdlm(Vector{UInt8}(s[z:n]), ','))

            sz1 = findlast("random seed", s)
            sz2 = findfirst('\n', s[sz1[end]:end])
            global srands = vcat(srands, parse(Int64, s[sz1[end]+2:sz1[end]+sz2-2]))
        end
        # println("Checked $i")
    end

    println("$(length(costs)) networks went into further training")

    p = sortperm(costs, rev=true)
    costs = costs[p]
    fnames = fnames[p]
    pvals = pvals[p,:]
    srands = srands[p]

    u = findall(costs .<= cost_threshold)
    println("\nFound $(length(u)) networks below cost threshold of $cost_threshold")

    matwrite("Solutions/solutions6.mat", Dict(
        "costs"=>costs[u],
        "fnames"=>fnames[u],
        "argnames"=>args,
        "paramvals"=>pvals[u,:],
        "randseeds"=>srands[u],
        "weightMatrixMap"=>makeWeightMatrix(mypars; asString=true)
        ))
end


##
