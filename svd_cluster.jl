# Set of functions for doing SVD clustering on neural trajectories

# I have no fucking idea why I need to include these, but load() won't work unless I do!
include("pro_anti.jl")
using HDF5

# load parameters and info for each farm in farmdir 
function load_farm_params(farm_id; farmdir="MiniFarms", verbose=true, verbose_every=50);
    if typeof(farmdir)==String; farmdir=[farmdir]; end
    results = Dict(); dirs=[]; files =[]; qs=[]; tcosts=[]; costs=[]; pars=[]; n=0;
    for dd in farmdir
        for f in filter(x -> startswith(x, "farm_" * farm_id * "_"), readdir(dd * "/"))
            n +=1
            myfile = dd * "/" * f; 
            args, params, traj3, cost = load(myfile, "args", "pars3", "traj3", "cost")
            files = [files ; myfile]; dirs = [dirs ; dd]
            tcosts = [tcosts; traj3[2,end]]; costs = [costs; cost]
            if length(pars)==0; pars=params'; else pars = [pars ; params']; end
            if verbose && rem(n, verbose_every)==0
                @printf("%s %g\n", myfile, tcosts[end])
            end

        end
    end

    results["dirs"] = dirs
    results["files"]  = files
    results["tcost"]  = tcosts
    results["cost"]   = costs
    results["params"] = pars

    return results
end

# returns avgerage external trajectory for each node in 4 conditions (hit/miss X pro/anti)
function run_farm(filename; testruns=400, overrideDict=Dict(),all_conditions=false)
    mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")
    numConditions = 1;
    if all_conditions
        numConditions = 3;
    end
    avgHP = [];
    avgMP = [];
    avgHA = [];
    avgMA = [];
    for period = 1:numConditions
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
        :opto_times=>reshape(extra_pars[:opto_periods][period,:], 1, 2),
        ))

        proVs, antiVs, proF, antiF = run_ntrials(testruns, testruns; merge(make_dict(args, pars3, these_pars), overrideDict)...);

        hitBP = find(proVs[1,:]  .> proVs[4,:])
        avghp = mean(proF[:,:,hitBP],3);

        missBP= find(proVs[1,:]  .< proVs[4,:])
        avgmp = mean(proF[:,:,missBP],3);

        hitBA = find(antiVs[4,:] .> antiVs[1,:])
        avgha = mean(antiF[:,:,hitBA],3);

        missBA= find(antiVs[4,:] .< antiVs[1,:])              
        avgma = mean(antiF[:,:,missBA],3);
    
        avgHP = [avghp[1,:,1]; avghp[2,:,1]; avghp[3,:,1]; avghp[4,:,1]];
        avgMP = [avgmp[1,:,1]; avgmp[2,:,1]; avgmp[3,:,1]; avgmp[4,:,1]];
        avgHA = [avgha[1,:,1]; avgha[2,:,1]; avgha[3,:,1]; avgha[4,:,1]];
        avgMA = [avgma[1,:,1]; avgma[2,:,1]; avgma[3,:,1]; avgma[4,:,1]];
    end
    return avgHP, avgMP, avgHA, avgMA 
end

function build_response_matrix(results)
    response_matrix = [];
    for i = 1:length(results["cost"])
        filename = results["files"][i];
        avgHP, avgMP, avgHA, avgMA  = run_farm(filename);
        if isempty(response_matrix)
            response_matrix = [avgHP' avgMP' avgHA' avgMA'];
        else
            response_matrix = [response_matrix; avgHP' avgMP' avgHA' avgMA'];
        end
        @printf("%g %s %g\n",i, filename, results["tcost"][i])
    end

    myfilename = "SVD_response_matrix.jld";
    save(myfilename, Dict("response"=>response_matrix, "results"=>results))

    return response_matrix
end


function SVD_analysis()
    # Load responses from all models
    response, results = load("SVD_response_matrix.jld", "response","results")


    # need to filter out NaN rows (not sure why some are bad! Maybe make no mistakes?)
    nanrows = any(isnan(response),2);
    r_all = response[!vec(nanrows),:];

    # Sort out responses in each condition
    # Each is 4 nodes x 76 timesteps = 304 
    # node 1, 2, 3, 4 is the same ordering as run_nTrials
    r_hp = r_all[:,1:304];
    r_mp = r_all[:,304*2+1:2*304];
    r_ha = r_all[:,304*2+1:3*304];
    r_ma = r_all[:,304*3+1:end];

    # F[:U], F[:S], F[:Vt]
    F_all = svdfact(r_all);

end

function plot_SVD_approx(rank, condition, F)
    figure()
    S = copy(F[:S]);
    S[rank+1:end] = 0;    
    low_rank = F[:U] * Diagonal(S)* F[:Vt];
    U = copy(F[:U])
    Vt = copy(F[:Vt])

    neural_dex = 1;
    if condition     == "hit-pro"
        cdex = 0;
    elseif condition == "miss-pro"
        cdex = 1;
    elseif condition == "hit-anti"
        cdex = 2;
    elseif condition == "miss-anti"
        cdex = 3;
    else
        @printf("Unknown Condition: %s \n", condition)
    end

    numts = Int(size(Vt,2)/4/4);
    numtsC= Int(size(Vt,2)/4);
    node1 = copy(Vt[neural_dex,cdex*numtsC+1:cdex*numtsC+numts]);
    node2 = copy(Vt[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts]);
    node3 = copy(Vt[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts]);
    node4 = copy(Vt[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts]);
    plot(-node1)
    plot(-node2)
    plot(-node3)
    plot(-node4)
    title(condition)
    return low_rank, U, S, Vt 
end


function plot_psth(r_all, neural_dex)
    numts = Int(size(r_all,2)/4/4);
    numtsC= Int(size(r_all,2)/4);
    figure()
    cdex = 0;
    plot(r_all[neural_dex,cdex*numtsC+1:cdex*numtsC+numts])
    plot(r_all[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts])
    figure()
    cdex = 1;
    plot(r_all[neural_dex,cdex*numtsC+1:cdex*numtsC+numts])
    plot(r_all[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts])
    figure()
    cdex = 2;
    plot(r_all[neural_dex,cdex*numtsC+1:cdex*numtsC+numts])
    plot(r_all[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts])
    figure()
    cdex = 3;
    plot(r_all[neural_dex,cdex*numtsC+1:cdex*numtsC+numts])
    plot(r_all[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts])

end
# BASIC PIPELINE:
# results  = load_farm_params("C17")
# response = build_response_matrix(results)
# SVD_analysis()


