
using JLD
using MAT


"""
    cluster_farms(farm_id; farmdir="MiniOptimized", cost_choice="cost", threshold=-0.0002,
                       nclusters=collect(2:7), opto_conditions = 3, use_reduced_SVD=true, further_params...)

Runs the Matlab code that computes the optimum number of clusters according to BIC and the indices of k-means clustering

# PARAMETERS
- farm_id           e.g. "C17"

- farmdir           eg "MiniOptimized"

- opto_conditions   either 1 for control only, or 3 for all opto conditions

- use_reduced_SVD   If TRUE, uses truncated SVD response matrix with shortened rule period

- cost_choice       One of "cost" or "tcost", for test cost versus training cost, respectively

- threshold         Only use runs with cost less than this

- nclusters         A vector of the number of clusters to test BIC for. Note default is collect(2:7) as
                    opposed to just 2:7 because the latter is a Range type that can't be used with matwrite, whereas
                    the former is a vector that works fine.
        
- further_params    Any further params will be passed directly on to Matlab. 

# RETURNED VALUES
ndim = dimensionality of response space
bic = BIC value for # of clusters from 2 to 5
ngroups = number of groups with minimum BIC value
idx = indices resulting from k-means clustering



"""

function cluster_farms(farm_id; farmdir="MiniOptimized", cost_choice="cost", threshold=-0.0002,
                       nclusters=collect(2:7), opto_conditions = 3, use_reduced_SVD=true, further_params...)


    # Load responses from all models
    if use_reduced_SVD
        reduced = "_reduced";
    else
        reduced = "";       
    end
    if opto_conditions > 1
        input_filename = farmdir*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*reduced;    
    else
        input_filename = farmdir*"_"*farm_id*"_SVD_response_matrix"*reduced;
    end
    response, results = load(input_filename*".jld", "response","results");

    dict_for_matlab = merge(Dict(further_params), 
          Dict("response"=>response, "results"=>results, "nclusters"=>nclusters,
               "cost_choice"=>cost_choice, "threshold"=>threshold))

    try
        # save responses as Matlab file
        matwrite("compute_clustering/data_for_matlab_temp.mat", dict_for_matlab);
    catch y
        error("Couldn't write matlab file, are you giving me parameters that matwrite() doesn't",
              "understand, e.g., Range objects like 2:7?\n", y)
        catch_stacktrace()
    end

    # launch Matlab function to perform clustering
    run(`matlab -nodisplay -nojvm -nosplash -r "run('cluster_farms.m');exit;" ">>" /dev/null`);

    # read data from Matlab
    d = matread("compute_clustering/data_for_julia_temp.mat")

    ndim    = d["ndim"];
    bic     = d["bic"];
    ngroups = d["ngroups"];
    idx     = d["idx"];
    files   = d["files"]

    # remove temporary files
    run(`rm compute_clustering/data_for_matlab_temp.mat`);
    run(`rm compute_clustering/data_for_julia_temp.mat`);

    save("compute_clustering/MarinoCode_"*farm_id*"_"*farmdir*".jld", 
         merge(dict_for_matlab, Dict("ndim"=>ndim, "bic"=>bic, "ngroups"=>ngroups, "idx"=>idx, "files"=>files)))

    return ndim,bic,ngroups,idx


end


# cluster_farms("C17")
# print(ndim)
# print(ngroups)


