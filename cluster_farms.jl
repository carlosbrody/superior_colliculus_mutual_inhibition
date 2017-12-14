
using JLD
using MAT


"""
    cluster_farms(farm_id; farmdir="MiniOptimized", opto_conditions = 3, use_reduced_SVD=true)

Runs the Matlab code that computes the optimum number of clusters according to BIC and the indices of k-means clustering

# PARAMETERS
- farm_id           e.g. "C17"

- farmdir           eg "MiniOptimized"

- opto_conditions   either 1 for control only, or 3 for all opto conditions

- use_reduced_SVD   If TRUE, uses truncated SVD response matrix with shortened rule period

# RETURNED VALUES
ndim = dimensionality of response space
bic = BIC value for # of clusters from 2 to 5
ngroups = number of groups with minimum BIC value
idx = indices resulting from k-means clustering



"""

function cluster_farms(farm_id; farmdir="MiniOptimized", opto_conditions = 3, use_reduced_SVD=true)


    # Load responses from all models
    if use_reduced_SVD
        reduced = "_reduced";
    else
        reduced = "";       
    end
    if opto_conditions > 1
        input_filename = farmdir*farm_id*"_SVD_response_matrix"*string(opto_conditions)*reduced;    
    else
        input_filename = farmdir*farm_id*"_SVD_response_matrix"*reduced;
    end
    response, results = load(input_filename*".jld", "response","results");


    # save responses as Matlab file
    matwrite("data_for_matlab_temp.mat", Dict("response"=>response,"results"=>results));


    # launch Matlab function to perform clustering
    run(`matlab -nodisplay -nojvm -nosplash -r "run('cluster_farms.m');exit;" ">>" /dev/null`);

    # read data from Matlab
    d = matread("data_for_julia_temp.mat")

    ndim = d[:"ndim"];
    bic = d[:"bic"];
    ngroups = d[:"ngroups"];
    idx = d[:"idx"];

    # remove temporary files
    run(`rm data_for_matlab_temp.mat`);
    run(`rm data_for_julia_temp.mat`);

    return ndim,bic,ngroups,idx


end


# cluster_farms("C17")
# print(ndim)
# print(ngroups)


