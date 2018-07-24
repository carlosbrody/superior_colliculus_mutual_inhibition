THIS FILE IS ALEX'S ATTEMPT (07/2018) AT HAVING SOME TOP LEVEL ORGANIZATION

IMPORTANT FILES FOR FITTING:
    pro_anti.jl
    rate_networks.jl    
    optimization_utils.jl
    general_utils.jl
    constrained_parabolic_minimization.jl
    gradient_utils.jl

IMPORTANT FILES FOR ANALYSIS:
    results_analysis.jl analysis code for farms, includes some snippets from svd_cluster
    svd_cluster.jl      analysis based on neural dynamics of the proanti model
    rule_encoding.jl    analysis based on neural encoding of task rule
    cluster_farms.jl    clusters farms based on model parameters, calls some matlab code
    cluster_farms.m     code that actually clusters farms

DATA FILES:
    for each <Farm>:
        <Farm> #Directory with each individual farm run
        <Farm>_results.jld
        <Farm>_encoding.jld
        <Farm>_hessians.jld
        <Farm>_SVD_response_matrix.jld
        <Farm>_SVD_response_matrix3.jld
        <Farm>_SVD_response_matrix_reduced.jld
        <Farm>_SVD_response_matrix3_reduced.jld

TO START A NEW FARM OPTIMIZATION:
    include("opto_reduced_farmC30.jl")          Will start a single fitting process
    include("spock_opto_reduced_farmC30.jl")    Starts a single process but takes in random seed number
    sbatch --array=0-9 ./start_farm_spockC30.sh Starts 10 fitting processes with different random seeds, for use on spock computing server
 
 




