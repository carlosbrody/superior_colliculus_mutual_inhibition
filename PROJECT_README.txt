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
    parameter_analysis.jl   analysis for examining parameters
    unilateral_analysis.jl  simulate unilateral analysis

DATA FILES:
    for each <Farm>:
        <Farm> #Directory with each individual farm run
        <Farm>_results.jld
        <Farm>_encoding.jld
        <Farm>_hessians.jld
        <Farm>_clusters.jld
        <Farm>_SVD_response_matrix.jld
        <Farm>_SVD_response_matrix3.jld
        <Farm>_SVD_response_matrix_reduced.jld
        <Farm>_SVD_response_matrix3_reduced.jld

TO START A NEW FARM OPTIMIZATION:
    include("opto_reduced_farmC30.jl")          Will start a single fitting process
    include("spock_opto_reduced_farmC30.jl")    Starts a single process but takes in random seed number
    sbatch --array=0-9 ./start_farm_spockC30.sh Starts 10 fitting processes with different random seeds, for use on spock computing server
  
ANALYSIS STEPS:
    1.copy files from spock to local dir

    2.make_mini_farm()        in results_analysis.jl
    > make_mini_farm("C30"; farmdirs ="Farms_C30", todir="MiniC30")

    3. update_farm()           in svd_cluster.jl
    cluster_farms()         in cluster_farms.jl
    histo_params()          in results_analysis.jl
    scatter_by_arg()        in parameter_analysis.jl 



