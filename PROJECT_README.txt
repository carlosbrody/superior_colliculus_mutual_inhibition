This file is Alex's attempt (07/2018) at having some top level organization

Important files for fitting:
    pro_anti.jl
    rate_networks.jl    
    optimization_utils.jl
    general_utils.jl
    constrained_parabolic_minimization.jl
    gradient_utils.jl

Important files for analysis:
    results_analysis.jl analysis code for farms, includes some snippets from svd_cluster
    svd_cluster.jl      analysis based on neural dynamics of the proanti model
    rule_encoding.jl    analysis based on neural encoding of task rule
    cluster_farms.jl    clusters farms based on model parameters, calls some matlab code
    cluster_farms.m     code that actually clusters farms

Data files:
    for each <Farm>:
        <Farm> #Directory with each individual farm run
        <Farm>_results.jld
        <Farm>_encoding.jld
        <Farm>_hessians.jld
        <Farm>_SVD_response_matrix.jld
        <Farm>_SVD_response_matrix3.jld
        <Farm>_SVD_response_matrix_reduced.jld
        <Farm>_SVD_response_matrix3_reduced.jld

To start a new farm optimization:
    include("opto_reduced_farmC30.jl")
 
        




