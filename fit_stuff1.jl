FarmName = "MA"

# Incorporate packages
using ForwardDiff
using DiffBase
using MAT

# Include helper functions
include("general_utils.jl")
include("constrained_parabolic_minimization.jl")
include("hessian_utils.jl")
include("pro_anti.jl")
include("rate_networks.jl")
include("pro_anti_opto.jl")

# Define core model parameters
model_params = Dict(
:dt     =>  0.002, 
:tau    =>  0.02, 
:vW     =>  -1.58,
:hW     =>  -0.05,
:sW     =>  0,
:dW     =>  0,
:nsteps =>  301, 
:noise  =>  [], 
:sigma  =>  0.08, 
:input  =>  0, 
:g_leak =>  1, 
:U_rest =>  0,
:theta  =>  0.05, 
:beta   =>  0.5, 
:constant_excitation      => 0, 
:anti_rule_strength       => 0.05,
:pro_rule_strength        => 0.05, 
:target_period_excitation => 0,
:right_light_excitation   => 0.6, 
:right_light_pro_extra    => 0,
:const_add                => 0, 
:init_add                 => 0, 
:rule_and_delay_period    => 0.4,
:target_period            => 0.2,
:post_target_period       => 0.1,
:const_pro_bias           => 0.0427,
:nPro                     => 100,
:nAnti                    => 100,
:theta1                   => 0.05,
:theta2                   => 0.15,
:opto_strength  => .9,
:opto_periods   => [-1  -1 ; 0   20 ;0    0.2 ; 0.2  0.4;0.4  20],  # set of opto conditions, in seconds, with 0 the start 
# of the trial (i.e. start of rule_and_delay_period), anything before 0 or after end of trial gets ignored.
# :opto_targets   => [.75 .73;.77 .58;.75 .74;.72 .66;.73 .75] 
:opto_targets => [.9 .7; .9 .5; .9 .7; .9 .5; .9 .7]  # first column is frachit Pro, next column is Anti, rows are conditions
# The "conditions" correspond to the rows of opto_periods.
);

# ======= ARGUMENTS AND SEED VALUES:
args = ["vW", "hW", "dW",  "right_light_excitation", "const_pro_bias", "sigma","opto_strength"];
seed = [-1.58,   -0.05,    0.001,     0.6,                 0.0427,          0.05, .9];   

# ======= BOUNDING BOX:
bbox = Dict(:vW=>[-3 3], :hW=>[-3 3], :dW=>[-3 3], 
:right_light_excitation=>[0.05 4], :const_pro_bias=>[-2 2],
:sigma=>[0.01 0.2],:opto_strength=>[0 1]);

# ======== SEARCH ZONE:
sbox = Dict(:vW=>[-.5 .5], :hW=>[-.5 .5], :dW=>[-.5 .5],
:right_light_excitation=>[0.15 .5], :const_pro_bias=>[-.5 .5], :sigma=>[0.02 0.19],:opto_strength=>[.7 .99]);

# define a few hyper parameters
cbetas = [0.04];
rule_and_delay_periods = [0.4];
post_target_periods    = [0.1];
num_eval_runs           = 1000;
num_optimize_iter       = 2000;
num_optimize_restarts   = 100;

# define base filename
fbasename = "FarmFields/farm_"*string(FarmName);

for i=1:num_optimize_restarts   # Iterate this many optimizations
for cb in cbetas                # Iterate over beta values, if there are multiple
    # figure out initial seed for random number generator
    sr = convert(Int64, round(time()))
    srand(sr);
    myseed = copy(seed);
    # get initial parameter values by sampling from sbox
    #myseed = ForwardDiffZeros(length(args),1);
    #for j=1:length(args)
    #    sym = Symbol(args[j])
    #    if haskey(sbox, sym)
    #        myseed[j] = sbox[sym][1] + diff(sbox[sym],2)[1]*rand();
    #    else
    #        myseed[j] = seed[j];
    #    end
    #end

    # define opto function with just value output
    func =  (;params...) -> JJ_opto(model_params[:nPro], model_params[:nAnti]; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)[1]
    
    # run optimization with all parameters
    @printf("Going with seed = "); print_vector_g(myseed); print("\n")
    pars, traj, cost, cpm_traj = bbox_Hessian_keyword_minimization(myseed, args, bbox, func, start_eta = 1, tol=1e-12, verbose=true, verbose_every=2, maxiter=num_optimize_iter)
    @printf("Came out with cost %g and pars = ", cost); print_vector_g(pars); print("\n\n")

    # get gradient and hessian at end of optimization 
    value, grad, hess = keyword_vgh(func, args, pars)

    # define function with all outputs, evaluate on training noise
    t_standard_func =  (;params...) -> JJ_opto(model_params[:nPro], model_params[:nAnti]; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)

    # run opto model with all outputs, evaluate on training noise
    t_opto_scost, t_opto_scost1, t_opto_scost2, t_opto_hitsP,t_opto_hitsA, t_opto_diffsP, t_opto_diffsA, t_opto_bP, t_opto_bA = t_standard_func(;make_dict(args, pars, model_params)...)
    
   ## evaluate for long form info
    # reset random number generator for testing purposes
    test_sr = convert(Int64, round(time()))
    srand(test_sr); 

    # define function with all outputs, evaluate on test noise
    standard_func =  (;params...) -> JJ_opto(num_eval_runs, num_eval_runs; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=test_sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)

    # run opto model with all outputs, evaluate on test noise
    opto_scost, opto_scost1, opto_scost2, opto_hitsP,opto_hitsA, opto_diffsP, opto_diffsA, opto_bP, opto_bA = standard_func(;make_dict(args, pars, model_params)...)

    # define non-opto model with all outputs, to check opto model, evaluate on test noise
    standard_func =  (;params...) -> JJ(num_eval_runs, num_eval_runs; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=test_sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)

    # run non-opto model with all outputs, evaluate on test noise
    scost, scost1, scost2, hitsP,hitsA, diffsP, diffsA = standard_func(;make_dict(args, pars, model_params)...)
 
   ## Save this run out to a file
    # get filename
    myfilename = next_file(fbasename, 4)
    myfilename = myfilename*".mat"
    # write file
    matwrite(myfilename, Dict("args"=>args, "myseed"=>myseed, "pars"=>pars, "traj"=>traj, "cost"=>cost, "cpm_traj"=>cpm_traj, "nPro"=>model_params[:nPro], "nAnti"=>model_params[:nAnti], "sr"=>sr, "cb"=>cb, "theta1"=>model_params[:theta1], "theta2"=>model_params[:theta2],"value"=>value,"grad"=>grad, "hess"=>hess, "scost"=>scost, "scost1"=>scost1, "scost2"=>scost2, "hitsP"=>hitsP,"hitsA"=>hitsA, "diffsP"=>diffsP, "diffsA"=>diffsA, "model_params"=>ascii_key_ize(model_params), "bbox"=>ascii_key_ize(bbox), "sbox"=>ascii_key_ize(sbox), "rule_and_delay_periods"=>rule_and_delay_periods, "post_target_periods"=>post_target_periods, "opto_scost"=>opto_scost, "opto_scost1"=>opto_scost1, "opto_scost2"=>opto_scost2, "opto_hitsP"=>opto_hitsP, "opto_hitsA"=>opto_hitsA, "opto_diffsP"=>opto_diffsP, "opto_diffsA"=>opto_diffsA,"test_sr"=>test_sr,"opto_bP"=>opto_bP, "opto_bA"=>opto_bA, "t_opto_scost"=>t_opto_scost, "t_opto_scost1"=>t_opto_scost1, "t_opto_scost2"=>t_opto_scost2, "t_opto_hitsP"=>t_opto_hitsP, "t_opto_hitsA"=>t_opto_hitsA, "t_opto_diffsP"=>t_opto_diffsP, "t_opto_diffsA"=>t_opto_diffsA,"t_opto_bP"=>t_opto_bP, "t_opto_bA"=>t_opto_bA  ))
end
end

