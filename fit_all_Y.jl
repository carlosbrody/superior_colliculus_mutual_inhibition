FarmName = "Y"

# Incorporate packages
using ForwardDiff
using DiffBase
using MAT

# Include helper functions
include("pro_anti.jl")
include("pro_anti_opto.jl")

# Define core model parameters
model_params = Dict(
:dt     =>  0.002,      # Changed hyper-parameter
:tau    =>  0.02,       # Changed hyper-parameter
:vW     =>  -1.7, #FIT parameter
:hW     =>  -1.7, #FIT parameter
:sW     =>  0.2,  #FIT parameter
:dW     =>  0,    #FIT parameter
:nsteps =>  2, 
:noise  =>  [], 
:sigma  =>  0.08, #FIT parameter
:input  =>  0,    #### Alex thinks this is out-dated 
:g_leak =>  1,          # Changed hyper-parameter
:U_rest =>  0,          # Changed hyper-parameter
:theta  =>  0.05,       # Changed hyper-parameter
:beta   =>  0.5,        # Changed hyper-parameter 
:sw     =>  0.2,  #### Alex thinks this is out-dated
:hw     =>  -1.7, #### Alex thinks this is out-dated
:vw     =>  -1.7, #### Alex thinks this is out-dated
:constant_excitation      => 0.19, #FIT  parameter
:anti_rule_strength       => 0.05,      # Changed hyper-parameter
:pro_rule_strength        => 0.05,      # Changed hyper-parameter
:target_period_excitation => 1,    #FIT  parameter
:right_light_excitation   => 0.5,  #FIT parameter
:right_light_pro_extra    => 0,
:const_add                => 0, 
:init_add                 => 0, 
:rule_and_delay_period    => 0.4,
:target_period            => 0.1,
:post_target_period       => 0.5,
:const_pro_bias           => 0,    #FIT parameter
:nPro                     => 200,
:nAnti                    => 200,
:theta1                   => 0.05,
:theta2                   => 0.15,
:start_pro                => [-0.5, -0.5, -0.5, -0.5],
:start_anti               => [-0.5, -0.5, -0.5, -0.5],
:opto_strength  => .7,             #FIT parameter
:opto_periods   => [-1.1 -1; 0 20;0 50; 50 100; 100 20],  
# The opto "conditions" correspond to the rows of opto_periods.
# all conditions are in seconds relative to start of the trial
# any value before 0 gets changed to 0
# any value after the end of trial gets changed to end of trial
# Special values allow for variable durations
#  20 codes for "end of trial"
#  -1 codes for "start of trial"
# 100 codes for "end of rule and delay period"
# 200 codes for "end of target period"
# -------------------------------------------------------------
# first column is frachit Pro, next column is Anti, rows are conditions
# Actual Opto targets
#:opto_targets   => [.75 .73;.77 .58;.75 .74; .72 .66;.73 .75] 
# Fake Targets
:opto_targets => [.9 .7; .9 .5;.9 .7; .9 .5; .9 .7]  
);

# ======= ARGUMENTS AND SEED VALUES:
args = ["sW", "vW", "hW", "dW", "constant_excitation", "right_light_excitation", "target_period_excitation", "const_pro_bias", "sigma","opto_strength"];
seed = [0.001,  -1.58,   -0.05,  0.001,    0.001,                0.6,                       0.001,                        0.0427,              0.05, .9];   
# ======= BOUNDING BOX:
bbox = Dict(:sW=>[0 3], :vW=>[-3 3], :hW=>[-3 3], :dW=>[-3 3], :constant_excitation=>[-2 2],
:right_light_excitation=>[0.05 4], :target_period_excitation=>[0 4], :const_pro_bias=>[-2 2],
:sigma=>[0.01 2],:opto_strength=>[0 1]);

# ======== SEARCH ZONE:
sbox = Dict(:sW=>[0 .5], :vW=>[-.5 .5], :hW=>[-.5 .5], :dW=>[-.5 .5],
:constant_excitation=>[-.5 .5], :right_light_excitation=>[0.05 .5], :target_period_excitation=>[0.001 .5],:const_pro_bias=>[-.5 .5], :sigma=>[0.02 1],:opto_strength=>[.7 .99]);

# define a few hyper parameters
cbetas = [0.04];
rule_and_delay_periods = [0.2];
post_target_periods    = [0.05];
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

    dista=rand();
    println("dista")
    println(dista)
    myseed1 = copy(seed);

    # get initial parameter values by sampling from sbox
    myseed = ForwardDiffZeros(length(args),1);
    for j=1:length(args)
        sym = Symbol(args[j])
        if haskey(sbox, sym)
#           myseed[j] = sbox[sym][1] + diff(sbox[sym],2)[1]*rand();
            myseed[j] = dista*(sbox[sym][1] + diff(sbox[sym],2)[1]*rand())+(1-dista)*myseed1[j];
	    if myseed[j] > bbox[sym][2]
	    	myseed[j] =bbox[sym][2] - 0.1;
	    elseif myseed[j] < bbox[sym][1]
	        myseed[j] =bbox[sym][1]+0.1;
	    end
        else
            myseed[j] = seed[j];
        end
    end

    # define opto function with just value output
    func =  (;params...) -> JJ_opto(model_params[:nPro], model_params[:nAnti]; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)[1]
    
    # run optimization with all parameters
    @printf("Going with seed = "); print_vector_g(myseed); print("\n")
    pars, traj, cost, cpm_traj = bbox_Hessian_keyword_minimization(myseed, args, bbox, func, start_eta = 1, tol=1e-12, verbose=true, verbose_every=10, maxiter=num_optimize_iter)
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

   ## Save this run out to a file
    # get filename
    myfilename = next_file(fbasename, 4)
    myfilename = myfilename*".mat"
    # write file
    matwrite(myfilename, Dict("args"=>args, "myseed"=>myseed, "pars"=>pars, "traj"=>traj, "cost"=>cost, "cpm_traj"=>cpm_traj, "nPro"=>model_params[:nPro], "nAnti"=>model_params[:nAnti], "sr"=>sr, "cb"=>cb, "theta1"=>model_params[:theta1], "theta2"=>model_params[:theta2],"value"=>value,"grad"=>grad, "hess"=>hess, "model_params"=>ascii_key_ize(model_params), "bbox"=>ascii_key_ize(bbox), "sbox"=>ascii_key_ize(sbox), "rule_and_delay_periods"=>rule_and_delay_periods, "post_target_periods"=>post_target_periods, "opto_scost"=>opto_scost, "opto_scost1"=>opto_scost1, "opto_scost2"=>opto_scost2, "opto_hitsP"=>opto_hitsP, "opto_hitsA"=>opto_hitsA, "opto_diffsP"=>opto_diffsP, "opto_diffsA"=>opto_diffsA,"test_sr"=>test_sr,"opto_bP"=>opto_bP, "opto_bA"=>opto_bA, "t_opto_scost"=>t_opto_scost, "t_opto_scost1"=>t_opto_scost1, "t_opto_scost2"=>t_opto_scost2, "t_opto_hitsP"=>t_opto_hitsP, "t_opto_hitsA"=>t_opto_hitsA, "t_opto_diffsP"=>t_opto_diffsP, "t_opto_diffsA"=>t_opto_diffsA,"t_opto_bP"=>t_opto_bP, "t_opto_bA"=>t_opto_bA,"dista"=>dista  ))
end
end

