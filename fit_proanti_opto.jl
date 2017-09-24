FarmName = "I"

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
:dt     =>  0.02, 
:tau    =>  0.1, 
:vW     =>  -1.7,
:hW     =>  -1.7,
:sW     =>  0.2,
:dW     =>  0,
:nsteps =>  2, 
:noise  =>  [], 
:sigma  =>  0.08, 
:input  =>  0, 
:g_leak =>  0.25, 
:U_rest =>  -1,
:theta  =>  1, 
:beta   =>  1, 
:sw     =>  0.2,
:hw     =>  -1.7,
:vw     =>  -1.7,
:constant_excitation      => 0.19, 
:anti_rule_strength       => 0.1,
:pro_rule_strength        => 0.1, 
:target_period_excitation => 1,
:right_light_excitation   => 0.5, 
:right_light_pro_extra    => 0,
:const_add                => 0, 
:init_add                 => 0, 
:rule_and_delay_period    => 0.4,
:target_period            => 0.1,
:post_target_period       => 0.5,
:const_pro_bias           => 0,
:nPro                     => 100,
:nAnti                    => 100,
:theta1                   => 0.05,
:theta2                   => 0.15,
:opto_strength  => .7,
:opto_periods   => [-1  -1 ; 0   20 ;0    1 ; 1  1.5;1.5  20],
#:opto_targets   => [.75 .73;.77 .58;.75 .74;.72 .66;.73 .75] 
:opto_targets => [.9 .7; .9 .5; .9 .7; .9 .5; .9 .7]
);

# ======= ARGUMENTS AND SEED VALUES:
args = ["sW", "vW", "hW", "dW", "constant_excitation", "right_light_excitation", "target_period_excitation", "const_pro_bias", "sigma"];
seed = [0.2,   1,   0.2,  1,    0.39,                0.15,                       0.1,                        0.1,              0.1];   

# ======= BOUNDING BOX:
bbox = Dict(:sW=>[0 3], :vW=>[-3 3], :hW=>[-3 3], :dW=>[-3 3], :constant_excitation=>[-2 2], :right_light_excitation=>[0.05 4], :target_period_excitation=>[0.05 4], :const_pro_bias=>[-2 2], :sigma=>[0.01 0.2]);

# ======== SEARCH ZONE:
sbox = Dict(:sW=>[0.1 .5], :vW=>[-.5 .5], :hW=>[-.5 .5], :dW=>[-.5 .5], :constant_excitation=>[-.5 .5], :right_light_excitation=>[0.15 .5], :target_period_excitation=>[0.15 .5],:const_pro_bias=>[-.5 .5], :sigma=>[0.02 0.19]);

# define a few hyper parameters
cbetas                  = [0.02];
rule_and_delay_periods  = [1.5];
post_target_periods     = [0.5];
num_eval_runs           = 1000;
num_optimize_iter       = 2000;
num_optimize_restarts   = 100;

# define base filename
fbasename = "FarmFields/farm_"*string(FarmName);


for i=1:num_optimize_restarts   # iterate many optimization runs
for cb in cbetas                # Iterate over beta values
    # figure out initial seed for random number generator
    sr = convert(Int64, round(time()))
    srand(sr);
    #myseed = copy(seed);
    # set initial parameters by sampling from sbox
    myseed = ForwardDiffZeros(length(args),1);
    for j=1:length(args)
        sym = Symbol(args[j])
        if haskey(sbox, sym)
            myseed[j] = sbox[sym][1] + diff(sbox[sym],2)[1]*rand();
        else
            myseed[j] = seed[j];
        end
    end

    # define non-opto model, returning just value
    func =  (;params...) -> JJ(model_params[:nPro], model_params[:nAnti]; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)[1]
    
    # run optimization of non-opto model
    @printf("Going with seed = "); print_vector_g(myseed); print("\n")
    pars, traj, cost, cpm_traj = bbox_Hessian_keyword_minimization(myseed, args, bbox, func, start_eta = 1, tol=1e-12, verbose=true, verbose_every=10, maxiter=num_optimize_iter)
    @printf("Came out with cost %g and pars = ", cost); print_vector_g(pars); print("\n\n")

    # get gradient and hessian at end of optimization 
    value, grad, hess = keyword_vgh(func, args, pars)
    
    # evaluate to get long form info of non-opto model with all outputs
    test_sr = convert(Int64, round(time()))
    srand(test_sr);
    standard_func =  (;params...) -> JJ(num_eval_runs,num_eval_runs; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods,  seedrand=test_sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)     
   
    # run non-opto model to get all outputs
    scost, scost1, scost2, hitsP,hitsA, diffsP, diffsA = standard_func(;make_dict(args, pars, model_params)...)

    # Run second optimization step on optogenetic inactivations
    opto_args=["opto_strength"];
    opto_bbox = Dict(:opto_strength=>[0 1]);
    opto_sbox = Dict(:opto_strength=>[.7 .99]);
   
    # get initial value of opto fraction by sampling from opto_sbox 
    srand(sr);
    opto_seed = ForwardDiffZeros(length(opto_args),1);
    sym = Symbol(opto_args[1])
    opto_seed[1] = opto_sbox[sym][1] + diff(opto_sbox[sym],2)[1]*rand();

    # define opto model with just value output
    @printf("Starting optogenetic optimization")
    opto_func =  (;params...) -> JJ_opto(model_params[:nPro], model_params[:nAnti]; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=sr, cbeta=cb, verbose=false, merge(make_dict(args,pars,model_params), Dict(params))...)[1]

    # run optimization on opto model, over just opto fraction
    opto_pars, opto_traj, opto_cost, opto_cpm_traj = bbox_Hessian_keyword_minimization(opto_seed, opto_args, opto_bbox, opto_func, start_eta = 1, tol=1e-12, verbose=true, verbose_every=10, maxiter=num_optimize_iter)

    # get gradient and hessian at final point of opto model
    opto_value, opto_grad, opto_hess = keyword_vgh(opto_func, opto_args, opto_pars)

    # define function to get long form info of opto model, evaluate on training noise
    t_standard_opto_func =  (;params...) -> JJ_opto(num_eval_runs,num_eval_runs; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=sr, cbeta=cb, verbose=false, merge(make_dict(args,pars,model_params), Dict(params))...)
 
    # run opto-model to get all outputs, evaluate on training noise
    t_opto_scost, t_opto_scost1, t_opto_scost2, t_opto_hitsP,t_opto_hitsA, t_opto_diffsP, t_opto_diffsA,t_opto_bP, t_opto_bA = t_standard_opto_func(;make_dict(opto_args, opto_pars, make_dict(args,pars,model_params))...)

    # define function to get long form info of opto model, evaluate on test noise
    standard_opto_func =  (;params...) -> JJ_opto(num_eval_runs,num_eval_runs; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=test_sr, cbeta=cb, verbose=false, merge(make_dict(args,pars,model_params), Dict(params))...)
 
    # run opto-model to get all outputs, evaluate on test noise
    opto_scost, opto_scost1, opto_scost2, opto_hitsP,opto_hitsA, opto_diffsP, opto_diffsA,opto_bP, opto_bA = standard_opto_func(;make_dict(opto_args, opto_pars, make_dict(args,pars,model_params))...)

    # Save this run out to a file
    # get filename
    myfilename = next_file(fbasename, 4)
    myfilename = myfilename*".mat"
    # write to file
    matwrite(myfilename, Dict("args"=>args, "myseed"=>myseed, "pars"=>pars, "traj"=>traj, "cost"=>cost, "cpm_traj"=>cpm_traj, "nPro"=>model_params[:nPro], "nAnti"=>model_params[:nAnti], "sr"=>sr, "cb"=>cb, "theta1"=>model_params[:theta1], "theta2"=>model_params[:theta2], "scost"=>scost,"scost1"=>scost1,"scost2"=>scost2, "value"=>value, "grad"=>grad, "hess"=>hess, "hitsP"=>hitsP,"hitsA"=>hitsA, "diffsP"=>diffsP, "diffsA"=>diffsA, "model_params"=>ascii_key_ize(model_params), "bbox"=>ascii_key_ize(bbox), "sbox"=>ascii_key_ize(sbox), "rule_and_delay_periods"=>rule_and_delay_periods, "post_target_periods"=>post_target_periods,"opto_args"=>opto_args, "opto_bbox"=>ascii_key_ize(opto_bbox),"opto_sbox"=>ascii_key_ize(opto_sbox),"opto_seed"=>opto_seed,"opto_pars"=>opto_pars, "opto_traj"=>opto_traj, "opto_cost"=>opto_cost, "opto_cpm_traj"=>opto_cpm_traj,"opto_value"=>opto_value, "opto_grad"=>opto_grad, "opto_hess"=>opto_hess, "opto_scost"=>opto_scost, "opto_scost1"=>opto_scost1, "opto_scost2"=>opto_scost2, "opto_hitsP"=>opto_hitsP, "opto_hitsA"=>opto_hitsA, "opto_diffsP"=>opto_diffsP, "opto_diffsA"=>opto_diffsA,"test_sr"=>test_sr,"opto_bP"=>opto_bP, "opto_bA"=>opto_bA, "t_opto_scost"=>t_opto_scost, "t_opto_scost1"=>t_opto_scost1, "t_opto_scost2"=>t_opto_scost2, "t_opto_hitsP"=>t_opto_hitsP, "t_opto_hitsA"=>t_opto_hitsA, "t_opto_diffsP"=>t_opto_diffsP, "t_opto_diffsA"=>t_opto_diffsA,"t_opto_bP"=>t_opto_bP, "t_opto_bA"=>t_opto_bA ))
end
end

