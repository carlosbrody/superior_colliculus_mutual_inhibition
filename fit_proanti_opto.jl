FarmName = "H"
FarmNum = 222;

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

#function optimize_proanti_opto(FarmName, FarmNum)

# load existing fit
fbasename = "FarmFields/farm_"*string(FarmName)*lpad(FarmNum,4,0);
loadname = fbasename*".mat"
FIT = matread(loadname)
args = FIT[:"args"];
pars = FIT[:"pars"];

# If needed include opto parameters
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
:opto_targets   => [.9 .7;.9 .5;.9 .7;.9 .5;.9 .7] 
);
model_params = merge(model_params, make_dict(args, pars));
rule_and_delay_periods = FIT[:"rule_and_delay_periods"];
post_target_periods    = FIT[:"post_target_periods"];
cb = FIT[:"cb"];

# Define which parameters to optimize
opto_args=["opto_strength"];
opto_bbox = Dict(:opto_strength=>[0 1]);
opto_sbox = Dict(:opto_strength=>[.1 .9]);

# figure out initial seed
temp_seed = [.5];
sr = convert(Int64, round(time()))
srand(sr);
opto_seed = ForwardDiffZeros(length(opto_args),1);
for i=1:length(opto_args)
    sym = Symbol(opto_args[i])
    if haskey(opto_sbox, sym)
        opto_seed[i] = opto_sbox[sym][1] + diff(opto_sbox[sym],2)[1]*rand();
    else
        opto_seed[i] = temp_seed[i];
    end
end

# Run second optimization step on optogenetic inactivations
@printf("Starting optogenetic optimization")
opto_func =  (;params...) -> JJ_opto(model_params[:nPro], model_params[:nAnti]; opto_targets=model_params[:opto_targets],opto_periods=model_params[:opto_periods],rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)[1]

# actually do the run
opto_pars, opto_traj, opto_cost, opto_cpm_traj = bbox_Hessian_keyword_minimization(opto_seed, opto_args, opto_bbox, opto_func, start_eta = 1, tol=1e-12, verbose=true, verbose_every=10, maxiter=2000)
@printf("Came out with opto_cost %g and opto_pars = ", cost); print_vector_g(pars); print("\n\n")  

opto_value, opto_grad, opto_hess = keyword_vgh(opto_func, opto_args, opto_pars)

# Save this run out to a file
myfilename = fbasename*"_opto.mat"

matwrite(myfilename, Dict("args"=>args, "myseed"=>myseed, "pars"=>pars, "traj"=>traj, "cost"=>cost, "cpm_traj"=>cpm_traj, "nPro"=>model_params[:nPro], "nAnti"=>model_params[:nAnti], "sr"=>sr, "cb"=>cb, "theta1"=>model_params[:theta1], "theta2"=>model_params[:theta2], "scost"=>scost,"scost1"=>scost1,"scost2"=>scost2, "value"=>value, "grad"=>grad, "hess"=>hess, "model_params"=>ascii_key_ize(model_params), "bbox"=>ascii_key_ize(bbox), "sbox"=>ascii_key_ize(sbox), "rule_and_delay_periods"=>rule_and_delay_periods, "post_target_periods"=>post_target_periods,"opto_args"=>opto_args, "opto_bbox"=>opto_bbox,"opto_sbox"=>opto_sbox,"opto_seed"=>opto_seed,"opto_pars"=>opto_pars, "opto_traj"=>opto_traj, "opto_cost"=>opto_cost, "opto_cpm_traj"=>opto_cpm_traj))


