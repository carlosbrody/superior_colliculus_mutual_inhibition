
FarmName = "G"
# Pkg.add("ForwardDiff")
# Pkg.add("MAT")
# Pkg.add("PyPlot")

# Incorporate packages
#using PyCall
#using PyPlot
using ForwardDiff
using DiffBase
using MAT

# Needed for ForwardDiff, apparently
import Base.convert
convert(::Type{Float64}, x::ForwardDiff.Dual) = Float64(x.value)
function convert(::Array{Float64}, x::Array{ForwardDiff.Dual}) 
    y = zeros(size(x)); 
    for i in 1:prod(size(x)) 
        y[i] = convert(Float64, x[i]) 
    end
    return y
end

# Include helper functions
include("../general_utils.jl")
include("../constrained_parabolic_minimization.jl")
include("../hessian_utils.jl")
include("pro_anti.jl")
include("rate_networks.jl")

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
);

# ======= ARGUMENTS AND SEED VALUES:
args = ["sW", "vW", "hW", "dW", "constant_excitation", "right_light_excitation", "target_period_excitation", "const_pro_bias", "sigma"];
seed = [0.2,   1,   0.2,  1,    0.39,                0.15,                       0.1,                        0.1,              0.1];   

# ======= BOUNDING BOX:
bbox = Dict(:sW=>[0 3], :vW=>[-3 3], :hW=>[-3 3], :dW=>[-3 3], :constant_excitation=>[-2 2],
:right_light_excitation=>[0.05 4], :target_period_excitation=>[0.05 4], :const_pro_bias=>[-2 2],
:sigma=>[0.01 0.2]);

# ======== SEARCH ZONE:
sbox = Dict(:sW=>[0.1 2.9], :vW=>[-2.9 2.9], :hW=>[-2.9 2.9], :dW=>[-2.9 2.9],
:constant_excitation=>[-1.9 1.9], :right_light_excitation=>[0.15 3.9], :target_period_excitation=>[0.15 3.9],
:const_pro_bias=>[-1.9 1.9], :sigma=>[0.02 0.19]);


# define a few hyper parameters
cbetas = [0.02, 0.04, 0.1];
rule_and_delay_periods = [0.4, 1.2];
post_target_periods    = [0.5, 1.5];

fbasename = "FarmFields/farm_"*string(FarmName);

# figure out initial seed
sr = convert(Int64, round(time()))
srand(sr);
#myseed = copy(seed);
myseed = ForwardDiffZeros(length(args),1);
for i=1:length(args)
    sym = Symbol(args[i])
    if haskey(sbox, sym)
        myseed[i] = sbox[sym][1] + diff(sbox[sym],2)[1]*rand();
    else
        myseed[i] = seed[i];
    end
end

# define data to evaluate against

# Iterate over beta values, using same initial seed?
for i=1:10
for cb in cbetas
    func =  (;params...) -> JJ(model_params[:nPro], model_params[:nAnti]; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)[1]
    
    # And at the standard cb=0.01 for comparison to other cbs
    standard_func =  (;params...) -> JJ(model_params[:nPro], model_params[:nAnti]; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods,  seedrand=sr, cbeta=0.01, verbose=false, merge(model_params, Dict(params))...)[1]        

    # run optimization                
    @printf("Going with seed = "); print_vector_g(myseed); print("\n")
    pars, traj, cost, cpm_traj = bbox_Hessian_keyword_minimization(myseed, args, bbox, func,
        start_eta = 0.01, tol=1e-12, verbose=true, verbose_every=10, maxiter=400)
    @printf("Came out with cost %g and pars = ", cost); print_vector_g(pars); print("\n\n")

    # get gradient and hessian at end of optimization (couldn't we get this out of the optimizer???)
    value, grad, hess = keyword_vgh(func, args, pars)
    
    # evaluate at standard beta just for comparison
#    scost = standard_func(;make_dict(args, pars, model_params)...)

    standard_func_split =  (;params...) -> JJ_split(model_params[:nPro], model_params[:nAnti]; rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=post_target_periods, seedrand=sr, cbeta=0.01, verbose=false, merge(model_params, Dict(params))...)        

    scost1, scost2 = standard_func_split(;make_dict(args, pars, model_params)...)
    scost = scost1 + scost2;

    # Run second optimization step on optogenetic inactivations
    
    # Save this run out to a file
    myfilename = next_file(fbasename, 4)
    myfilename = myfilename*".mat"
    matwrite(myfilename, Dict("args"=>args, "myseed"=>myseed, "pars"=>pars, "traj"=>traj,
    "cost"=>cost, "cpm_traj"=>cpm_traj, "nPro"=>model_params[:nPro], "nAnti"=>model_params[:nAnti], "sr"=>sr, "cb"=>cb,
    "theta1"=>model_params[:theta1], "theta2"=>model_params[:theta2],
    "scost"=>scost,"scost1"=>scost1,"scost2"=>scost2, "value"=>value, "grad"=>grad, "hess"=>hess,
    "model_params"=>ascii_key_ize(model_params), "bbox"=>ascii_key_ize(bbox), "sbox"=>ascii_key_ize(sbox),
    "rule_and_delay_periods"=>rule_and_delay_periods, "post_target_periods"=>post_target_periods))
end
end

