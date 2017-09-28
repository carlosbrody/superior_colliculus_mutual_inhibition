

using MAT


# include("general_utils.jl")
# include("hessian_utils.jl")
# # include("pro_anti.jl")
# include("marino_pro_anti.jl")
# include("rate_networks.jl")
# include("pro_anti_opto.jl")


include("general_utils.jl")
include("hessian_utils.jl")
include("pro_anti.jl")
include("marino_pro_anti.jl")
include("rate_networks.jl")
include("pro_anti_opto_marino.jl")


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
:sigma  =>  0, 
:input  =>  0, 
:g_leak =>  1, 
:U_rest =>  0, 
:theta  =>  0.05, 
:beta   =>  0.5, 
# :sw     =>  0.2,
# :hw     =>  -1.7,
# :vw     =>  -1.7,
:constant_excitation      => 0, 
:anti_rule_strength       => 0.05,
:pro_rule_strength        => 0.05, 
:target_period_excitation => 0,
:right_light_excitation   => 0.6, 
:right_light_pro_extra    => 0,
:const_add                => 0, 
:init_add                 => 0, 
:rule_and_delay_period    => 0.2,
:target_period            => 0.1,
:post_target_period       => 0,
:const_pro_bias           => 0.0427,
:nPro                     => 100,
:nAnti                    => 100,
:theta1                   => 0.05,
:theta2                   => 0.15,
:opto_strength  => .9,
:rule_and_delay_periods   => [0.4],
:target_periods   => [0.2],
:post_target_periods   => [0],
:opto_periods   => [-1  -1 ; 0   20 ;0    0.2 ; 0.2  0.4;0.4  20],
#:opto_targets   => [.75 .73;.77 .58;.75 .74;.72 .66;.73 .75] 
:opto_targets => [.9 .7; .9 .5; .9 .7; .9 .5; .9 .7]
);




jjj=1

# load existing fit
# fbasename="FarmFields/farm_JJ"*lpad(jjj,4,0)
# loadname = fbasename*".mat"
# loadname=fbasename
# FIT = matread(loadname)
# args = FIT[:"args"];
# pars = FIT[:"pars"];

# pars = [0.697416, -0.667894, 0.133318, 0.89205, -0.00809659, 1.51729, 0.564566, 0.0445603, 0.01, 0.95];
args=[]
pars=[]

# model_params = merge(model_params, make_dict(args, pars));


# @printf("Starting optogenetic optimization")
println(jjj)




# define function to get long form info of opto model, evaluate on test noise
# standard_func =  (;params...) -> JJ_opto_marino(10, 10; rule_and_delay_periods=FIT[:"rule_and_delay_periods"], theta1=FIT[:"theta1"], theta2=FIT[:"theta2"], post_target_periods=FIT[:"post_target_periods"], seedrand=FIT[:"sr"]+1, cbeta=FIT[:"cb"], verbose=false, merge(make_dict(args,pars,model_params), Dict(params))...)
standard_func =  (;params...) -> JJ_opto_marino(10, 10;  merge(make_dict(args,pars,model_params), Dict(params))...)

# run opto model with all outputs, evaluate on test noise
dif1_opto,dif2_opto,runs_pro,runs_anti, opto_fractions, pro_inputs, anti_inputs, weightss = standard_func(;make_dict(args, pars, model_params)...)




# # define function to get long form info of opto model, evaluate on test noise
# standard_func =  (;params...) -> JJ_marino(1000, 1000; rule_and_delay_periods=FIT[:"rule_and_delay_periods"], theta1=FIT[:"theta1"], theta2=FIT[:"theta2"], post_target_periods=FIT[:"post_target_periods"], seedrand=FIT[:"sr"]+1, cbeta=FIT[:"cb"], verbose=false, merge(make_dict(args,pars,model_params), Dict(params))...)


# dif1,dif2 = standard_func(;make_dict(args, pars, model_params)...)





# fbasename="farm_JJ_runs_"*lpad(jjj,4,0)
# savename = fbasename*".mat"



# write file
# matwrite(savename, Dict("args"=>args, "pars"=>pars, "dif1_opto"=>dif1_opto, "dif2_opto"=>dif2_opto ))


# write file
matwrite("model_opto_test_runs.mat", Dict("args"=>args, "pars"=>pars, "runs_pro"=>runs_pro, "runs_anti"=>runs_anti, "opto_fractions"=>opto_fractions, "pro_inputs"=>pro_inputs, "anti_inputs"=>anti_inputs, "weightss"=>weightss ))











