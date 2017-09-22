

using MAT


include("pro_anti.jl")




model_params = Dict(
:dt     =>  0.02,    # timestep, in secs
:tau    =>  0.1,     # tau, in ms
:vW     =>  -1.7,    # vertical weight
:hW     =>  -1.7,    # horizontal weight
:sW     =>  0.2,     # self-connection weight
:dW     =>  0,       # diagonal weight
:nsteps =>  2,       # number of timesteps in the simulation
:noise  =>  [],      # noise added during simulation. Can be empty matrix, or an nunits-by-nsteps matrix
:sigma  =>  0.08,    # standard deviation of Gaussian noise added (will be scaled by sqrt(dt) to be relatively dt-INsensitive)
:input  =>  0,       # input current. Can be scalar, nunits-by-1, or nunits-by-nsteps matrix
:g_leak =>  0.25,    # leak conductance
:U_rest =>  -1,      # resting membrane potential
:theta  =>  1,       # inverse slope of g() function
:beta   =>  1,       # offset to g() function
:constant_excitation      => 0.19,   # constant input, added to all units at all timesteps
:anti_rule_strength       => 0.1,    # input added only to anti units during rule_and_delay_period in Anti trials
:pro_rule_strength        => 0.1,    # input added only to pro units during rule_and_delay_period in Pro trials
:const_pro_bias           => 0,      # input added only to pro units during all times in all trial types
:target_period_excitation => 1,      # input added to all units during target_period
:right_light_excitation   => 0.5,    # input added to the Anti and the Pro unit on one side during the target_period
:right_light_pro_extra    => 0,      # input added to the right side Pro unit alone during the target_period
:rule_and_delay_period    => 0.4,    # duration of rule_and_delay_period, in secs
:target_period            => 0.1,    # duration of target_period, in secs
:post_target_period       => 0.5,    # duration of post_target_period, in secs
:const_add => 0,  # from rate_networks.jl, unused here
:init_add  => 0,  # from rate_networks.jl, unused here 
)

# for jjj in [1:6;]

    # load existing fit
    fbasename = "/Users/mpagan/animals/farm_AB_0001"
    loadname = fbasename*".mat"
    FIT = matread(loadname)
    args = FIT[:"args"];
    # pars = FIT[:"pars"];
    pars = [0.697416, -0.667894, 0.133318, 0.89205, -0.00809659, 1.51729, 0.564566, 0.0445603, 0.175292];


    model_params = merge(model_params, make_dict(args, pars));


    @printf("Starting optogenetic optimization")

    opto_func =  (;params...) -> JJ(1000,1000; rule_and_delay_periods=FIT[:"rule_and_delay_periods"], theta1=FIT[:"theta1"], theta2=FIT[:"theta2"], post_target_periods=FIT[:"post_target_periods"], seedrand=FIT[:"sr"], cbeta=FIT[:"cb"], verbose=false, merge(model_params, Dict(params))...)


    bao=opto_func(;make_dict(args,pars)...)



    matwrite("caco.mat", Dict("cost1"=>bao[1], "cost2"=>bao[2], "dif1"=>bao[3], "dif2"=>bao[4],
    "hit1"=>bao[5], "hit2"=>bao[6]))

# end





