# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


# ======= ARGUMENTS AND SEED VALUES:
args = ["sW", "vW", "hW", "dW", "constant_excitation", "right_light_excitation", "target_period_excitation"]
seed = [0.2,   1,   0.2,  1,    0.39,                0.15,                       0.1]

args = [args ; ["const_pro_bias", "sigma"]]
seed = [seed ; [0.1,               0.1]]


# ======= BOUNDING BOX:
bbox = Dict(:sW=>[0 3], :vW=>[-3 3], :hW=>[-3 3], :dW=>[-3 3], :constant_excitation=>[-2 2],
:right_light_excitation=>[0.05 4], :target_period_excitation=>[0.05 4], :const_pro_bias=>[-2 2],
:sigma=>[0.01 0.2])

model_params = merge(model_params, Dict(:post_target_period=>0.5))
# seed = [0.0840597,  -1.32677,  -0.437334,  -0.324835,  0.567997, 0.712216,  0.0500075,  0.0858569,  0.25]


# ======== SEARCH ZONE:

sbox = Dict(:sW=>[0.001 0.5], :vW=>[-0.5 0.5], :hW=>[-0.5 0.5], :dW=>[-0.5 0.5],
:constant_excitation=>[-0.5 0.5], :right_light_excitation=>[0.1 0.5], :target_period_excitation=>[0.1 0.5],
:const_pro_bias=>[0 0.2], :sigma=>[0.02 0.19])

cbetas = [0.02, 0.04]

fbasename = "FarmFields/farm_F_"

while true
    myseed = seed;
    sr = convert(Int64, round(time()))

    myseed = copy(seed);
    for i=1:length(args)
        sym = Symbol(args[i])
        if haskey(sbox, sym)
            myseed[i] = sbox[sym][1] + diff(sbox[sym],2)[1]*rand()
        end
    end
    nPro=100; nAnti=100

    rule_and_delay_periods = [0.4, 1.2]
    post_target_periods    = [0.5, 1.5]

    theta1 = 0.15; theta2 = 0.25

    for cb in cbetas
        func =  (;params...) -> JJ(nPro, nAnti; rule_and_delay_periods=rule_and_delay_periods,
            theta1=theta1, theta2=theta2,
            post_target_periods=post_target_periods,
            seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)[1]
        
        # And at the standard cb=0.01 for comparison to other cbs
        standard_func =  (;params...) -> JJ(nPro, nAnti; rule_and_delay_periods=rule_and_delay_periods,
            theta1=theta1, theta2=theta2,
            post_target_periods=post_target_periods,
            seedrand=sr, cbeta=0.01, verbose=false, merge(model_params, Dict(params))...)[1]        
                
        @printf("Going with seed = "); print_vector_g(myseed); print("\n")
        pars, traj, cost, cpm_traj = bbox_Hessian_keyword_minimization(myseed, args, bbox, func,
            start_eta = 0.01, tol=1e-12, verbose=true, verbose_every=10, maxiter=400)
        @printf("Came out with cost %g and pars = ", cost); print_vector_g(pars); print("\n\n")

        value, grad, hess = keyword_vgh(func, args, pars)
        scost = standard_func(;make_dict(args, pars, model_params)...)
        
        myfilename = next_file(fbasename, 4)

        matwrite(myfilename, Dict("args"=>args, "myseed"=>myseed, "pars"=>pars, "traj"=>traj,
        "cost"=>cost, "cpm_traj"=>cpm_traj, "nPro"=>nPro, "nAnti"=>nAnti, "sr"=>sr, "cb"=>cb,
        "theta1"=>theta1, "theta2"=>theta2,
        "scost"=>scost, "value"=>value, "grad"=>grad, "hess"=>hess,
        "model_params"=>ascii_key_ize(model_params), "bbox"=>ascii_key_ize(bbox), "sbox"=>ascii_key_ize(sbox),
        "rule_and_delay_periods"=>rule_and_delay_periods, "post_target_periods"=>post_target_periods))

    end
end



