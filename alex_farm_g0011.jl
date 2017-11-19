# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb. Look there for further documentation and examples of running the code.


include("pro_anti.jl")

A = matread("FarmFields/farm_G_1.mat0011")

model_params = symbol_key_ize(A["model_params"])
for k in keys(A); print(k); print("  "); end

A


# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb. Look there for further documentation and examples of running the code.


bbox = Dict(:sW=>[0 3], :vW=>[-3 3], :hW=>[-3 3], :dW=>[-3 3], :constant_excitation=>[-2 2],
:right_light_excitation=>[0.05 4], :target_period_excitation=>[0.05 4], :const_pro_bias=>[-2 2],
:sigma=>[0.01 0.2]);

model_params = symbol_key_ize(A["model_params"])

    rule_and_delay_periods = [0.4, 1.2]
    post_target_periods    = [0.5, 1.5]

    standard_func =  (;params...) -> JJ(model_params[:nPro], model_params[:nAnti]; 
    rule_and_delay_periods=rule_and_delay_periods, theta1=model_params[:theta1], theta2=model_params[:theta2], 
    post_target_periods=post_target_periods,  seedrand=A["sr"], cbeta=0.01, verbose=true, 
    merge(model_params, Dict(params))...)[1]        

value, grad, hess = keyword_vgh(standard_func, A["args"], A["myseed"])



# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb. Look there for further documentation and examples of running the code.


pars, traj, cost, cpm_traj = bbox_Hessian_keyword_minimization(A["myseed"], A["args"], bbox, standard_func,
        start_eta = 0.03, tol=1e-18, verbose=true, verbose_every=1, maxiter=4000)



