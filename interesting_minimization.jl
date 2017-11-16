# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


include("pro_anti.jl")


# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


pygui(true)
d = load("FarmFields/farm_C3_0005.jld")
mypars     = d["mypars"]
extra_pars = d["extra_pars"] 

seed = d["seed"]
args = d["args"]
bbox = d["bbox"]

extra_pars[:opto_conditions] = []
extra_pars[:plot_list] = [1:10;]
extra_pars[:plot_condition] = 1

# This function will get the output of new_J() at each iteration, and will return "true", stopping
# the minimization, if hBP and hBA are both above a threshold.
function stopping_func(;cost=0, func_out=[], ignored_extra_params...)
    costf, dP, dA, hBP, hBA = func_out
    return hBP[1]>=0.6 && hBA[1] >= 0.6
end

func =  (;params...) -> new_J(10, 10; verbose=true, merge(merge(mypars, extra_pars), Dict(params))...)

t_pars, traj, cost, cpm_traj, ftraj = bbox_Hessian_keyword_minimization(seed, args, bbox, func, 
    stopping_function = stopping_func, 
start_eta = 0.1, tol=1e-12, verbose=true, verbose_every=1, maxiter=2000)



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


seed = t_pars

func =  (;params...) -> JJ(mypars[:nPro], mypars[:nAnti]; verbose=false, 
    merge(merge(mypars, extra_pars), Dict(params))...)


pars, traj, cost, cpm_traj, ftraj = bbox_Hessian_keyword_minimization(seed, args, bbox, func, 
    start_eta = 0.1, tol=1e-12, verbose=true, verbose_every=1, maxiter=1000)



