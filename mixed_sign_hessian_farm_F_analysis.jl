# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb. Look there for further documentation and examples of running the code.


"""
walled, free = id_pars_at_walls(bbox, args, pars; tol=1e-3)

Returns a list of elements in pars that are within tol of their bounds described in the bbox Dict(),
and a list of those not within tol of their bounds.  The union of walled and free will be equal to 1:length(pars).

"""
function id_pars_at_walls(bbox, args, pars; tol=1e-3)
    walled = Array{Int}(0)
    free   = Array{Int}(0)
    i = 1; while i<=length(args)
        if typeof(args[i])<:String && haskey(bbox, Symbol(args[i]))
            range = bbox[Symbol(args[i])]
            if any(abs(pars[i]-range).<tol)
                walled = [walled; i]
            else
                free = [free; i]
            end    
        elseif typeof(args[i])<:String
            free = [free; i]
        elseif haskey(bbox, Symbol(args[i][1]))
            range = bbox[Symbol(args[i][1])]
            for j=1:args[i][2]
                if any(abs(pars[i+j-1]-range).<tol)
                    walled = [walled; i+j-1]
                else
                    free = [free; i+j-1]
                end    
            end
            i = i+args[i][2]-1
        else
            for j=1:args[i][2]; free = [free; i+j-1]; end
            i = i+args[i][2]-1
        end
    i=i+1; end    
    
    return walled, free
end


walled, free = id_pars_at_walls(Dict(:a=>[1, 2], :b=>[0.1, 0.2]), ["c", "a", ["b" 2]], [10, 1.9999, 0.1, 0.15])

walled, free = id_pars_at_walls(bbox, args, pars; tol=1e-1)


# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb. Look there for further documentation and examples of running the code.


# ========  This is copied from the code for farm F, telling us how it was run ==============================

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

# ========  END  --   This is copied from the code for farm F, telling us how it was run ===================

# Now read all the files to get all there costs and iteration numbers

fnames = filter(x -> startswith(x, "farm_F"), readdir("FarmFields/"))
scosts = zeros(size(fnames))
niters = zeros(size(fnames))
for i=1:length(fnames) 
    A = matread("FarmFields/" * fnames[i])  
    scosts[i] = A["scost"]
    niters[i] = size(A["traj"],2)
end

pygui(true)
figure(1); clf();
plot(scosts, niters, ".")
xlabel("cost at standard beta=0.01")
ylabel("number of iterations run")


# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb. Look there for further documentation and examples of running the code.



# #################  FINALLY TO THE ANALYSIS PART   ################

goodies = find((niters.<400)) #  & (scosts.<-0.007))
mixed_sign = Array{String}(0)

for i=1:length(goodies)
    A = matread("FarmFields/" * fnames[goodies[i]])

    args = A["args"]
    pars = A["pars"]
    nPro = A["nPro"]
    nAnti = A["nAnti"]
    rule_and_delay_periods = A["rule_and_delay_periods"]
    post_target_periods    = A["post_target_periods"]
    theta1 = A["theta1"]
    theta2 = A["theta2"]
    sr = A["sr"]
    cb = A["cb"]
    sbox = symbol_key_ize(A["sbox"])
    bbox = symbol_key_ize(A["bbox"])
    model_params = symbol_key_ize(A["model_params"])


    func =  (;params...) -> JJ(nPro, nAnti; rule_and_delay_periods=rule_and_delay_periods,
        theta1=theta1, theta2=theta2,
        post_target_periods=post_target_periods,
        seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)[1]

    standard_func =  (;params...) -> JJ(nPro, nAnti; rule_and_delay_periods=rule_and_delay_periods,
        theta1=theta1, theta2=theta2,
        post_target_periods=post_target_periods,
        seedrand=sr, cbeta=0.01, verbose=false, merge(model_params, Dict(params))...)[1]


    value, grad, hess = keyword_vgh(func, args, pars)

    walled, free = id_pars_at_walls(bbox, args, pars)
    L, V = eig(hess[free,free])
    @printf("%d/%d: file %s had %d free parameters, with hessian eigenvalues ", 
        i, length(goodies), fnames[goodies[i]], length(free))
    print_vector_g(L); print("\n")
    
    if any(L.<0)
        @printf("\n    *** File %s had off-wall mixed sign eigenvalues!\n\n", fnames[goodies[i]])
        mixed_sign = [mixed_sign ; fnames[goodies[i]]]
    end
end




