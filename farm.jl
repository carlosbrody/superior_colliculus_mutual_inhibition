
using PyCall
using PyPlot
using ForwardDiff
using DiffBase
using MAT

# pygui(true)

import Base.convert
convert(::Type{Float64}, x::ForwardDiff.Dual) = Float64(x.value)
function convert(::Array{Float64}, x::Array{ForwardDiff.Dual}) 
    y = zeros(size(x)); 
    for i in 1:prod(size(x)) 
        y[i] = convert(Float64, x[i]) 
    end
    return y
end

include("general_utils.jl")
include("hessian_utils.jl")



######################################################
#                                                    #
#         FORWARD AND BACKWARD MODELS                #
#                                                    #
######################################################



"""
o = g(z)    squashing tanh function, running from 0 to 1, is equal to 0.5 when input is 0.
"""
function g(z)
    return 0.5*tanh.(z)+0.5
end
    
"""
z = g^-1(o)    inverse of squashing tanh function, input must be in (0, 1), output is zero when passed 0.5.
"""
function ginverse(z)
    return 0.5*log.(z./(1-z))
end


"""
forwardModel(startU; dt=0.01, tau=0.1, nsteps=100, input=[0.1, 0], noise=[], W=[0 -5;-5 0], 
init_add=0, start_add=0, const_add=0, sigma=0, gleak=1, U_rest=0, 
    do_plot=false, nderivs=0, difforder=0, clearfig=true, fignum=1, dUdt_mag_only=false,
    warn_if_unused_params=false)

Runs a tanh() style-network forwards in time, given its starting point, using simple Euler integration
    tau dU/dt = -U + W*V + I
    V = 0.5*tanh(U)+ 0.5

**PARAMETERS:**

startU     A column vector, nunits-by-1, indicating the values of U at time zero


**OPTIONAL PARAMETERS**

dt      Scalar, timestep size

tau     Scalar, in seconds

gleak   
        dUdt will have a term equal to gleak*(U_rest - U)
U_rest

nsteps  Number of timesteps to run, including time=0.

input   Either an nunits-by-1 vector, in which case inputs to each unit are constant
        across time, or a matrix, nunits-by-nsteps, indicating input for each unit at each timepoint.

W       Weight matrix, nunits-by-nunits

init_add    DEPRECATED: Vector or scalar that gets added to the input current at very first timestep.
            Deprecated because this made it dt-dependent. Replaced by start_add.

start_add   Vector or scalar that gets added, once, to the initial U[:,1], before the integration process begins.

const_add   Scalar that gets added to U after every timestep

sigma       After each timestep, add sigma*sqrt(dt)*randn() to each element of U

do_plot   Default false, if true, plots V of up to the first two dimensions

fignum     Figure number on which to plot

clrearfig  If true, the figure is first cleared, otherwise any plot ois overlaid

nderivs, difforder     Required for making sure function can create its own arrays and 
                       still be differentiated

dUdt_mag_only  If true, returns |dUdt|^2 from the first timestep only, then stops.

warn_if_unused_params     If true, pronts out a warning of some of the passed parameters are not used.

** RETURNS:**

Uend Vend       nunits-by-1 vectors representing the final values of U and V that were found.
U, V            nunits-by-nsteps matrices containing the full trajectories

"""
function forwardModel(startU; dt=0.01, tau=0.1, nsteps=100, input=[], noise=[], W=[0 -5;-5 0], 
    init_add=0, start_add=0, const_add=0, do_plot=false, nderivs=0, difforder=0, clearfig=true, fignum=1,
    dUdt_mag_only=false, sigma=0, g_leak=1, U_rest=0, theta=0, beta=1, 
    warn_if_unused_params=false, other_unused_params...)

    if warn_if_unused_params && length(other_unused_params)>0
        @printf("\n\n=== forwardModel warning, had unused params ")
        for k in keys(Dict(other_unused_params))
            @printf("%s, ", k)
        end
    end
    
    my_input = ForwardDiffZeros(size(input,1), size(input,2), nderivs=nderivs, difforder=difforder)
    for i=1:prod(size(input)); my_input[i] = input[i]; end
    input = my_input;
    
    nunits = length(startU)
    if size(startU,2) > size(startU,1)
        error("startU must be a column vector")
    end
    
    # --- formatting input ---
    if ~(typeof(input)<:Array) || prod(size(input))==1  # was a scalar
        input = input[1]*(1+ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder))
    elseif length(input)==0 # was the empty matrix
        input = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    elseif size(input,2)==1     # was a column vector
        input = input*(1+ForwardDiffZeros(1, nsteps, nderivs=nderivs, difforder=difforder))
    end    
    # --- formatting noise ---
    if ~(typeof(noise)<:Array) || prod(size(noise))==1  # was a scalar
        noise = noise*(1+ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder))
    elseif length(noise)==0 # was the empty matrix
        noise = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    elseif size(noise,2)==1     # was a column vector
        noise = noise*(1+ForwardDiffZeros(1, nsteps, nderivs=nderivs, difforder=difforder))
    end    
    
    U = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    V = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    
    if ~(typeof(W)<:Array); W = [W]; end

    W     = reshape(W, nunits, nunits)
    U     = reshape(U, nunits, nsteps)
    V     = reshape(V, nunits, nsteps)
    input = reshape(input, nunits, nsteps)
    noise = reshape(noise, nunits, nsteps)

    input[:,1] += init_add
    input      += const_add

    #@printf("size(U) is (%d,%d), and size(startU) is (%d,%d) and size(noise) is (%d,%d)", 
    #    size(U,1), size(U,2), size(startU,1), size(startU,2), size(noise,1), size(noise,2))
    # @printf("U[1]=%g, noise[1]=%g\n", startU, noise[1])
    U[:,1] = startU + noise[:,1] + start_add; # @printf("Resulting U=%g\n", U[1])
    V[:,1] = g((U[:,1]-theta)/beta); # @printf("Resulting V=%g\n", V[1])
    
    for i=2:nsteps
        dUdt = g_leak*(U_rest -U[:,i-1]) + W*V[:,i-1] + input[:,i-1]
        if dUdt_mag_only; return sum(dUdt.*dUdt); end;
        # @printf("dUdt=%g\n", dUdt[1])
        # @printf("i=%g\n", i)
        # @printf("noise[2]=%g\n", noise[2])
        U[:,i] = U[:,i-1] + (dt/tau)*dUdt + noise[:,i] + sigma*sqrt(dt)*randn(size(U,1),1)
        # @printf("Resulting U[2]=%g\n", U[2])
        V[:,i] = g((U[:,i]-theta)/beta)
        # @printf("Resulting V[2]=%g\n", V[2])
    end

    if do_plot
        figure(fignum)
        if length(startU)==1
            if clearfig; clf(); end;
            t = (0:nsteps-1)*dt
            plot(t, V[1,:], "b-")
            plot(t[1], V[1,1], "g.")
            plot(t[end], V[1,end], "r.")
            xlabel("t"); ylabel("V1"); ylim([-0.01, 1.01])
        elseif length(startU)>=2
            if clearfig; clf(); end;
            plot(V[1,:], V[2,:], "b-")
            plot(V[1,1], V[2,1], "g.")
            plot(V[1,end], V[2,end], "r.")
            xlabel("V1"); ylabel("V2"); 
            xlim([-0.01, 1.01]); ylim([-0.01, 1.01])
        end
    end

    return U[:,end], V[:,end], U, V
end


"""
backwardsModel(endU; dt=0.01, tau=0.1, nsteps=100, input=[0],noise=[],  W=[0 -5;-5 0], 
    do_plot=false, nderivs=0, difforder=0, clearfig=true, fignum=1, tol=1e-15, start_eta=10)

Runs a tanh() style-network BACKWARDS in time, given its ending point, by making a backwards
guess at each timepoint and then using Hessian minimization to find the backwards vector that correctly
leads to the current timestep value.  Uses forwardModel() . The forwards equations are:

    tau dU/dt = -U + W*V + I
    V = 0.5*tanh(U)+ 0.5

**PARAMETERS:**

endU     A column vector, nunits-by-1, indicating the values of U at time=end


**OPTIONAL PARAMETERS:**

dt      Scalar, timestep size

tau     Scalar, in seconds

nsteps  Number of timesteps to run, including time=0.

input   Either an nunits-by-1 vector, in which case inputs to each unit are constant
        across time, or a matrix, nunits-by-nsteps, indicating input for each unit at each timepoint.

W       Weight matrix, nunits-by-nunits

do_plot   Default false, if true, plots V of up to the first two dimensions

tol       Tolerance in the minimization procedure for finding each backwards timestep. Passed on
          to trust_region_Hessian_minimization()

start_eta   Passed on to trust_region_Hessian_minimization()

fignum     Figure number on which to plot

clrearfig  If true, the figure is first cleared, otherwise any plot ois overlaid

nderivs, difforder     Required for making sure function can create its own arrays and 
                       still be differentiated



** RETURNS:**

Ustart Vstart   nunits-by-1 vectors representing the starting values of U and V that were found.
U, V            nunits-by-nsteps matrices containing the full trajectories
costs           1-by-nsteps vector with the final cost from the minimization procedure for each
                timestep. This is the squared difference between the U[t+1] produced by the U[t] 
                guess and the actual U[t+1]

"""
function backwardsModel(endU; nsteps=100, start_eta=10, tol=1e-15, maxiter=400, 
    do_plot=false, init_add=0, start_add=0, dt=0.01, 
    input=[], noise=[], nderivs=0, difforder=0, clearfig=false, fignum=1, params...)    

    nunits = length(endU)

    # --- formatting input ---
    if ~(typeof(input)<:Array) || prod(size(input))==1  # was a scalar
        input = input[1]*(1+ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder))
    elseif length(input)==0 # was the empty matrix
        input = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    elseif size(input,2)==1     # was a column vector
        input = input*(1+ForwardDiffZeros(1, nsteps, nderivs=nderivs, difforder=difforder))
    end    
    # --- formatting noise ---
    if ~(typeof(noise)<:Array)  # was a scalar
        noise = noise*(1+ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder))
    elseif length(noise)==0 # was the empty matrix
        noise = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    elseif size(noise,2)==1     # was a column vector
        noise = noise*(1+ForwardDiffZeros(1, nsteps, nderivs=nderivs, difforder=difforder))
    end    
    
    function J(U1, U2; nderivs=0, difforder=0, noise=[], inputs=[], pars...)
        U2hat = forwardModel(U1; nsteps=2, noise=noise, input=input, nderivs=nderivs, difforder=difforder, pars...)[1]
        U2hat = U2hat
        DU = U2hat - U2
    
        return sum(DU.*DU)
    end
    
    if length(noise)==0
        noise = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    end

    U = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    U = reshape(U, nunits, nsteps)
    costs = ForwardDiffZeros(nsteps, 1, nderivs=nderivs, difforder=difforder)    
    
    U[:,end] = endU
    for i=(nsteps-1):-1:1
        if i==1
            my_init_add = init_add
            my_start_add = start_add
        else
            my_init_add = 0
            my_start_add = 0
        end
                
        U[:,i], costs[i] = trust_region_Hessian_minimization(U[:,i+1], 
            (x) -> J(x, U[:,i+1]; nderivs=length(endU), difforder=2, 
            input=input[:,i:i+1], noise = noise[:,i:i+1], 
            init_add=my_init_add, start_add=my_start_add, params...); 
            verbose=false, start_eta=start_eta, tol=tol, maxiter=maxiter)
        if i>1; U[:,i] += noise[:,i]; end
    end
    
    
    V = g(U)
    
    if do_plot
        figure(fignum)   
        if typeof(params)<:Array; params = Dict(params); end;
        if haskey(params, :dt);     dt     = params[:dt];     end
        if haskey(params, :nsteps); nsteps = params[:nsteps]; end
        if length(endU)==1
            if clearfig; clf(); end;
            t = (0:nsteps-1)*dt
            plot(t, V[1,:], "m-")
            plot(t[1], V[1,1], "go")
            plot(t[end], V[1,end], "ro")            
            ylim([-0.01, 1.01])
        elseif length(endU)>=2
            if clearfig; clf(); end;            
            plot(V[1,:], V[2,:], "m-")
            plot(V[1,1], V[2,1], "go")
            plot(V[1,end], V[2,end], "ro")
            xlim([-0.01, 1.01]); ylim([-0.01, 1.01])
        end
    end
    
    return U[:,1], V[:,1], U, V, costs
end



######################################################
#                                                    #
#         PLOT_PA                                    #
#                                                    #
######################################################



"""
    plot_PA(t, U, V; fignum=1, clearfig=true, rule_and_delay_period=1, target_period=1, post_target_period=1,
        other_unused_params...)

Helper function for plotting ProAnti results
"""
function plot_PA(t, U, V; fignum=1, clearfig=true, rule_and_delay_period=1, target_period=1, post_target_period=1,
    other_unused_params...)
    figure(fignum)
    if clearfig; clf(); end
    
    ax1 = subplot(3,1,1)
    h = plot(t, V'); 
    setp(h[1], color=[0, 0, 1])
    setp(h[2], color=[1, 0, 0])
    setp(h[3], color=[1, 0.5, 0.5])
    setp(h[4], color=[0, 1, 1])
    ylabel("V")

    ax = gca()
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([V[:];oldlims[1]])-0.1, maximum([V[:];oldlims[2]])+0.1)
    yl = [ylim()[1], ylim()[2]]
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            -0.05, 1.05, linewidth=2)
    if yl[1]<0.02
        yl[1] = -0.02
    end
    if yl[2]>0.98
        yl[2] = 1.02
    end
    ylim(yl)
    grid(true)
    remove_xtick_labels(ax1)
        
    ax2 = subplot(3,1,2)
    hu = plot(t, U')
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([U[:];oldlims[1]])-0.1, maximum([U[:];oldlims[2]])+0.1)
    setp(hu[1], color=[0, 0, 1])
    setp(hu[2], color=[1, 0, 0])
    setp(hu[3], color=[1, 0.5, 0.5])
    setp(hu[4], color=[0, 1, 1])
    ylabel("U"); 
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            ylim()[1], ylim()[2], linewidth=2)
    remove_xtick_labels(ax2)

    grid(true)
    
    subplot(3,1,3)
    delta = V[1,:] - V[4,:]
    hr = plot(t, delta)
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([delta[:];oldlims[1]])-0.1, maximum([delta[:];oldlims[2]])+0.1)
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            ylim()[1], ylim()[2], linewidth=2)
    xlabel("t"); ylabel("Pro R - Pro L")
    grid(true)
        
end


######################################################
#                                                    #
#      MODEL SETUP, MAKE_INPUT(), RUN_TRIALS()       #
#                                                    #
######################################################




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
:const_add => 0, 
:init_add  => 0, 
:rule_and_delay_period    => 0.4,
:target_period            => 0.1,
:post_target_period       => 0.5,
:const_pro_bias           => 0,
)


function make_input(trial_type; dt=0.02, nderivs=0, difforder=0, constant_excitation=0.19, anti_rule_strength=0.1, 
    pro_rule_strength=0.1, target_period_excitation=1, right_light_excitation=0.5, right_light_pro_extra=0, 
    rule_and_delay_period=0.4, target_period=0.1, post_target_period=0.4, const_pro_bias=0,
    other_unused_params...)

    T = rule_and_delay_period + target_period + post_target_period
    t = 0:dt:T
    nsteps = length(t)

    input = constant_excitation + ForwardDiffZeros(4, nsteps, nderivs=nderivs, difforder=difforder)
    if trial_type=="Anti"
        input[2:3, t.<rule_and_delay_period] += anti_rule_strength
    elseif trial_type=="Pro"
        input[[1,4], t.<rule_and_delay_period] += pro_rule_strength
    else
        error("make_input: I don't recognize input type \"" * trial_type * "\"")
    end
    
    input[:,     (rule_and_delay_period.<=t) & (t.<rule_and_delay_period+target_period)] += target_period_excitation
    input[1:2,   (rule_and_delay_period.<=t) & (t.<rule_and_delay_period+target_period)] += right_light_excitation
    input[1,     (rule_and_delay_period.<=t) & (t.<rule_and_delay_period+target_period)] += right_light_pro_extra
    
    input[[1,4],:] += const_pro_bias
    
    return input, t, nsteps
end


function run_ntrials(nPro, nAnti; plot_list=[], nderivs=0, difforder=0, model_params...)
    pro_input,  t, nsteps = make_input("Pro" ; model_params...)
    anti_input, t, nsteps = make_input("Anti"; model_params...)

    model_params = Dict(model_params)
    sW = model_params[:sW]
    hW = model_params[:hW]
    vW = model_params[:vW]
    dW = model_params[:dW]
    model_params = make_dict(["nsteps", "W"], [nsteps, [sW vW dW hW; vW sW hW dW; dW hW sW vW; hW dW vW sW]], 
        model_params)
    model_params = make_dict(["nderivs", "difforder"], [nderivs, difforder], model_params)
    
    proVs  = ForwardDiffZeros(4, nPro, nderivs=nderivs, difforder=difforder)
    antiVs = ForwardDiffZeros(4, nAnti, nderivs=nderivs, difforder=difforder)

    # --- PRO ---
    figure(1); clf();
    model_params = make_dict(["input"], [pro_input], model_params)
    for i=1:nPro
        startU = [-0.3, -0.7, -0.7, -0.3]
        Uend, Vend, U, V = forwardModel(startU, do_plot=false; model_params...)
        proVs[:,i] = Vend
        if any(plot_list.==i) 
            plot_PA(t, U, V; fignum=1, clearfig=false, model_params...)
            subplot(3,1,1); title("PRO")
        end
    end

    # --- ANTI ---
    figure(2); clf();
    model_params = make_dict(["input"], [anti_input], model_params)
    for i=1:nAnti
        startU = [-0.7, -0.3, -0.3, -0.7]
        Uend, Vend, U, V = forwardModel(startU, do_plot=false; model_params...)
        antiVs[:,i] = Vend
        if any(plot_list.==i) 
            plot_PA(t, U, V; fignum=2, clearfig=false, model_params...)
            subplot(3,1,1); title("ANTI")
        end
    end
    
    return proVs, antiVs
end

nPro = 10; nAnti = 5;
proVs, antiVs = @time(run_ntrials(nPro, nAnti; plot_list=[1:5;], model_params...))

@printf("Pro %% correct = %g%%\n", 100*length(find(proVs[1,:].>proVs[4,:]))/nPro)
@printf("Anti %% correct = %g%% \n", 100*length(find(antiVs[1,:].<antiVs[4,:]))/nAnti)



######################################################
#                                                    #
#      JJ() COST FUNCTION                            #
#                                                    #
######################################################


function JJ(nPro, nAnti; pro_target=0.9, anti_target=0.7, 
    theta1=0.025, theta2=0.035, cbeta=0.003, verbose=false, 
    pre_string="", zero_last_sigmas=0, seedrand=NaN, 
    rule_and_delay_periods = [0.4], target_periods = [0.1], post_target_periods = [0.5],
    nderivs=0, difforder=0, model_params...)

    nruns = length(rule_and_delay_periods)*length(target_periods)*length(post_target_periods)
    
    cost1s = ForwardDiffZeros(1, nruns, nderivs=nderivs, difforder=difforder)
    cost2s = ForwardDiffZeros(1, nruns, nderivs=nderivs, difforder=difforder)

    if ~isnan(seedrand); srand(seedrand); end
    
    n = totHitsP = totHitsA = totDiffsP = totDiffsA = 0
    for i in rule_and_delay_periods
        for j in target_periods
            for k = post_target_periods
                n += 1
                
                my_params = make_dict(["rule_and_delay_period", "target_period", "post_target_period"],
                [i, j, k], Dict(model_params))
    
                # print("model params is " ); print(model_params); print("\n")
                proVs, antiVs = run_ntrials(nPro, nAnti; nderivs=nderivs, difforder=difforder, my_params...)

                hitsP  = 0.5*(1 + tanh.((proVs[1,:]-proVs[4,:,])/theta1))
                diffsP = tanh.((proVs[1,:,]-proVs[4,:])/theta2).^2
                hitsA  = 0.5*(1 + tanh.((antiVs[4,:]-antiVs[1,:,])/theta1))
                diffsA = tanh.((antiVs[4,:,]-antiVs[1,:])/theta2).^2

                if nPro>0 && nAnti>0
                    cost1s[n] = (nPro*(mean(hitsP) - pro_target).^2  + nAnti*(mean(hitsA) - anti_target).^2)/(nPro+nAnti)
                    cost2s[n] = -cbeta*(nPro*mean(diffsP) + nAnti*mean(diffsA))/(nPro+nAnti)
                elseif nPro>0
                    cost1s[n] = (mean(hitsP) - pro_target).^2
                    cost2s[n] = -cbeta*mean(diffsP)
                else
                    cost1s[n] = (mean(hitsA) - anti_target).^2
                    cost2s[n] = -cbeta*mean(diffsA)
                end

                totHitsP  += mean(hitsP);  totHitsA  += mean(hitsA); 
                totDiffsP += mean(diffsP); totDiffsA += mean(diffsA);
            end
        end
    end
    
    cost1 = mean(cost1s)
    cost2 = mean(cost2s)

    hitsP = totHitsP/n; hitsA = totHitsA/n; diffsP = totDiffsP/n; diffsA = totDiffsA/n
    
    if verbose
        @printf("%s", pre_string)
        @printf("     -- cost=%g,   cost1=%g, cost2=%g\n", 
            convert(Float64, cost1+cost2), convert(Float64, cost1), convert(Float64, cost2))
        if nPro>0 && nAnti>0
            @printf("     -- mean(hitsP)=%g, mean(diffsP)=%g mean(hitsA)=%g, mean(diffsA)=%g\n", 
                convert(Float64, mean(hitsP)), convert(Float64, mean(diffsP)),
                convert(Float64, mean(hitsA)), convert(Float64, mean(diffsA)))
        elseif nPro>0
            @printf("     -- mean(hitsP)=%g, mean(diffsP)=%g (nAnti=0)\n", 
                convert(Float64, mean(hitsP)), convert(Float64, mean(diffsP)))
        else
            @printf("     -- (nPro=0) mean(hitsA)=%g, mean(diffsA)=%g\n", 
                convert(Float64, mean(hitsA)), convert(Float64, mean(diffsA)))
        end        
    end
    
    return cost1 + cost2
end



JJ(2, 10; plot_list=1:5, verbose=true, model_params...)


######################################################
#                                                    #
#      STANDARDIZING COST AS CB=0.01                 #
#                                                    #
######################################################

"""
cost = standard_cost(args, pars, sr, model_params)

Computes the cost as if cb had been 0.01
"""
function standard_cost(args, pars, sr, model_params)
    cb = 0.01
    nPro=100; nAnti=100

    rule_and_delay_periods = [0.4, 1.2]
    post_target_periods    = [0.5, 1.5]

    func = (;params...) -> JJ(nPro, nAnti; rule_and_delay_periods=rule_and_delay_periods,
            post_target_periods=post_target_periods,
            seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...)
    
    return func(;make_dict(args, pars)...)
end
    
"""
cost = standard_cost(filename)

Returns the standard cost (at cb=0.01) and inserts it into the file with key "scost" if it wasn't there
already
"""
function standard_cost(filename; verbose=false)
    A = matread(filename)
    if !haskey(A, "scost")
        get!(A, "scost", standard_cost(A["args"], A["pars"], A["sr"], symbol_key_ize(A["model_params"])))
        if verbose
            @printf("File %s did not have scost, adding its value %g\n", filename, A["scost"])
        end
        matwrite(filename, A)
    end
    return A["scost"]
end


######################################################
#                                                    #
#          THE ACTION                                #
#                                                    #
######################################################


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

basename = "farm_E_"

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
        @printf("Going with seed = "); print_vector_g(myseed); print("\n")
        pars, traj, cost, cpm_traj = bbox_Hessian_keyword_minimization(myseed, args, bbox, 
            (;params...) -> JJ(nPro, nAnti; rule_and_delay_periods=rule_and_delay_periods,
	    theta1=theta1, theta2=theta2,
            post_target_periods=post_target_periods,
            seedrand=sr, cbeta=cb, verbose=false, merge(model_params, Dict(params))...),
            start_eta = 0.01, tol=1e-12, verbose=true, verbose_every=10, maxiter=400)
        @printf("Came out with cost %g and pars = ", cost); print_vector_g(pars); print("\n\n")

        myfilename = next_file(basename, 4)

        matwrite(myfilename, Dict("args"=>args, "myseed"=>myseed, "pars"=>pars, "traj"=>traj,
        "cost"=>cost, "cpm_traj"=>cpm_traj, "nPro"=>nPro, "nAnti"=>nAnti, "sr"=>sr, "cb"=>cb,
	"theta1"=>theta1, "theta2"=>theta2,
        "model_params"=>ascii_key_ize(model_params), "bbox"=>ascii_key_ize(bbox), "sbox"=>ascii_key_ize(sbox),
        "rule_and_delay_periods"=>rule_and_delay_periods, "post_target_periods"=>post_target_periods))

        standard_cost(myfilename)  # make sure the cost at cb=0.01 is included, for comparison of results across cb values
    end
end
