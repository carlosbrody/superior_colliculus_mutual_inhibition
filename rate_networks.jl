# DON'T MODIFY THIS FILE -- the source is in file Rate Networks.ipynb


using PyCall
using PyPlot
using ForwardDiff
using DiffBase
using MAT

# pygui(true)

include("general_utils.jl")
include("constrained_parabolic_minimization.jl")
include("hessian_utils.jl")





# DON'T MODIFY THIS FILE -- the source is in file Rate Networks.ipynb



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


