module RateNetworks
# include("optimization_utils.jl")

using PyCall
if !occursin("spock", String(chomp(read(`hostname`, String))))
    using PyPlot
end
using ForwardDiff
using DiffResults
using MAT
using GradientUtils

# pygui(true)


export forwardModel


"""
forwardModel(startU; dt=0.01, tau=0.1, nsteps=100, W=[0 -5;-5 0],
    sigma=0, g_leak=1, U_rest=0, theta=0, beta=1,
    input=[], noise=[], init_add=0, start_add=0, const_add=0,
    do_plot=false, clearfig=true, fignum=1,
    opto_strength=1, opto_units=[], opto_times=zeros(0,2),
    Iperturb_strength=0, Iperturb_units=[], Iperturb_times=zeros(0,2),
    dUdt_mag_only=false,
    warn_if_unused_params=false, other_unused_params...)

Runs a tanh() style-network forwards in time, given its starting point, using simple Euler integration
    tau dU/dt = -U + W*V + I
    V = 0.5*tanh(U)+ 0.5

# PARAMETERS:

- startU     A column vector, nunits-by-1, indicating the values of U at time zero


# OPTIONAL PARAMETERS

- dt      Scalar, timestep size

- tau     Scalar, in seconds

- gleak   dUdt will have a term equal to gleak*(U_rest - U)

- U_rest  dUdt will have a term equal to gleak*(U_rest - U)

- nsteps  Number of timesteps to run, including time=0.

- input   Either an nunits-by-1 vector, in which case inputs to each unit are constant
        across time, or a matrix, nunits-by-nsteps, indicating input for each unit at each timepoint.

- noise   Meant as a frozen noise input. Same format as input.

- W       Weight matrix, nunits-by-nunits

- init_add    DEPRECATED: Vector or scalar that gets added to the input current at very first timestep.
            Deprecated because this made it dt-dependent. Replaced by start_add.

- start_add   Vector or scalar that gets added, once, to the initial U[:,1], before the integration process begins.

- const_add   Scalar that gets added to U after every timestep

- sigma       After each timestep, add sigma*sqrt(dt)*randn() to each element of U

- opto_strength    The outputs V, after being computed, will get multiplied by this number. opto_strength should *EITHER* be a scalar, in which case optional params opto_units and opto_times below are also relevant; *OR* it should be an nunits-by-nsteps matrix, completely specifying how much each unit's V should be multiplied by at each timestep, in which case opto_times and opto_units are irrelevant

- opto_units       A list of the unit numbers that will have their V multiplied by opto_strength. For example, [1,3] would affect only units 1 and 3.  Can be the empty matrix (equivalent to no opto effect). Irrelevant if opto_strength = 1

- opto_times    An n-by-2 matrix, where each row lists t_start_of_opto_effect, t_end_of_opto_effect. For example,
                [1 3 ; 6 8]  would mean "have an opto effect during both 1 <= t <=3 and 6 <= t <= 8]. With the
                code as currently configured, this would mean the same opto_strength and opto_units across all
                the relevant time intervals in a run.

_ Iperturb_strength   Magnitude of a current that will be added to specified units at specified times. Default is zero

- Iperturb_units      A list of the unit numbers that will have Iperturb_strength added to them.
                      For example, [1,3] would affect only units 1 and 3.
                      Can be the empty matrix (equivalent to no opto effect).
                      Irrelevant if Iperturb_strength = 1. Default is 1:4

- Iperturb_times    An n-by-2 matrix, where each row lists t_start_of_Iperturb_effect, t_end_of_Iperturb_effect.
                For example, [1 3 ; 6 8]  would mean "have an Iperturb effect during both 1 <= t <=3 and
                6 <= t <= 8]. With the code as currently configured, this would mean the same
                Iperturb_strength and Iperturb_units across all the relevant time intervals in a run.

- do_plot   Default false, if true, plots V of up to the first two dimensions

- fignum     Figure number on which to plot

- clearfig  If true, the figure is first cleared, otherwise any plot ois overlaid

- dUdt_mag_only  If true, returns |dUdt|^2 from the first timestep only, then stops.

- warn_if_unused_params     If true, pronts out a warning of some of the passed parameters are not used.



** RETURNS:**

- Uend Vend       nunits-by-1 vectors representing the final values of U and V that were found.

- U, V            nunits-by-nsteps matrices containing the full trajectories

- t               A time vector, so one could things like plot(t, U[1,:])

```jldoctest
    forwardModel([1.1]; do_plot=true, Iperturb_strength=-5, Iperturb_units=1,
        Iperturb_times=[0.04, 0.06],  dt=0.001, W=[-2], nsteps=100, start_add=-1.9, noise=0)
```

"""
function forwardModel(startU; dt=0.01, tau=0.1, nsteps=100, W=[0 -5;-5 0],
    sigma=0, g_leak=1, U_rest=0, theta=0, beta=1,
    input=[], noise=[], init_add=0, start_add=0, const_add=0,
    do_plot=false, clearfig=true, fignum=1,
    opto_strength=1, opto_units=[], opto_times=zeros(0,2),
    Iperturb_strength=0, Iperturb_units=[], Iperturb_times=zeros(0,2),
    dUdt_mag_only=false,
    warn_if_unused_params=false, other_unused_params...)

    @assert size(W,1)==size(W,2)       "W must be square"
    @assert length(startU)==size(W,1)  "length(start(U)) must match size(W)"

    # All the variables that we MIGHT choose to differentiate w.r.t. go into this bag -- further down
    # we'll use get_eltype(varbag) to check for any of them being ForwardDiff.Dual.
    # That is how we'll tell whether new matrices should be regular numbers of ForwardDiff.Dual's.
    # *** if you add a new variable you'll want to differentiate w.r.t., it should be added here too ***
    varbag = (startU, opto_strength, opto_times, dt, tau, input, noise, W,
        init_add, start_add, const_add, sigma, g_leak, U_rest, theta,
        beta, Iperturb_strength, Iperturb_times)

    """
    o = g(z)    squashing tanh function, running from 0 to 1, is equal to 0.5 when input is 0.
    """
    function g(z)
        return 0.5*tanh.(z).+0.5
    end

    if warn_if_unused_params && length(other_unused_params)>0
        println("\n\n=== forwardModel warning, had unused params ")
        println(keys(Dict(other_unused_params)))
    end

    if length(size(opto_times))==1
        opto_times = reshape(opto_times, 1, 2)
    end
    if length(size(Iperturb_times))==1
        Iperturb_times = reshape(Iperturb_times, 1, 2)
    end

    nunits = length(startU)
    if size(startU,2) > size(startU,1)
        error("startU must be a column vector")
    end

    # We copy the input so as not to mess with the original -- remember, Julia passes arrays by reference, not by value
    # -------------------------------------------------------------------------------------------
    #
    #    Formatting and declaring arrays for ForwardDiff version >= 0.6  (Julia 0.6 and onwards)
    #
    # -------------------------------------------------------------------------------------------

    my_input = zeros(get_eltype(varbag), size(input,1), size(input,2))
    for i=1:prod(size(input)); my_input[i] = input[i]; end
    input = my_input;

    # --- formatting input ---
    if ~(typeof(input)<:Array) || prod(size(input))==1  # was a scalar
        input = input[1]*ones(get_eltype(varbag), nunits, nsteps)
    elseif length(input)==0 # was the empty matrix
        input = zeros(get_eltype(varbag), nunits, nsteps)
    elseif size(input,2)==1     # was a column vector
        input = input*ones(get_eltype(varbag), 1, nsteps)
    end
    # --- formatting noise ---
    if ~(typeof(noise)<:Array) || prod(size(noise))==1  # was a scalar
        noise = noise[1]*ones(get_eltype(varbag), nunits, nsteps)
    elseif length(noise)==0 # was the empty matrix
        noise = zeros(get_eltype(varbag), nunits, nsteps)
    elseif size(noise,2)==1     # was a column vector
        noise = noise*ones(get_eltype(varbag), 1, nsteps)
    end
    # --- formatting opto fraction ---
    if typeof(opto_strength)<:Array
        if size(opto_strength,1) != nunits || size(opto_strength,2) != nsteps
            error("opto_strength must be either a scalar or an nunits-by-nsteps matrix")
        end
        opto_matrix = copy(opto_strength)
    else # We assume that if opto_strength is not an Array, then it is a scalar
        opto_matrix = ones(get_eltype(varbag), nunits, nsteps)
        time_axis = dt*(0:nsteps-1)
        for i=1:size(opto_times,1)
            opto_matrix[opto_units, (opto_times[i,1] .<= time_axis) .&
                (time_axis .<= opto_times[i,2])] .= opto_strength
        end
    end
    # --- formatting Iperturb  ---
    if typeof(Iperturb_strength)<:Array
        if size(Iperturb_strength,1) != nunits || size(Iperturb_strength,2) != nsteps
            error("Iperturb_strength must be either a scalar or an nunits-by-nsteps matrix")
        end
        Iperturb_matrix = copy(Iperturb_strength)
    else # We assume that if opto_strength is not an Array, then it is a scalar
        Iperturb_matrix = zeros(get_eltype(varbag), nunits, nsteps)
        time_axis = dt*(0:nsteps-1)
        for i=1:size(Iperturb_times,1)
            Iperturb_matrix[Iperturb_units, (Iperturb_times[i,1] .<= time_axis) .& (time_axis .<= Iperturb_times[i,2])] = Iperturb_strength
        end
    end

    U = zeros(get_eltype(varbag), nunits, nsteps)
    V = zeros(get_eltype(varbag), nunits, nsteps)


    if ~(typeof(W)<:Array); W = [W]; end

    W     = reshape(W, nunits, nunits)
    U     = reshape(U, nunits, nsteps)
    V     = reshape(V, nunits, nsteps)
    input = reshape(input, nunits, nsteps)
    noise = reshape(noise, nunits, nsteps)

    input[:,1] .+= init_add
    input      .+= const_add

    #@printf("size(U) is (%d,%d), and size(startU) is (%d,%d) and size(noise) is (%d,%d)",
    #    size(U,1), size(U,2), size(startU,1), size(startU,2), size(noise,1), size(noise,2))
    # @printf("U[1]=%g, noise[1]=%g\n", startU, noise[1])
    U[:,1] = startU .+ noise[:,1] .+ start_add; # @printf("Resulting U=%g\n", U[1])
    V[:,1] = g((U[:,1] .- theta)/beta);
#    @printf("U[1U[1,1])
    V[:,1] .*= opto_matrix[:,1]

    for i=2:nsteps
        dUdt = g_leak*(U_rest .- U[:,i-1]) .+ W*V[:,i-1] .+ input[:,i-1] .+ Iperturb_matrix[:,i-1]
        if dUdt_mag_only; return sum(dUdt.*dUdt); end;
        # @printf("dUdt=%g\n", dUdt[1])
        # @printf("i=%g\n", i)
        # @printf("noise[2]=%g\n", noise[2])
        U[:,i] = U[:,i-1] .+ (dt/tau)*dUdt .+ noise[:,i] .+ sigma*sqrt(dt)*randn(size(U,1),1)
        # @printf("Resulting U[2]=%g\n", U[2])
        V[:,i] = g((U[:,i] .- theta)/beta)
        V[:,i] .*= opto_matrix[:,i]
        # @printf("Resulting V[2]=%g\n", V[2])
    end

    if do_plot
        figure(fignum)
        myV = get_value(V);
        if length(startU)==1
            if clearfig; clf(); end;
            t = (0:nsteps-1)*dt
            plot(t, myV[1,:], "b-")
            plot(t[1], myV[1,1], "g.")
            plot(t[end], myV[1,end], "r.")
            xlabel("t"); ylabel("V1"); ylim([-0.01, 1.01])
        elseif length(startU)>=2
            if clearfig; clf(); end;
            plot(myV[1,:], myV[2,:], "b-")
            plot(myV[1,1], myV[2,1], "g.")
            plot(myV[1,end], myV[2,end], "r.")
            xlabel("V1"); ylabel("V2");
            xlim([-0.01, 1.01]); ylim([-0.01, 1.01])
        end
    end

    return U[:,end], V[:,end], U, V, (0:nsteps-1)*dt
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

    if FDversion()>=0.6
        # All the variables that we MIGHT choose to differentiate w.r.t. go into this bag -- further down
        # we'll use get_eltype(varbag) to check for any of them being ForwardDiff.Dual.
        # That is how we'll tell whether new matrices should be regular numbers of ForwardDiff.Dual's.
        # *** if you add a new variable you'll want to differentiate w.r.t., it should be added here too ***
        varbag = (endU, dt, init_add, start_add)
    end

    """
    o = g(z)    squashing tanh function, running from 0 to 1, is equal to 0.5 when input is 0.
    """
    function g(z)
        return 0.5*tanh.(z)+0.5
    end

    nunits = length(endU)

    if FDversion() < 0.6
        # -----------------------------------------------------------------------------------
        #
        #    Formatting and declaring arrays for ForwardDiff version < 0.6  (Julia 0.5.2)
        #
        # -----------------------------------------------------------------------------------

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

        if length(noise)==0
            noise = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
        end

        U     = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
        costs = ForwardDiffZeros(nsteps, 1, nderivs=nderivs, difforder=difforder)

    else
        # -------------------------------------------------------------------------------------------
        #
        #    Formatting and declaring arrays for ForwardDiff version >= 0.6  (Julia 0.6 and onwards)
        #
        # -------------------------------------------------------------------------------------------

        # --- formatting input ---
        if ~(typeof(input)<:Array) || prod(size(input))==1  # was a scalar
            input = input*ones(get_eltype(varbag), nunits, nsteps)
        elseif length(input)==0 # was the empty matrix
            input = zeros(get_eltype(varbag), nunits, nsteps)
        elseif size(input,2)==1     # was a column vector
            input = input*ones(get_eltype(varbag), 1, nsteps)
        end
        # --- formatting noise ---
        if ~(typeof(noise)<:Array) || prod(size(noise))==1  # was a scalar
            noise = noise*ones(get_eltype(varbag), nunits, nsteps)
        elseif length(noise)==0 # was the empty matrix
            noise = zeros(get_eltype(varbag), nunits, nsteps)
        elseif size(noise,2)==1     # was a column vector
            noise = noise*ones(get_eltype(varbag), 1, nsteps)
        end

        if length(noise)==0
            noise = zeros(get_eltype(varbag), nunits, nsteps)
        end

        U     = zeros(get_eltype(varbag), nunits, nsteps)
        costs = zeros(get_eltype(varbag), nsteps, 1)

    end

    function J(U1, U2; nderivs=0, difforder=0, noise=[], inputs=[], pars...)
        U2hat = forwardModel(U1; nsteps=2, noise=noise, input=input, nderivs=nderivs, difforder=difforder, pars...)[1]
        U2hat = U2hat
        DU = U2hat - U2

        return sum(DU.*DU)
    end


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


    V = g(U)  # REALLY???? HOW ABOUT THETA AND BETA?

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


end  #  END MODULE
