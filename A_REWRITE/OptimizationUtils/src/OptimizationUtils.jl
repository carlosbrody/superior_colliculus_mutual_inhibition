module OptimizationUtils

using Printf
if !occursin("spock", String(chomp(read(`hostname`, String))))
    using PyPlot
    PyPlot.rc("font", family="Helvetica Neue")
end
using JLD
using LinearAlgebra
using MAT
using Dates


using GeneralUtils
using GradientUtils
using ConstrainedParabolicMinimization


export bbox_Hessian_keyword_minimization, inverse_wall



"""
pdict = wallwrap(bdict, pdict)
Given bdict, a dictionary of symbols to [minval, maxval] vectors, and pdict, a dictionary of symbols
to values (or, alternatively, an Array of (Symbol, value) tuples], goes through each of the symbols in
bdict and modifies the corresponding value in pdict putting it through a tanh so the final output lies
within the limits in bdict.  Returns the new pdict.  Makes a copy of pdict so as not to modify the original.
"""
function wallwrap(bdict, epdict)
    local pdict = deepcopy(epdict)  # Must be very careful here! I got bit by the bug of forgetting that without
    # an explicit copy() call, Julia does not make copies of the contents of arrays or dictionaries, making it
    # easy to inadvertently modify something one did not intend to perturb.  Note the deepcopy() call,
    # necessary to make sure we don't mess up the content of the caller's dictionary.

    pdict = Dict(pdict)  # Make sure it's a Dict, even if it came in as a Named Tuple

    allkeys = keys(bdict)

    for k in allkeys
        local bbox = bdict[k]
        d = 0.5*(bbox[2] - bbox[1])
        m = 0.5*(bbox[2] + bbox[1])

        pdict[k] = bbox[1] + d*(tanh((pdict[k]-m)/d)+1)
    end
    return pdict
end


"""
params = vector_wrap(bbox, args, eparams)
Given bdict, a dictionary of symbols to [minval, maxval] vectors, args, an array of strings representing
symbols, and params, an array of values corresponding to the args list, puts each param that has an entry
in bdict through the tanh-walling mechanism, and returns the result. Does not modify the contents of the
original params vector (or bdict or args).
"""
function vector_wrap(bbox, args, eparams)
    local params = deepcopy(eparams)
    pdict = wallwrap(bbox, make_dict(args, params))
    i=1; j=1
    for i=1:length(args)
        if typeof(args[i])<:Array
            params[j:j+args[i][2]-1] = pdict[Symbol(args[i][1])]
            j += args[i][2]-1
        else
            params[j] = pdict[Symbol(args[i])]
        end
    j = j+1
    end
    return params
end


"""
params = inverse_wall(bdict, args, wparams)
Given bdict, a dictionary of symbols to [minval, maxval] vectors, args, an array of strings representing
symbols, and wparams, an array of values corresponding to the args list where each param that has an entry
in bdict has alreadt been through the tanh-walling mechanism, UNwalls the ones that have a bdict entry and
returns the result. Does not modify the contents of the original params vector (or bdict or args).
"""
function inverse_wall(bdict, args, wparams)
    local params = deepcopy(wparams)
    pdict = inverse_wall(bdict, make_dict(args, params))
    i=1; j=1
    for i=1:length(args)
        if typeof(args[i])<:Array
            params[j:j+args[i][2]-1] = pdict[Symbol(args[i][1])]
            j += args[i][2]-1
        else
            params[j] = pdict[Symbol(args[i])]
        end
        j = j+1
    end
    return params
end


"""
pdict = inverse_wall(bdict, wdict)
Given bdict, a dictionary of symbols to [minval, maxval] vectors, and wdict, a dictionary of symbols to values
(or vectors of values)  UNwalls the ones that have a bdict entry and
returns the result. Does not modify the contents of any dictionaries.
"""
function inverse_wall(bdict, wdict)
    local pdict = deepcopy(wdict)

    allkeys = keys(bdict)
    for k in allkeys
        local bbox = bdict[k]
        d = 0.5*(bbox[2] - bbox[1])
        m = 0.5*(bbox[2] + bbox[1])

        pdict[k] = m + d*0.5*log((pdict[k]-bbox[1])./(2*d - pdict[k] + bbox[1]))
    end
    return(pdict)
end



# DON'T MODIFY THIS FILE -- the source is in file Optimization Utilities.ipynb. Look there for further documentation and examples of running the code.



"""
function adaptive_gradient_minimization(seed, func; start_eta=0.1, tol=1e-6, maxiter=400, verbose=false)

NEEDS DOCMENTING
"""

function adaptive_gradient_minimization(seed, func; start_eta=0.1, tol=1e-6, maxiter=400,
    verbose=false)

    params = seed
    eta = start_eta

    out = DiffBase.GradientResult(params)
    ForwardDiff.gradient!(out, func, params)
    cost = DiffBase.value(out)
    grad = DiffBase.gradient(out)

    for i in [1:maxiter;]
        new_params = params - eta*grad

        ForwardDiff.gradient!(out, func, new_params)
        new_cost = DiffBase.value(out)
        new_grad = DiffBase.gradient(out)

        if abs(new_cost - cost) < tol
            break
        end

        if new_cost >= cost
            eta = eta/2
        else
            eta = eta*1.1
            params = new_params
            cost = new_cost
            grad = new_grad
        end

        if verbose
            @printf "%d: eta=%.3f cost=%.4f ps=[" i eta cost
            for p in [1:length(params);]
                @printf "%.3f" params[p]
                if p<length(params) @printf ", "; end
            end
            @printf "]\n"
        end
    end

    return params
end


#############################################################################
#                                                                           #
#                   TRUST_REGION_HESSIAN_MINIMIZATION                       #
#                                                                           #
#############################################################################



"""
function trust_region_Hessian_minimization(seed, func; start_eta=10, tol=1e-6, maxiter=400,
    verbose=false)

(below, x stands for delta_x, the step from the current x=x0 position at which the cost = const)

cost = 0.5*x'*H*x + grad*x + const

dcost/dx = H*x + grad  ;   dcost/dx = 0  ==> x =  - inv(H)*grad

Trust-region says have a parameter lambda, and replace H with hat{H} = H +  I/eta.
When eta is very large, this is equivalent to a straight Newton method jump,
because hat{H} ~= H.  But when eta is small, this is more like a small gradient
descent step, because for small eta inv(hat{H}) ~= eta and therefore the delta x is like
-eta*grad.  So, if the cost function is going down, make eta larger, and if it is going
up, make eta a lot smaller. Just like we do in other adaptive methods

PARAMETERS:
===========

seed        column vector, representing the starting value of the parameters.

func        Function that takes a vector and returns a scalar.  If you want to
            work with a function that tales more parameterrs and returns more than one
            output, you can use something like

                    x -> orig_func(x, other_params)[1]

            You only need the "[1]" part if the orig_func returns more outputs than a scalar.

OPTIONAL PARAMETERS:
====================

start_eta=10    Starting value of eta.  It's good to start with somethibg biggish, if it is
                too much, it'll quickly get cut down.

tol=1e-15       Numerical tolerance. If a proposed jump produces a change in func that is less than
               this, the minimization stops.

maxiter=400    Maximum number of iterations to do before stopping

verbose=false   If true, print out a report on each iteration of iteration number, radius size (eta),
                what type jump was proposed ("Newton" means going straight to global min, "constrained" means jump has
                norm eta, failed means that finding the minimum at a given radius somehow didn't work). Will also
                print out the cosine of the angle between the proposed jump and the gradient.

RETURNS:
========

params       A vector the size of seed that has the last values of the minimizing parameters for func

"""
function trust_region_Hessian_minimization(seed, func; start_eta=10, tol=1e-15, maxiter=400,
    verbose=false, verbose_level=1)

    params = seed
    eta = start_eta

    cost, grad, hess = vgh(func, params)
    if verbose && verbose_level >= 2
        @printf("Initial cost, grad, hess:\n")
        println(cost)
        println(grad)
        println(hess)
    end


    for i in [1:maxiter;]
        hathess    = hess + eye(length(grad), length(grad))/eta
        new_params = params - inv(hathess)*grad
        new_cost, new_grad, new_hess = vgh(func, new_params)

        if abs(new_cost - cost) < tol
            break
        end

        if new_cost >= cost
            eta = eta/2
            costheta = NaN
        else
            eta = eta*1.1
            costheta = dot(new_params-params, grad)/(norm(new_params-params)*norm(grad))

            params = new_params
            cost = new_cost
            grad = new_grad
            hess = new_hess
        end

        if verbose
            @printf "%d: eta=%.3f cost=%.4f costheta=%.3f ps=" i eta cost  costheta
            print(params)
            @printf "\n"
        end
    end

    return params, cost
end







######################################################
#                                                    #
#         BBOX_HESSIAN_KEYWORD_MINIMIZATION          #
#                                                    #
######################################################




"""
function bbox_Hessian_keyword_minimization(seed, args, bbox, func; start_eta=0.1,
    tol=1e-6, maxiter=400, frac_cost_threshold = 0.5, stopping_function = nothing,
    verbose=false, verbose_level=1, verbose_every=1, verbose_file=stdout,
    verbose_timestamp = false, start_iter_num = 1,
    softbox=true, hardbox=false, report_file="")

Like constrained_Hessian_minimization, but uses keyword_hessian!().

# PARAMETERS:

- seed        column vector, representing the starting value of the parameters.

- args        List of strings identifying parameters for differentiation, e.g., ["const_E", "w_self]

- bbox        If softbox=true (the default), should then be a Dict of Symbol=>[minval maxval] entries. An entry
            in this Dict indicates that the corresponding parameter is to be bounded, as indicated by the associated
            [minval maxval] vector. The bbox dictionary can have fewer entries than the number of parameters, and its
            default value is Dict(), indicating an unbounded search.
                If softbox=false, then bbox should be an nargs-by-2 matrix indicating the range for each argument,
            with the minima (first column) and maxima (second column), and entries for ALL parameters.

- func      func must take only optional keyword args, and must
            take nderivs=0, difforder=0  and declare any new matrices using ForwardDiffZeros() instead of zeros().
            [THE ABOVE PART OF THE DOCUMENTATION MUST BE UPDATED FOR JULIA 0.6]
            func must either return a scalar, or the first output it returns must be a scalar.
            That scalar is what will be minimized. The trajectory across the minimization of
            any further outputs that f() returns will be available in ftraj (see RETURNS below)


# OPTIONAL PARAMETERS:

- start_eta    Starting value of the radius.  It's good to start with somethibg biggish, if it is
             too much, it'll quickly get cut down.

- tol=1e-6     Numerical tolerance. If a proposed jump produces a change in func that is less than
             this, the minimization stops.

- maxiter=400  Maximum number of iterations to do before stopping

- start_iter_num=1  Number to start iterations count at

- frac_cost_threshold   When the algorithm is going to take a step, it first checks whether this will reduce the cost.
                If the answer is "no" then the step is not taken, and the step size is halved.
                Small step sizes will make the algorithm gradient-descent-like, so if the step size keeps
                getting smaller, eventualluy we get to simple gradient descent. If the cost will be reduced,
                another check is done, comparing to the expected cost change, from the quadratic approximation.
                If the ratio of the actual cost reduction / expected cost reduction is less than "frac_cost_threshold",
                then the step is not taken and the step size is halved. Otherwise, the  step is taken, and step
                size goes up by 1.1.  A frac_cost_threshold=0 makes both of these checks equivalent.

- stopping_function   If present, this should be a function that returns a boolean, that can take as
                keyword-value pairs the same args, pars that are given to func, and also takes keyword-value "cost"
                and keyword_value "func_out". The value of cost will be the latest scalar value of func, and the
                value of "func_out" will be a tuple with any further returns from func. `stopping_function()` is run
                each iteration, and if it returns true, the minimization is stopped.

- verbose=false   If true, print out a report on each iteration of iteration number, radius size (eta),
                what type jump was proposed ("Newton" means going straight to global min, "constrained" means jump has
                norm eta, failed means that finding the minimum at a given radius somehow didn't work). Will also
                print out the cosine of the angle between the proposed jump and the gradient.

- verbose_timestamp   If true, adds a time and date stamp to the verbose report

- verbose_level   If less than 2, regular verbose output, if 2 or greater, very verbose, for debugging.

- verbose_file  If other than stdout, should be a string indicating a filename that verbose
                output will be written to.

- softbox       If true, then bbox must be a Dict() and we use the tanh() mechanism for putting a fixed limit
                on the parameters. NO LONGER SUPPORTING ANYTHING OTHER THAN softbox=true (which is the default)

- report_file   If non-empty, at each iteration timestep will write into this file outputs trajectory,
                (which contains eta, cost, and parameters), cpm_traj, and ftraj (which contains gradient, hessian,
                and further cost function outputs).  The file must be a JLD file, and so will end with a .jld extension.
                To load the saved dictionary, simply do D = load(filename)  (we have already called "using JLD" for you.)



# RETURNS:

- params       A vector the size of seed that has the last values of the minimizing parameters for func
- trajectory   A (2+length(params))-by-nsteps matrix. Each column corresponds to an iteration step, and contains
                 the value of eta used, the cost, and the value of the parameters at that iteration
- cost         Final value of objective function
- cpm_traj     A 4-by-nsteps matrix, containing reports from the contrained parabolic minimization at each timestep.
             The first row is niters (how many iterations cpm's 1-d minimization ran for) and the second row is
            Dlambda, the last change in the parameter being minimized in cpm's internal search,
            the third row is the squared difference between the returned and desired radius (should be very small),
            and the fourth row is cost change expected under the quadratic approximation

- ftraj     Further components for the trajectory, will be an Array{Any}(3, nsteps). First row is gradient,
            second row is Hessian, third row is second-and-further outputs of func, each one at each step of
            the minimization. **NOTE** that if these further outputs contain variables that are being minimized,
            they'll come out as ForwardDiff Duals, which you might not want!  So, for example, you might want to
            convert vectors and matrices into Float64s before returning them in those extra outputs. E.g.,
            if you want to return sum(err.*err) as the scalar to be minimized, and also return err, in your
            cost function you would write   " return sum(err.*err), Array{Float64}(err) ".   That way the first,
            scalar output can still be differentiated, for minimization, and the second one comes out in readable form.



# EXAMPLE:  (see also a more complete example in Optimization Utilities.ipynb)

```
function tester(;x=5, y=10, z=20, nderivs=0, difforder=0)
    return x^2*y + z/tanh(y)
end

params, trajectory = bbox_Hessian_keyword_minimization([0.5, 0.5], ["x", "y"], [1.1 2 ; 1.1 4], tester,
    verbose=true, tol=1e-12, start_eta=1);
```


"""
function bbox_Hessian_keyword_minimization(seed, args, bbox, func; start_eta=0.1, tol=1e-6, maxiter=400,
    frac_cost_threshold = 0.5, stopping_function = nothing,
    verbose=false, verbose_level=1, verbose_every=1, verbose_file=stdout,
    verbose_timestamp = false, start_iter_num = 1,
    softbox=true, hardbox=false, report_file="")

    # --- check that saving will be done to a .jld file ---
    if length(report_file)>0 && splitext(report_file)[2] != ".jld"
        if splitext(report_file)[2] == ""
            report_file = report_file * ".jld"
        else
            error("Sorry, report_file can only write to JLD files, the extension has to be .jld")
        end
    end


    # --------- Initializing the trajectory trace and function wrapper--------

    traj_increment = 100
    params = 0  # Make sure to have this here so that params stays defined beyond the try/catch
    if ( !(typeof(bbox)<:Dict) ); error("Currently only supporting softbox=true, bbox must be a Dict"); end;
    try
        params = copy(inverse_wall(bbox, args, seed))
    catch
        error("Were all initial param values within the indicated walls?")
    end
    eta = start_eta
    trajectory = zeros(2+length(params), traj_increment); cpm_traj = zeros(4, traj_increment)

    ftraj = Array{Any}(undef,3,0)  # will hold gradient, hessian, and further_out,  per iteration

    further_out =[];  # We define this variable here so it will be available for stashing further outputs from func
    stopping_func_out = false;   # Default value of stopping_func()

    # Now we define a wrapper around func() to do three things: (a) wallwrap parameters using the softbox method;
    # (b) return as the desired scalar the first output of func; (c) stash in further_out any further outputs of func
    internal_func = (;pars...) -> begin
        fresults = func(;wallwrap(bbox, pars)...)   # note use of bbox external to this begin...end
        if typeof(fresults)<:Tuple
            answer = fresults[1]
            further_out = fresults[2:end]
        else
            answer = fresults
        end
        return answer  # we assume that the first output of func() will always be a scalar, and that's what we return for ForwardDiff
    end

    # --------- END Initializing the trajectory trace --------

    if verbose
        if verbose_file==stdout; ostr=stdout else ostr=open(verbose_file, "a"); end
        @printf ostr "%d: eta=%g ps=" start_iter_num-1 eta
        println(ostr, vector_wrap(bbox, args, params))
        if verbose_file!=stdout; close(ostr); end
    end

    if softbox
        if !(typeof(bbox)<:Dict); error("bhm: If softbox=true, then bbox must eb a Dict"); end
        cost, grad, hess = keyword_vgh(internal_func, args, params)  # further_out will be mutated
    elseif hardbox
        error("Sorry, no longer supporting hardbox=true")
    else
        error("Sorry, no longer supporting softbox=false")
    end

    chessdelta = zeros(size(params))
    expected_cost_delta = 0   # this will be the quadratically expected cost change, either the Newton prediction or the constrained prediction
    hess_cost_delta     = 0   # this will hold the Newtorn prediction
    chess_cost_delta    = 0   # this will hold the constrained prediction

    my_iter=new_params=new_cost=new_grad=new_hess=0  # here so these variables are available outside the loop
    for i in [start_iter_num:(maxiter+start_iter_num);]
        my_iter = i
        while i > size(trajectory, 2)
            trajectory = [trajectory zeros(2+length(params), traj_increment)]
            cpm_traj   = [cpm_traj   zeros(4, traj_increment)]
        end
        trajectory[1:2, i]   = [eta;cost]
        trajectory[3:end, i] = vector_wrap(bbox, args, params)
        ftraj = [ftraj [grad, hess, further_out]]

        if length(report_file)>0
            save(report_file, Dict("traj"=>trajectory[:,1:i], "cpm_traj"=>cpm_traj[:,1:i], "ftraj"=>ftraj))
        end

        hessdelta  = - inv(hess)*grad
        hess_cost_delta = 0.5 * hessdelta' * hess * hessdelta + grad'*hessdelta   # Netwon prediction for how much cost should change
        try
            if verbose && verbose_level >= 2
                if verbose_file==stdout; ostr=stdout else ostr=open(verbose_file, "a"); end
                @printf(ostr, "bhm: about to try cpm with grad : "); print(ostr, grad); print(ostr, "\n")
                @printf(ostr, "bhm:   hess :"); print(ostr, hess[:]); print(ostr, "\n");
                if verbose_file!=stdout; close(ostr); end
            end
            if verbose && verbose_level >= 2
                cpm_out = constrainedParabolicMinimization(hess, grad'', eta,
                    maxiter=500, tol=1e-20, do_plot=true, verbose=true)
            else
                cpm_out = constrainedParabolicMinimization(hess, grad'', eta, maxiter=500, tol=1e-20)
            end
            chess_cost_delta = cpm_out[2]
            chessdelta = cpm_out[1];

            cpm_traj[1,i] = cpm_out[5]; cpm_traj[2,i] = cpm_out[6]; #  niters and Dlambda
            cpm_traj[3,i] = cpm_out[4]; cpm_traj[4,i] = cpm_out[2]; #  radius error and quadratically predicted change in J
            jumptype = "not failed"
        catch y
            if isa(y, InterruptException); throw(InterruptException()); end  # External interrupts should not be catchable
            jumptype = "failed"
            if verbose
                if verbose_file==stdout; ostr=stdout else ostr=open(verbose_file, "a"); end
                @printf ostr "Constrained parabolic minimization failed with error %s\n" y
                @printf ostr "\n"
                @printf ostr "eta was %g\n" eta
                @printf ostr "grad was\n"
                print(ostr, grad)
                @printf ostr "\n\nhess was\n"
                for k in [1:length(grad);]
                    print(ostr, hess[k,:])
                    @printf ostr "\n"
                end
                @printf ostr "\n"
                matwrite("error_report.mat", Dict("grad"=>grad, "hess"=>hess, "eta"=>eta))
                if verbose_file!=stdout; close(ostr); end
            end
            break
        end

        if norm(hessdelta) <= eta
            new_params = params + hessdelta
            jumptype = "Newton"
            expected_cost_delta = hess_cost_delta
        elseif jumptype != "failed"
            new_params = params + chessdelta
            jumptype  = "constrained"
            expected_cost_delta = chess_cost_delta
        end

        if jumptype != "failed"
            new_cost, new_grad, new_hess = keyword_vgh(internal_func, args, new_params)   # further_out may mutate
            if stopping_function != nothing
                stopping_func_out = stopping_function(; cost=new_cost, func_out=further_out,
                    make_dict(args, new_params)...)
            end
            if verbose && verbose_level >=2
                if verbose_file==stdout; ostr=stdout else ostr=open(verbose_file, "a"); end
                @printf(ostr, "bhm: had new_params = : "); print(ostr, vector_wrap(bbox, args, params)); print(ostr, "\n");
                @printf(ostr, "bhm: and my bbox was : "); print(ostr, bbox); print(ostr, "\n")
                @printf(ostr, "bhm: and my wallwrap output was : "); print(ostr, wallwrap(bbox, make_dict(args, new_params))); print(ostr, "\n")
                @printf(ostr, "bhm: and this produced new_grad : "); print(ostr, new_grad); print(ostr, "\n")
                @printf(ostr, "bhm:   new_hess :"); print(ostr, new_hess[:]); print(ostr, "\n");
                if verbose_file!=stdout; close(ostr); end
            end

            if abs(new_cost - cost) < tol || eta < tol || stopping_func_out
                if verbose
                    if verbose_file==stdout; ostr=stdout else ostr=open(verbose_file, "a"); end
                    @printf(ostr, "About to break -- stop_func_out = %s, tol=%g, new_cost-cost=%g, eta=%g\n",
                        stopping_func_out, tol, new_cost-cost, eta)
                if verbose_file!=stdout; close(ostr); end
                end
                break
            end
        end

        if jumptype == "failed" || cost <= new_cost || (new_cost - cost)/expected_cost_delta <= frac_cost_threshold
            if verbose
                if verbose_file==stdout; ostr=stdout else ostr=open(verbose_file, "a"); end
                @printf(ostr, "eta going down: ")
                if jumptype=="failed"; @printf(ostr, "jtype=failed");
                else                   @printf(ostr, "cost (new-old)/expect = %.3f", (new_cost - cost)/expected_cost_delta)
                end
                @printf(ostr, " new_cost-cost=%g and jumptype='%s'\n", new_cost-cost, jumptype)
                if verbose_level >= 2
                    nwp = vector_wrap(bbox, args, new_params); wp = vector_wrap(bbox, args, params)
                    @printf(ostr, "   vvv: proposed new params were : "); print(ostr, nwp); print(ostr, "\n")
                    @printf(ostr, "   vvv: proposed delta params was : "); print(ostr, nwp-wp); print(ostr, "\n")
                    @printf(ostr, "   vvv: grad was : "); print(ostr, grad); print(ostr, "\n")
                    costheta = dot(new_params-params, grad)/(norm(new_params-params)*norm(grad))
                    @printf(ostr, "   vvv: costheta of proposed jump was %g\n", costheta)
                end
                if verbose_file!=stdout; close(ostr); end
            end
            eta = eta/2
            costheta = NaN
            if eta < tol || stopping_func_out
                if verbose
                    if verbose_file==stdout; ostr=stdout else ostr=open(verbose_file, "a"); end
                    @printf(ostr, "About to break -- stop_func_out = %s, tol=%g, new_cost-cost=%g, eta=%g\n",
                        stopping_func_out, tol, new_cost-cost, eta)
                    if verbose_file!=stdout; close(ostr); end
                end
                break
            end
        else
            eta = eta*1.1
            costheta = dot(new_params-params, grad)/(norm(new_params-params)*norm(grad))

            params = new_params
            cost = new_cost
            grad = new_grad
            hess = new_hess
        end

        if verbose
            if rem(i, verbose_every)==0
                if verbose_file==stdout; ostr=stdout else ostr=open(verbose_file, "a"); end
                if verbose_timestamp; @printf ostr "[%s] " Dates.format(now(), "e, dd u yyyy HH:MM:SS"); end
                @printf ostr "%d: eta=%g cost=%g jtype=%s costheta=%.3f ps=" i eta cost jumptype costheta
                print(ostr, vector_wrap(bbox, args, params))
                @printf ostr "\n"
                if verbose_level >= 3
                    @printf ostr "    At this point, grad is ="
                    print(ostr, grad)
                    @printf ostr "\n"
                end
                if verbose_file!=stdout; close(ostr); end
            end
        end
    end
    i = my_iter

    if stopping_func_out
        # The last params were the good ones!
        params = new_params
        cost = new_cost
        grad = new_grad
        hess = new_hess

        i = i+1
        trajectory[1:2, i]   = [eta;cost]
        trajectory[3:end, i] = vector_wrap(bbox, args, params)
        ftraj = [ftraj [grad, hess, further_out]]
    end


    trajectory = trajectory[:,1:i]; cpm_traj = cpm_traj[:,1:i]
    if length(report_file)>0
        JLD.save(report_file, Dict("traj"=>trajectory, "cpm_traj"=>cpm_traj, "ftraj"=>ftraj))
    end

    return vector_wrap(bbox, args, params), trajectory, cost, cpm_traj, ftraj
end


using Random

@doc """
    testMe()

    Run this function (no parameters, no returns) to test OptimizationUtils
""" function testMe()
    sr = Int64(round(time()*1000))
    # sr = 1510162002784 # For these values of sr
    # sr = 1509561656447 # when start_eta=1, the threshold quickly goes very positive and the minimization gets stuck
    #
    # sr = 1510164239381   # For this value, it gets stuck at a very small inverse slope param... don't know why
    Random.seed!(sr)


    npoints = 1000; # srand(400)
    args = ["baseline", "amplitude", "threshold", "slope"]

    # Generating values for our four params:
    params = [1 5 0.5 0.8]

    # Make some points and plot them
    x = rand(npoints, 1)*6 .- 3
    y = params[1] .+ params[2]*0.5*(tanh.((x.-params[3])./params[4]).+1) .+ randn(npoints,1)*2
    figure(1); clf();
    plot(x, y, ".")

    # Starting values for the four params. Plot the corresponding curve they generate
    seed = [8, 3.1, 0, 0.02]
    xx = -3:0.01:3
    plot(xx, seed[1] .+ seed[2]*0.5*(tanh.((xx.-seed[3])/seed[4]).+1), "g-")

    # Cost function.  First output is the scalar
    # that will be minimized, and we also returns a second output whose trajectory will be stashed
    # by bbox in ftraj as a diagnostic during the minimization.
    function myJJ(x, y; baseline=0, amplitude=1, threshold=0, slope=1, do_plot=false, fignum=1, clearfig=true)

        if do_plot
            figure(fignum);
            if clearfig; clf(); end;
            xx = -3:0.01:3; x2=zeros(get_eltype((baseline,amplitude,threshold,slope)), size(xx,1), size(xx,2))
            for i=1:length(xx); x2[i]=xx[i]; end; xx= x2

            plot(x, y, ".")
            plot(xx, baseline .+ amplitude*0.5*(tanh.((xx.-threshold)/slope).+1), "r-")
        end

        yhat =  baseline .+ amplitude*0.5*(tanh.((x.-threshold)./slope).+1)
        err = yhat .- y
        return sum(err.*err), get_value(err)    # Note first output, the scalar to be minimized,
        # may be ForwardDiff Duals during the minimization, which is fine, so it can be differentiated.
        # The second one we use get_value to turn into regular Float64 array so it comes out readable.
    end



    if ~isdir("Trash"); mkdir("Trash"); end;  # we're going to put the iteration-step by iteration-step report file there

    bbox = Dict(:baseline=>[-2, 10], :slope=>[0.001 5])
    func = (;pars...) -> myJJ(x, y; do_plot=false, pars...)

    stopping_func = (;cost=0, func_out=[], pars...) -> return cost<1500;   # Make that a high number and it'll stop early

    opars, traj, cost, cpm_traj, ftraj = bbox_Hessian_keyword_minimization(seed, args, bbox, func,
        frac_cost_threshold = 0.5, stopping_function = stopping_func,
        verbose=true, verbose_level=1, verbose_file="tester.txt",
        softbox=true, start_eta=0.1, report_file="Trash/example_report.jld")

    # Note that the gradient at step i of the minimization will be available as ftraj[1,i], the hessian will be
    # in ftraj[2,i], and the error vector, which is the first of the extra outputs of myJJ(), will be in ftraj[3,i][1].
    # In our example myJJ() produced only one extra output; a second extra output would be in ftraj[3,i][2], and so on.

    # Plot the resulting curve, and report both final and generating params
    figure(1);
    plot(xx, opars[1] .+ opars[2]*0.5*(tanh.((xx.-opars[3])/opars[4]).+1), "r-")
    [opars' ; params]
    xlabel("x"); ylabel("y"); title("green is sigmoid with starting params, red is end")


    figure(2); clf();
    ax1 = subplot(2,1,1)
    plot(cpm_traj[4,:], ".-")
    plot(traj[2,2:end] - traj[2,1:end-1], ".-")
    grid("on")
    gca().set_xticklabels("")
    legend(["expected cost change", "actual cost change"])

    subplot(2,1,2)
    plot((traj[2,2:end] - traj[2,1:end-1])./cpm_traj[4,1:end-1], ".-")
    plot(traj[1,2:end]./traj[1,1:end-1], ".-")
    grid("on")

    legend(["actual/expected cost change", "fractional change in eta"])

    println("Final costs were: ", traj[2,end-3:end]);
end

end # ===== END MODULE =========
