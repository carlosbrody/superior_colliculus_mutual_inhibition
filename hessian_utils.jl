# DON'T MODIFY THIS FILE -- the source is in file Cost Function Minimization and Hessian Utilities.ipynb


"""
function value, gradient, hessian = vgh(func, x0)

Wrapper for ForwardDiff.hessian!() that computes and returns all three of a function's value, gradient, and hessian.

EXAMPLE:
========

function tester(x::Vector)

    return sum(x.*x)
end

value, grad, hess = vgh(tester, [10, 3.1])
"""
function vgh(func, x0)
    out = DiffBase.HessianResult(x0)             
    ForwardDiff.hessian!(out, func, x0)
    value    = DiffBase.value(out)
    gradient = DiffBase.gradient(out)
    hessian  = DiffBase.hessian(out)
    
    return value, gradient, hessian    
end




"""
function x, cost, iters_used, last_Delta_x = one_d_minimizer(seed, func; tol=1e-5, maxiter=100, start_eta=0.1)

Minimizes a 1-d function using constrained Hessian minimization. 
We don't trust the long-range info from the Hessian too much, meaning that there's
a given (adaptive) step size. If Newton's method suggests a step smaller than that step
size, we take it. Otherwise, we only move in the direction of the gradient by the stepsize.

Adaptive step size: Every step that the cost function gets smaller, the step grows by a factor 
of 1.1. If a step would have led to a larger cost function, the step is not taken, and
the size falls by a factor of 2.

PARAMETERS:
===========

seed       A float, the starting point for the minimization

func       A one dimensional function, takes a float and returns a float that represents the current cost.


OPTIONAL PARAMETERS:
====================

tol=1e-5   If a step would lead to a change in cost that is smaller in magnitude than tol, stop the minimization.

maxiter=100   Maximum number of iterations for the minimization

start_eta=0.1  The starting value of the step size


RETURNS:
========

x0     the value of x that minimizes f(x)

cost   The value of f(x0) 

niters  The number of iterations done

Dparam  The last change to the parameter value

"""
function one_d_minimizer(seed, func; tol=1e-5, maxiter=100, start_eta=0.1)
    eta = start_eta;
    lambdavec = [seed]


    out = DiffBase.HessianResult(lambdavec)
    ForwardDiff.hessian!(out, func, lambdavec)
    cost = DiffBase.value(out)
    grad = DiffBase.gradient(out)
    hess = DiffBase.hessian(out)

    i = delta_lambda = 0;  # declare it so we can use it after the if
    for i in [1:maxiter;]
        h_delta = - grad./hess;
        if abs(h_delta[1]) < eta
            new_lambdavec = lambdavec + h_delta
        else
            new_lambdavec = lambdavec - eta*sign(grad)
        end

        delta_lambda = new_lambdavec[1] - lambdavec[1]
        if abs(delta_lambda) < tol
            break
        end
        
        ForwardDiff.hessian!(out, func, new_lambdavec)
        new_cost = DiffBase.value(out)
        new_grad = DiffBase.gradient(out)
        new_hess = DiffBase.hessian(out)

        if new_cost .< cost
            lambdavec[1] = new_lambdavec[1]
            cost = new_cost;
            grad = new_grad;
            hess = new_hess;

            eta = eta*1.1
        else
            eta = eta/2
        end

        # @printf "%d: cost=%.3f lambda=%.3f\n" i cost lambdavec[1]
    end
    
    return lambdavec[1], cost, i, delta_lambda
end
    


"""
function constrained_parabolic_minimization(H, G, r; tol=1e-6, min_only=true, do_plot=false, 
    verbose=false, efactor=3.0, max_efactor_tries=10, 
    lambdastepsize=0.003, minimum_tol=1e-24, tol_delta=1e-3, maxiter=200)

Given a Hessian matrix, a gradient vector, and a desired radius from the origin, finds the vector 
that minimizes the parabola defined by the Hessian and the gradient, subject to the constraint that the
vector's length equals the desired radius.

PARAMETERS:
===========

H      A square symmetric matrix. It should have all positive eigenvalues.

G      A vector, length equal to the size(H,2)

r      desired radius

OPTIONAL PARAMETERS:
====================

tol=1e-6        Numerical tolerance on the computations

min_only=true   Return only the minimum, or, if false, all xs for all lambdas that match x'*x = r^2

efactor=3       The initial exploration of lambdas will go from -efactor(max(absolute(eig(H)))) to +efactor(max(absolute(eig(H))))
                If that does not produce a solution within the indicates tolerance, then efactor is multiplied by efactor_growth
                and we try again, up to max_efactor_tries.

efactor_growth  see efactor

max_efactor_tries   see efactor

minimum_tol     Each candidate value of lambda, identified as a minimum in the grid scan and 
                proposed as satisfying |x|^2=r^2, is then used to seed a Newton's method one-d minimization with a certain tolerance.
                If the results of that search do not satisfy |x|^2=r^2 within the desired tolerance, then the tolerance
                for the 1-d search is reduced (to make the one-d search more exact), by a factor of tol_delta, 
                up to minimu_tol

tol_delta       see minimum_tol

max_iter        Maximum number of iterations in the 1-d search (see minimum_tol).


lambdastepsize=0.003    The step size for initial exploration of lambdas, un units of efactor. It sshould
                probably scale with the smallest difference in the eigenvalues of H; that has not been implemented yet.

do_plot         If true, produces plot of the initial gridscans of lambda versus (x'*x - r^2)^2

verbose         If true, produces diagnostic info as it goes.


RETURNS:
========

x        The vector that minimizes 0.5*x'*H*x + x'*G subject to x'*x = r
 
J        0.5*x'*H*x + x'*G at the returned x

lambda   value of the Lagrange multiplier at which the radius constraint is satisfied

c        The squared difference between the length of x and r. Should be small, otherwise somthing went wrong!

niters   The number of iterations used in the one-d minimiation that identified the output lambdas.

Dlambda  The last change in lambda during the one-d minimization. If the one-d minimization did not reach
         its maxiters, then this will be smaller than the one-d minimization's tolerance.

"""
function constrained_parabolic_minimization(H, G, r; tol=1e-6, min_only=true, 
    do_plot=false, verbose=false, efactor=3.0, efactor_growth=4, max_efactor_tries=10, 
    lambdastepsize=0.003, minimum_tol=1e-24, tol_delta=1e-3, maxiter=200)

    #  --- First a couple of helper functions ----
    
    """
    function x_of_lambda(lambda)

    Given square matrix H, vector G, and passed scalar lambda, returns the vector x that minimizes
    
    0.5 x'*H*x + x'*G - lambda *x'*x

    """
    function x_of_lambda(lambda)
        return inv(H - lambda*eye(size(H,1)))*(-G)
    end
    
    
    """
    function q(lambda, r)

    Returns the squared difference between r and the norm of x_of_lambda(lambda).
    """
    function q(lambda, r)
        return (r - norm(x_of_lambda(lambda)))^2
    end


    # efactor is the factor that multiplies the biggest eigenvalue of H, to determine the range over which we'll
    # look for a lambda that satisfies the norm(x)==r requirement. If we don't find a solution, we iteratively 
    # increase efactor to try to get there, for a maximum of max_efactor_tries
    for m=1:max_efactor_tries
        # First scan lambda to find good candidates for minimizing the parabolic 
        # surface under the x'*x = r^2 constraint
        L = eig(H)[1]
        L0 = maximum(abs(L))
        lambdas = L0*efactor*[-1.0:lambdastepsize:1.0;]
        costs = zeros(size(lambdas))
        for i in [1:length(lambdas);]
            try 
                costs[i] = q(lambdas[i], r)
            catch
                costs[i] = Inf
            end
        end

        if do_plot
            figure(2); clf();
            plot(lambdas, costs, "b.-")
            xlabel("lambda")
            ylabel("cost")
        end

        # Take all candidates where the derivative of costs changes sign 
        # from negative to positive (those would be minima),
        # plus the smallest and the largest lambdas tested, as candidates
        g = append!(prepend!(find(diff(sign(diff(costs))) .> 0.99), [1]), [length(lambdas)])
        if verbose
            @printf("cpm: g (candidate indices) are : ");           print_vector_g(g);        print("\n")
            @printf("cpm: and their corresponding costs are : ");   print_vector(costs[g]);   print("\n");
            @printf("cpm: and their corresponding lambdas are : "); print_vector(lambdas[g]); print("\n");
            @printf("cpm: the minimum cost was : %g\n", minimum(costs[g]))
        end
        # found_it_flag = 0  # A flag for when we've found at least one lambda that satisfies the r constraint
        mytol = tol

        while mytol > minimum_tol
            lambdas_out = zeros(size(g))
            costs_out   = zeros(size(g))
            niters_out  = zeros(size(g))
            Dlambda_out = zeros(size(g))
            for i in [1:length(g);]
                lambdas_out[i], costs_out[i], niters_out[i], Dlambda_out[i] = 
                one_d_minimizer(lambdas[g[i]], x -> q(x[1], r), start_eta=1, tol=mytol, maxiter=maxiter)
            end

            # Eliminate any lambdas where x'*x doesn't match our desired value r
            I = find(costs_out .< tol)
            lambdas_out = lambdas_out[I]; costs_out = costs_out[I];
            niters_out  = niters_out[I];  Dlambda_out = Dlambda_out[I]

            if length(I) > 0; break; end

            mytol *= tol_delta
        end
        if verbose
            @printf("%d : After searching for lambdas with efactor=%g, we found these : ", m, efactor)
            print_vector_g(lambdas_out); print("\n")
        end
        if length(lambdas_out) > 0; break; end;
        efactor = efactor*efactor_growth
    end
    
    # Eliminate any repeated lambdas, to within the specified numerical tolerance.
    I = setdiff(1:length(lambdas_out), find(diff(lambdas_out) .< tol))
    lambdas_out = lambdas_out[I]; costs_out = costs_out[I];
    niters_out  = niters_out[I];  Dlambda_out = Dlambda_out[I]
    
    # Find the parabolic estimate of the cost function at these points
    J  = zeros(size(lambdas_out))
    xs = zeros(length(G), length(lambdas_out))
    for i in [1:length(J);]
        xs[:,i] = x_of_lambda(lambdas_out[i])
        J[i] = (0.5*xs[:,i]'*H*xs[:,i] + xs[:,i]'*G)[1]
    end

    # Find and return only the x that has the smallest J
    if min_only
        I = indmin(J)    
    else
        I = 1:length(J)
    end
    return xs[:,I], J[I], lambdas_out[I], costs_out[I], niters_out[I], Dlambda_out[I]
end


    

    


"""
function constrained_Hessian_minimization(seed, func; start_eta=10, tol=1e-6, maxiter=400,
    verbose=false)

Minimizes a function using. At each step, it computes the Hessian approximation to the function, 
and then asks, according to the corresponding parabolic approximation: is the global minimum within a
radius eta? If so, try to jump to it. If not, try to jump to the point at a distance eta from the 
current point that would minimize the function.  If the attempted jump leads to an increase in the function,
then the jump is rejected and eta is reduced by a factor of 2. If the attempted jump reduces the function,
then it is accepted and eta is increased by a factor of 1.1.  This proceeds until the change in the function
after a proposed jump would be less than tol, or the iteration number has reached maxiter, whichever happens
first.  Returns the minimizing parameters.

PARAMETERS:
===========

seed        column vector, representing the starting value of the parameters.

func        Function that takes a vector and returns a scalar.  If you want to
            work with a function that tales mpre parameterrs and returns more than one 
            output, you can use something like

                    x -> orig_func(x, other_params)[1]

            You only need the "[1]" part of the orig_func returns more outputs than a scalar. 

OPTIONAL PARAMETERS:
====================

start_eta=10 Starting value of the radius.  It's good to start with somethibg biggish, if it is
             too much, it'll quickly get cut down.

tol=1e-6     Numerical tolerance. If a proposed jump produces a change in func that is less than
             this, the minimization stops.

maxiter=400  Maximum number of iterations to do before stopping

verbose=false   If true, print out a report on each iteration of iteration number, radius size (eta),
                what type jump was proposed ("Newton" means going straight to global min, "constrained" means jump has 
                norm eta, failed means that finding the minimum at a given radius somehow didn't work). Will also
                print out the cosine of the angle between the proposed jump and the gradient.

RETURNS:
========

params       A vector the size of seed that has the last values of the minimizing parameters for func

"""
function constrained_Hessian_minimization(seed, func; start_eta=10, tol=1e-6, maxiter=400,
    verbose=false)

    params = seed
    eta = start_eta

    out = DiffBase.HessianResult(params)
    ForwardDiff.hessian!(out, func, params)
    cost = DiffBase.value(out)
    grad = DiffBase.gradient(out)
    hess = DiffBase.hessian(out)

    chessdelta = zeros(size(params))

    for i in [1:maxiter;]
        hessdelta  = - inv(hess)*grad
        try
            chessdelta = constrained_parabolic_minimization(hess, grad'', eta)[1]
            jumptype = "not failed"
        catch
            jumptype = "failed"
        end

        if norm(hessdelta) <= eta
            new_params = params + hessdelta
            jumptype = "Newton"
        elseif jumptype != "failed" 
            new_params = params + chessdelta
            jumptype  = "constrained"
        end

        if jumptype != "failed"
            ForwardDiff.hessian!(out, func, new_params)
            new_cost = DiffBase.value(out)
            new_grad = DiffBase.gradient(out)
            new_hess = DiffBase.hessian(out)
            
            if abs(new_cost - cost) < tol
                break
            end
        end

        if jumptype == "failed" || new_cost >= cost
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
            @printf "%d: eta=%.3f cost=%.4f jtype=%s costheta=%.3f ps=" i eta cost jumptype costheta
            print_vector(params)
            @printf "\n"
        end
    end
    
    return params
end



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
        print_vector_g(:cost)
        print_vector_g(:grad)
        print_vector_g(:hess)
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
            print_vector(params)
            @printf "\n"
        end
    end
    
    return params, cost
end


#############################################################################
#                                                                           #
#                   KEYWORD GRADIENTS AND HESSIANS                          #
#                                                                           #
#############################################################################
        
"""
dict = make_dict(argstrings, x, [starting_dict=Dict()] )

Given a list of strings, and a list of values, makes a dictionary of Symbols to values, with the Symbols 
corresponding to each of the strings.  Mostly used to pass arguments as a keyword-value set into a function.
If one of the elements of argstrings is *not* a string, but is instead a 2-long list, the first element of that 
list should be a string, and the second element of that list should be a positive integer. This will be 
interpreted as "don't take only one value, take this number of values and this parameter will be a vector"

PARAMS:
=======

argstrings       A list of strings. Each element may also be a two-long list of a string, positive integer,
                 e.g., ["this" 3]

x                A vector of numeric values. Its length must be such that all the strings in argstrings
                 can take their corresponding element(s), sequentially, from x

starting_dict    An initial dictionary to work with.  Any key in this starting dictionary matching an argstring
                 will be replaced by the new value. Keys not matched will remain.

RETURNS:
========

dict             The symbol dictionary.


EXAMPLES:
=========

>> make_dict(["this", "that", ["there", 2]], [10, 20, 3, 4])

Dict{Any,Any} with 3 entries:
  :this  => 10
  :that  => 20
  :there => [3,4]

>> make_dict(["doo", "gaa"], [10, 20], Dict(:blob=>100, :gaa=>-44))

Dict{Symbol,Int64} with 3 entries:
  :gaa  => 20
  :blob => 100
  :doo  => 10

"""
function make_dict(args, x, starting_dict=Dict())
    kwargs = starting_dict;
    i = 1; j=1
    while i<=length(args)
        if typeof(args[i])==String
            kwargs = merge(kwargs, Dict(Symbol(args[i]) => x[j]))
        else
            if length(args[i]) == 2
                extra = args[i][2]-1
                kwargs = merge(kwargs, Dict(Symbol(args[i][1]) => x[j:(j+extra)]))
                j = j+extra
            else
                error("Each element of the args vector must be either a string, or a 2-long vector, first element a string, second integer")
            end            
        end
        i = i+1; j=j+1
    end
    return kwargs
end 



"""
function M = ForwardDiffZeros(m, n; nderivs=0, difforder=0)

Use instead of zeros(). Creates a matrix of zeros, of size m rows by n columns, with elements appropriate for 
differentiation by ForwardDiff. If nderivs==0 or difforder==0 then the elements will be regular
Float64, not ForwardDiff types.

PARAMETERS:
===========

m        Integer, number of rows

n        Integer, number of columns


OPTIONAL PARAMETERS:
====================

nderivs=0       The number of variables that we'll be differentiating with respect to. In other
                words, this number is equal to the length of the gradient. If this is left as zero (the default) then 
                the data type will be regular Float64

difforder=0     The order of the derivative we will want to take.  Zero means nothing, stick with
                regular Float64, 1 means gradient, 2 means hessian

RETURNS:
========

An m-by-n matrix of zeros that can be used with Forward Diff.

"""
function ForwardDiffZeros(m, n; nderivs=0, difforder=0)
    if nderivs == 0 || difforder == 0
        return zeros(m, n)
    elseif difforder == 1
        return zeros(ForwardDiff.Dual{nderivs, Float64}, m , n)
    elseif difforder == 2
        return zeros(ForwardDiff.Dual{nderivs, ForwardDiff.Dual{nderivs, Float64}}, m, n)
    else
        error("Don't know how to do that order of derivatives!", nderivs)
    end
end
      

"""
function keyword_gradient(func, args, x0)

Same as ForwardDiff.gradient except that func() must be a function taking only optional 
keyword arguments, and the derivative is taken with respect to an arbitrarily chosen set of 
those, indicated by a list of strings.

In addition, func *MUST* take optional keyword args nderivs=0 and difforder=0, and within it,
if matrices or vectors of zeros are declared, use ForwardDiffZeros() instead of zeros().

PARAMETERS:
===========

func        A scalar function taking only optional keyword arguments, including nderivs=0 and difforder=0

args        A list of strings indicating which keyword arguments to differentiate. These strings must
            match the keyword names in func()   For example, func(;this=10, that=20) would mean that 
            "this" and "that" are allowable elements in args.

x0          A vector of floats, same length as args, representing the values of these args at which the
            derivatives will be taken.

RETURNS:
========

grad        The gradient of func w.r.t. args


EXAMPLE:
========

function tester(;a=10, b=20, c=30, nderivs=0, difforder=0)
    M = ForwardDiffZeros(3, 3; nderivs=nderivs, difforder=difforder)
    M[1,1] = a^2*10
    M[2,2] = b*20
    M[3,3] = a*sqrt(c)*30.1
    return trace(M)
end

grad_a_c = keyword_gradient((;pars...) -> tester(;pars...), ["a", "c"], [10, 3.1])  # note initial values must be floats

grad_b_c = keyword_gradient((;pars...) -> tester(;pars...), ["b", "c"], [10, 3.1]) 

"""
function keyword_gradient(func, args, x0)
    
    ans = ForwardDiff.gradient(x -> func(;nderivs=length(x0), difforder=1, make_dict(args, x)...), x0)
    
    return ans
end


"""
function keyword_gradient!(out, func, args, x0)

Same as keyword_gradient, but puts the result in mutable out. See keyword_gradient() for documentation.

EXAMPLE:
========

function tester(;a=10, b=20, c=30, nderivs=0, difforder=0)
    M = ForwardDiffZeros(3, 3; nderivs=nderivs, difforder=difforder)
    M[1,1] = a^2*10
    M[2,2] = b*20
    M[3,3] = a*sqrt(c)*30.1
    return trace(M)
end

out = DiffBase.GradientResult([10, 3.1])  # out must be same length as whatever we will differentiate w.r.t.
keyword_gradient!(out, (;pars...) -> tester(;pars...), ["a", "c"], [10, 3.1])  # note initial values must be floats
grad_a_c = DiffBase.gradient(out)
value    = DiffBase.value(out)

out = DiffBase.GradientResult([10, 3.1, 20])  # out must be same length as whatever we will differentiate w.r.t.
keyword_gradient!(out, (;pars...) -> tester(;pars...), ["a", "b", "c"], [10, 20, 3.1])  # note initial values must be floats
grad_a_b_c = DiffBase.gradient(out)

"""
function keyword_gradient!(out, func, args, x0)

    if length(args) != length(x0)
        error("Oy! args and x0 must be the same length!")
    end

    ForwardDiff.gradient!(out, x -> func(;nderivs=length(x0), difforder=1, make_dict(args, x)...), x0)
    
    return 
end


"""
function keyword_hessian(func, args, x0)

Same as ForwardDiff.hessian except that func() must be a function taking only optional 
keyword arguments, and the derivative is taken with respect to an arbitrarily chosen set of 
those, indicated by a list of strings.

In addition, func *MUST* take optional keyword args nderivs=0 and difforder=0, and within it,
if matrices or vectors of zeros are declared, use ForwardDiffZeros() instead of zeros().

PARAMETERS:
===========

func        A scalar function taking only optional keyword arguments, including nderivs=0 and difforder=0

args        A list of strings indicating which keyword arguments to differentiate. These strings must
            match the keyword names in func()   For example, func(;this=10, that=20) would mean that 
            "this" and "that" are allowable elements in args.

x0          A vector of floats, same length as args, representing the values of these args at which the
            derivatives will be taken.

RETURNS:
========

grad        The gradient of func w.r.t. args


EXAMPLE:
========

function tester(;a=10, b=20, c=30, nderivs=0, difforder=0)
    M = ForwardDiffZeros(3, 3; nderivs=nderivs, difforder=difforder)
    M[1,1] = a^2*10
    M[2,2] = b*20
    M[3,3] = a*sqrt(c)*30.1
    return trace(M)
end

hess_b_c = keyword_hessian((;pars...) -> tester(;pars...), ["b", "c"], [10, 3.1])  # note initial values must be floats

hess_a_b_c = keyword_hessian((;pars...) -> tester(;pars...), ["a", "b", c"], [10, 2, 3.1]) 

"""
function keyword_hessian(func, args, x0)

    if length(args) != length(x0)
        error("Oy! args and x0 must be the same length!")
    end

    ans = ForwardDiff.hessian(x -> func(;nderivs=length(x0), difforder=2, make_dict(args, x)...), x0)
    
    return ans
end


"""
function keyword_hessian!(out, func, args, x0)

Same as keyword_hessian, but puts the result in mutable out. See keyword_hessian() for documentation.

EXAMPLE:
========

function tester(;a=10, b=20, c=30, nderivs=0, difforder=0)
    M = ForwardDiffZeros(3, 3; nderivs=nderivs, difforder=difforder)
    M[1,1] = a^2*10
    M[2,2] = b*20
    M[3,3] = a*sqrt(c)*30.1
    return trace(M)
end

out = DiffBase.HessianResult([10, 3.1])  # out must be same length as whatever we will differentiate w.r.t.
keyword_hessian!(out, (;pars...) -> tester(;pars...), ["a", "c"], [10, 3.1])  # note initial values must be floats
hess_a_c = DiffBase.hessian(out)
grad_a_c = DiffBase.gradient(out)
value    = DiffBase.value(out)

out = DiffBase.HessianResult([10, 3.1, 20])  # out must be same length as whatever we will differentiate w.r.t.
keyword_hessian!(out, (;pars...) -> tester(;pars...), ["a", "b", "c"], [10, 20, 3.1])  # note initial values must be floats
hess_a_b_c = DiffBase.hessian(out)

"""
function keyword_hessian!(out, func, args, x0)
    nargs = 0
    for i in [1:length(args);]
        if typeof(args[i])==String
            nargs += 1
        else
            nargs += args[i][2]
        end
    end
    if nargs != length(x0)
        error("Oy! args and x0 must be the same length!")
    end

    ForwardDiff.hessian!(out, x -> func(;nderivs=length(x0), difforder=2, make_dict(args, x)...), x0)
    
    return 
end




"""
function value, gradient, hessian = keyword_vgh(func, args, x0)

Wrapper for keyword_hessian!() that computes and returns all three of a function's value, gradient, and hessian.

EXAMPLE:
========

function tester(;a=10, b=20, c=30, nderivs=0, difforder=0)
    M = ForwardDiffZeros(3, 3; nderivs=nderivs, difforder=difforder)
    M[1,1] = a^2*10
    M[2,2] = b*20
    M[3,3] = a*sqrt(c)*30.1
    return trace(M)
end

value, grad, hess = keyword_vgh(tester, ["a", "c"], [10, 3.1])
"""
function keyword_vgh(func, args, x0)
    out = DiffBase.HessianResult(x0)
    keyword_hessian!(out, func, args, x0)
    value    = DiffBase.value(out)
    gradient = DiffBase.gradient(out)
    hessian  = DiffBase.hessian(out)
    
    return value, gradient, hessian    
end




# DON'T MODIFY THIS FILE -- the source is in file Cost Function Minimization and Hessian Utilities.ipynb




"""
pdict = wallwrap(bdict, pdict)
Given bdict, a dictionary of symbols to [minval, maxval] vectors, and pdict, a dictionary of symbols
to values (or, alternatively, an Array of (Symbol, value) tuples], goes through each of the symbols in 
bdict and modifies the corresponding value in pdict putting it through a tanh so the final output lies 
within the limits in bdict.  Returns the new pdict.  Makes a copy of pdict so as not to modify the original.
"""
function wallwrap(bdict, epdict)
    local pdict = two_level_copy(epdict)  # Must be very careful here! I got bit by the bug of forgetting that without
    # an explicit copy() call, Julia does not make copies of the contents of arrays or dictionaries, making it
    # easy to inadvertently modify something one did not intend to perturb.  Note the two_level_copy() call, 
    # necessary to make sure we don't mess up the content of the caller's dictionary.
    
    if typeof(pdict)<:Array
        pdict = Dict(pdict)
    end

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
    local params = two_level_copy(eparams)
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
    local params = two_level_copy(wparams)
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
    local pdict = two_level_copy(wdict)

    allkeys = keys(bdict)
    for k in allkeys
        local bbox = bdict[k]
        d = 0.5*(bbox[2] - bbox[1])
        m = 0.5*(bbox[2] + bbox[1])

        pdict[k] = m + d*0.5*log((pdict[k]-bbox[1])./(2*d - pdict[k] + bbox[1]))
    end
    return(pdict)
end
  


# DON'T MODIFY THIS FILE -- the source is in file Cost Function Minimization and Hessian Utilities.ipynb



######################################################
#                                                    #
#         BBOX_HESSIAN_KEYWORD_MINIMIZATION          #
#                                                    #
######################################################




"""
function bbox_Hessian_keyword_minimization(seed, args, bbox, func; wallwidth=NaN, start_eta=10, tol=1e-6, 
    maxiter=400, verbose=false)

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

- func        func must take only optional keyword args, and must 
            take nderivs=0, difforder=0  and declare any new matrices using ForwardDiffZeros() instead of zeros()


# OPTIONAL PARAMETERS:

- start_eta    Starting value of the radius.  It's good to start with somethibg biggish, if it is
             too much, it'll quickly get cut down.

- tol=1e-6     Numerical tolerance. If a proposed jump produces a change in func that is less than
             this, the minimization stops.

- maxiter=400  Maximum number of iterations to do before stopping

- verbose=false   If true, print out a report on each iteration of iteration number, radius size (eta),
                what type jump was proposed ("Newton" means going straight to global min, "constrained" means jump has 
                norm eta, failed means that finding the minimum at a given radius somehow didn't work). Will also
                print out the cosine of the angle between the proposed jump and the gradient.

- verbose_level   If less than 2, regular verbose output, if 2 or greater, very verbose, for debugging.

- softbox         If true, then bbox must be a Dict() and we use the tanh() mechanism for putting a fixed limit
                on the parameters.

- hardbox=false   If true, ignores wallwidth, and just rests parameter values to the bounding box if they go outside it.
                If false, adds cost function "walls" to implement the bounding box.

- walldith=NaN     Used for putting up cost function "walls" that implement the bounding box limits. Can be NaN.
                If it is NaN, then the wallwidth is a constant factor of the range width for each argument. If not NaN, must
                be an nargs-long vector that indicates the actual wall widths.

- wallwidth_factor=0.18   Only relevant if wallwidth is NaN, otherwise ignored. For each arg, the wall width
                is going to be wall_width_factor*(bbox[i,2] - bbox[i,1])


# RETURNS:

- params       A vector the size of seed that has the last values of the minimizing parameters for func
- trajectory   A (2+length(params))-by-nsteps matrix. Each column corresponds to an iteration step, and contains
                 the value of eta used, the cost, and the value of the parameters at that iteration
- cost         Final value of objective function
- cpm_traj     A 2-by-nsteps matrix, containing reports from the contrained parabolic minimization at each timestep.
             The first row is niters (how many iterations cpm's 1-d minimization ran for) and the second row is
             Dlambda, the last change in the parameter being minimized in cpm's internal search


# EXAMPLE:

```
function tester(;x=5, y=10, z=20, nderivs=0, difforder=0)
    return x^2*y + z/tanh(y)
end

params, trajectory = bbox_Hessian_keyword_minimization([0.5, 0.5], ["x", "y"], [1.1 2 ; 1.1 4], tester, 
    verbose=true, tol=1e-12, start_eta=1);
```


"""
function bbox_Hessian_keyword_minimization(seed, args, bbox, func; start_eta=0.1, tol=1e-6, maxiter=400,
    verbose=false, verbose_level=1, verbose_every=1, 
    softbox=true, hardbox=false, wallwidth=NaN, wallwidth_factor=0.18)

      
    """
    Given args, a list of string representing the arguments of interest, a bounding box for each,
    and a Symbol=>value dictionary with the corresponding parameters, computes and returns a high cost for 
    being outside the bounding box
    """
    function wall_cost(args, bbox; wallwidth=NaN, nderivs=0, difforder=0, pars...) 
        myparams = ForwardDiffZeros(length(pars), 1, nderivs=nderivs, difforder=difforder)
        pars2 = Dict()
        for i in [1:length(pars);]
            pars2[string(pars[i][1])] = pars[i][2]
        end
        for i in [1:length(args);]
            myparams[i] = pars2[args[i]]
        end
        
        if isnan(wallwidth)
            # We know that we're going to be taking hessian for params, so declare zeros accordingly:
            wallwidth = ForwardDiffZeros(length(myparams), 1, nderivs=nderivs, difforder=difforder)

            for i in [1:length(myparams);]
                wallwidth[i] = wallwidth_factor*(bbox[i,2]-bbox[i,1])
            end
        end

        retval = 0
        for i in [1:length(myparams);]
            if myparams[i]<bbox[i,1]
                retval += cosh((bbox[i,1]-myparams[i])/wallwidth[i])-1.0
            elseif bbox[i,2] < myparams[i]
                retval += cosh((myparams[i]-bbox[i,2])/wallwidth[i])-1.0                
            end
        end

        return 2*retval
    end

    traj_increment = 100
    params = 0  # Make sure to have this here so that params stays defined beyond the try/catch
    if ( !(typeof(bbox)<:Dict) ); error("Currently only supporting softbox=true, bbox must be a Dict"); end;
    try
        params = copy(inverse_wall(bbox, args, seed))
    catch
        error("Were all initial param values within the indicated walls?")
    end
    eta = start_eta
    trajectory = zeros(2+length(params), traj_increment); cpm_traj = zeros(2, traj_increment)

    if verbose
        @printf "%d: eta=%g ps=" 0 eta 
        print_vector(vector_wrap(bbox, args, params))
        @printf "\n"
    end
    
    if softbox
        if !(typeof(bbox)<:Dict); error("bhm: If softbox=true, then bbox must eb a Dict"); end
        cost, grad, hess = keyword_vgh((;pars...)->func(;wallwrap(bbox, pars)...), args, params)
    elseif hardbox
        cost, grad, hess = keyword_vgh((;pars...) -> func(;pars...), args, params)
    else
        cost, grad, hess = keyword_vgh((;pars...) -> func(;pars...) + wall_cost(args, bbox; wallwidth=wallwidth, pars...),
            args, params)        
    end
        
    chessdelta = zeros(size(params))
    
    i=0  # here so variable i is available outside the loop
    for i in [1:maxiter;]
        if i > size(trajectory, 2)
            trajectory = [trajectory zeros(2+length(params), traj_increment)]
            cpm_traj   = [cpm_traj   zeros(2, traj_increment)]
        end
        trajectory[1:2, i]   = [eta;cost]
        trajectory[3:end, i] = vector_wrap(bbox, args, params)
        
        hessdelta  = - inv(hess)*grad
        try
            if verbose && verbose_level >= 2
                @printf("bhm: about to try cpm with grad : "); print_vector_g(grad); print("\n")
                @printf("bhm:   hess :"); print_vector_g(hess[:]); print("\n");
            end
            if verbose && verbose_level >= 2
                cpm_out = constrained_parabolic_minimization(hess, grad'', eta, 
                    maxiter=500, tol=1e-20, do_plot=true, verbose=true)                
            else
                cpm_out = constrained_parabolic_minimization(hess, grad'', eta, maxiter=500, tol=1e-20)
            end
            chessdelta = cpm_out[1]; cpm_traj[1,i] = cpm_out[5]; cpm_traj[2,i] = cpm_out[6]
            jumptype = "not failed"
        catch y
            jumptype = "failed"
            if verbose
                @printf "Constrained parabolic minimization failed with error %s\n" y
                @printf "\n"
                @printf "eta was %g\n" eta
                @printf "grad was\n"
                print_vector(grad)
                @printf "\n\nhess was\n"
                for k in [1:length(grad);]
                    print_vector(hess[k,:])
                    @printf "\n"
                end
                @printf "\n"
                matwrite("error_report.mat", Dict("grad"=>grad, "hess"=>hess, "eta"=>eta))
            end
            break
        end

        if norm(hessdelta) <= eta
            new_params = params + hessdelta
            jumptype = "Newton"
        elseif jumptype != "failed" 
            new_params = params + chessdelta
            jumptype  = "constrained"
        end

        if jumptype != "failed"
            if softbox
                new_cost, new_grad, new_hess = 
                    keyword_vgh((;pars...) -> func(;wallwrap(bbox, pars)...), args, new_params)
                if verbose && verbose_level >=2
                    @printf("bhm: had new_params = : "); print_vector_g(vector_wrap(bbox, args, params)); print("\n");
                    @printf("bhm: and my bbox was : "); print(bbox); print("\n")
                    @printf("bhm: and my wallwrap output was : "); print(wallwrap(bbox, make_dict(args, new_params))); print("\n")
                    @printf("bhm: and this produced new_grad : "); print_vector_g(new_grad); print("\n")
                    @printf("bhm:   new_hess :"); print_vector_g(new_hess[:]); print("\n");                                        
                end
            elseif hardbox
                for p in [1:length(new_params);]
                    if new_params[p] < bbox[p,1]; new_params[p] = bbox[p,1]; end
                    if bbox[p,2] < new_params[p]; new_params[p] = bbox[p,2]; end
                 end        
                
                new_cost, new_grad, new_hess = keyword_vgh((;pars...) -> func(;pars...), args, new_params)
            else
                new_cost, new_grad, new_hess = keyword_vgh((;pars...) -> func(;pars...) + 
                        wall_cost(args, bbox; wallwidth=wallwidth, pars...),
                    args, new_params)                
            end
            
            if abs(new_cost - cost) < tol || eta < tol
                if verbose
                    @printf("About to break -- tol=%g, new_cost-cost=%g, eta=%g\n", tol, new_cost-cost, eta)
                end
                break
            end
        end

        if jumptype == "failed" || new_cost >= cost  
            if verbose
                @printf("eta going down: new_cost-cost=%g and jumptype='%s'\n", new_cost-cost, jumptype)
                if verbose_level >= 2
                    nwp = vector_wrap(bbox, args, new_params); wp = vector_wrap(bbox, args, params)
                    @printf("   vvv: proposed new params were : "); print_vector_g(nwp); print("\n")
                    @printf("   vvv: proposed delta params was : "); print_vector_g(nwp-wp); print("\n")
                    @printf("   vvv: grad was : "); print_vector_g(grad); print("\n")
                    costheta = dot(new_params-params, grad)/(norm(new_params-params)*norm(grad))
                    @printf("   vvv: costheta of proposed jump was %g\n", costheta)
                end
            end
            eta = eta/2
            costheta = NaN
            if eta < tol
                if verbose
                    @printf("About to break -- tol=%g, new_cost-cost=%g, eta=%g\n", tol, new_cost-cost, eta)
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
                @printf "%d: eta=%g cost=%g jtype=%s costheta=%.3f ps=" i eta cost jumptype costheta
                print_vector_g(vector_wrap(bbox, args, params))
                @printf "\n"
                if verbose_level >= 3
                    @printf "    At this point, grad is ="
                    print_vector_g(grad)
                    @printf "\n"                
                end
            end
        end
    end

    trajectory = trajectory[:,1:i]; cpm_traj = cpm_traj[:,1:i]
    return vector_wrap(bbox, args, params), trajectory, cost, cpm_traj
end


