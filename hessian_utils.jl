# DON'T MODIFY THIS FILE -- the source is in file Cost Function Minimization and Hessian Utilities.ipynb. Look there for further documentation and examples of running the code.


include("constrained_parabolic_minimization.jl")

"""
dict = make_dict(argstrings, x, [starting_dict=Dict()] )

Given a list of strings, and a list of values, makes a dictionary of Symbols to values, with the Symbols 
corresponding to each of the strings.  Mostly used to pass arguments as a keyword-value set into a function.
If one of the elements of argstrings is *not* a string, but is instead a 2-long list, the first element of that 
list should be a string, and the second element of that list should be a positive integer. This will be 
interpreted as "don't take only one value, take this number of values and this parameter will be a vector"

# PARAMS:

* argstrings     A list of strings. Each element may also be a two-long list of a string, positive integer, e.g., ["this" 3]

* x              A vector of numeric values. Its length must be such that all the strings in argstrings
                 can take their corresponding element(s), sequentially, from x

* starting_dict  An optional initial dictionary to work with.  Any key in this starting dictionary matching an argstring
                 will be replaced by the new value. Keys not matched will remain.

# RETURNS:

dict             The symbol dictionary.


# EXAMPLES:

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
    # For error diagnostics, check that the length of the param vector specified in args matches the length of x
    nargs = 0
    for i in [1:length(args);]
        if typeof(args[i])==String # if the entry in args is a string, then there's one corresponding scalar entry in x0
            nargs += 1
        else
            nargs += args[i][2]    # otherwise, the entry in args should be a  [varnamestring, nvals] vector, 
            # indicating that the next nvals entries in x0 are all a single vector, belonging to variable
            # with name varnamestring. 
        end
    end
    if nargs != length(x)
        error("Oy! args and x must indicate the same total number of variables!")
    end

    
    # ---- done error-checking, now main function
    
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


using ForwardDiff


# Here we're going to define a closure over x so that when this code runs, it sets the local variable
# x to report ForwardDiff's verison number; then we export the function FDversion, that simply returns
# x.  When we call FDversion(), it simply returns the value of x, stored locally inside the let block.
# So it is extremely fast.
let x 
    try
        x = Pkg.installed("ForwardDiff").major + 0.1*Pkg.installed("ForwardDiff").minor
    catch
        error("Is ForwardDiff really installed???")
    end

    global FDversion

    @doc """
    v = FDversion()

    Return the installed version of ForwardDiff as a floating point, major.minor
    """ function FDversion()
        return x
    end
end


if FDversion() < 0.6
    # --------------------------------------------------------------
    #
    #               FOR FORWARDDIFF < 0.6   (Julia 0.5.2)
    #
    # --------------------------------------------------------------
    
    
    # """
    # We define functions to convert Duals, the variable types used by ForwardDiff, 
    # to Floats. This is useful if we want to print out the value of a variable 
    # (since print doesn't know how to Duals). Note that after being converted to a Float, no
    # differentiation by ForwardDiff can happen!  e.g. after
    #     x = convert(Float64, y)
    # ForwardDiff can still differentiate y, but it can't differentiate x
    # """

    import Base.convert
    convert(::Type{Float64}, x::ForwardDiff.Dual) = Float64(x.value)
    function convert(::Array{Float64}, x::Array{ForwardDiff.Dual}) 
        y = zeros(size(x)); 
        for i in 1:prod(size(x)) 
            y[i] = convert(Float64, x[i]) 
        end
        return y
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
          
else
    
    # --------------------------------------------------------------
    #
    #         FOR FORWARDDIFF >= 0.6   (Julia 0.6 and onwards)
    #
    # --------------------------------------------------------------
    
    
    @doc """
    e = get_eltype(vars)
    
    vars should be a tuple of variables. If any of them is a ForwardDiff Dual, this function returns
    the typeof of that one (the first one encountered); otherwise it returns Float64.
    """ function get_eltype(vars)
        # print("vars is "); print(vars); print("\n")
        # print("typeof(vars) is "); print(typeof(vars)); print("\n")
        if ! (typeof(vars)<:Tuple)
            error("vars must be a Tuple (of variables)")
        end
        for v in vars
            # print("About to check v = "); print(v); print("\n")
            # If it's an array, check its elements
            if typeof(v)<:Array && length(v)> 0 && typeof(v[1])<:ForwardDiff.Dual
                return typeof(v[1])
            end
            # If it's an array of tuples, check each one
            if typeof(v)<:Array && length(v)> 0 && typeof(v[1])<:Tuple
                for tup in v; if typeof(tup[2])<:ForwardDiff.Dual; return typeof(tup[2]); end; end
            end
            # If it's a Pair, check its value
            if typeof(v)<:Pair && typeof(v[2])<:ForwardDiff.Dual
                return typeof(v[2])
            end
            # If it's a tuple, try to turn it into a Dict
            if typeof(v)<:Tuple 
                try; myv = Dict(v); catch; error("Sorry don''t know how to deal with that kind of Tuple"); end
            else
                myv = v
            end
            # If it's a Dict, check all the key contents
            if typeof(myv)<:Dict
                for k in keys(myv)
                    # print("about to check key "); print(k); print(" with value "); print(myv[k]); print("\n")
                    if get_eltype(Tuple(myv[k]))<:ForwardDiff.Dual
                        return get_eltype(Tuple(myv[k]))
                    end
                end
            end
            if typeof(myv)<:ForwardDiff.Dual
                return typeof(myv)
            end
        end
        return Float64
    end


end


# We're deprecating automatic converts from Dual to Float64; 
# instead, let's use an explicit call to a get_value function.
# That means that any accidental casts into Float64s will cause errors, alerting us of 
# problems rather than letting the code run but producing derivatives that are zero
    
    
@doc """
v = get_value(x)
    
If you're going to @print something that might be a ForwardDiff Dual, use this function. It'll
return the value of the number if it is a Dual and just the number if it was not a Dual, suitable for
printing, e.g.
    
    @printf("%g\n", get_value(x))
    
will work regardless of whether x is a ForwardDiff Dual, a Float64, or an Int64    
""" function get_value(x)
    if typeof(x)<:Array && length(x)>0 
        if typeof(x[1])<:ForwardDiff.Dual
            y = zeros(size(x)); 
            if typeof(x[1].value)<:ForwardDiff.Dual  # nested, taking 2nd derivative
                for i=1:length(y); y[i] = x[i].value.value; end;
            else # not nested, 1st derivative
                for i=1:length(y); y[i] = x[i].value; end;
            end
            return y
        elseif typeof(x[1])<:Int64 || typeof(x[1])<:Float64
            return x
        else
            error(@sprintf("Don't know how to get the value of type %s", typeof(x[1])))
        end
    elseif typeof(x)<:Array
        return zeros(size(x))
    elseif typeof(x)<:ForwardDiff.Dual
        if typeof(x.value)<:ForwardDiff.Dual; return x.value.value; else; return x.value; end;
    elseif typeof(x)<:Int64 || typeof(x)<:Float64
        return x
    else
        error(@sprintf("Don't know how to get the value of type %s", typeof(x)))
    end
end    



# DON'T MODIFY THIS FILE -- the source is in file Cost Function Minimization and Hessian Utilities.ipynb. Look there for further documentation and examples of running the code.


using ForwardDiff

if FDversion() < 0.6
    # --------------------------------------------------------------
    #
    #               FOR FORWARDDIFF < 0.6   (Julia 0.5.2)
    #
    # --------------------------------------------------------------



    """
    function value, gradient, hessian = keyword_vgh(func, args, x0)

    Wrapper for vgh() that computes and returns all three of a function's value, gradient, and hessian, but now
    uses make_dict() to apply it to a function that only takes keyword-value pairs. 

    *Note that func MUST also take the keyword parameters nderivs and difforder*. If you declare any vectors or 
    matrices inside func() (or inside any function inside func()), use ForwardDiffZeros with these two parameters, 
    do NOT use zeros(). Your gradients will come out as zero is you use zeros().

    # PARAMETERS

    * func    A function that takes keyword-value pairs only, including nderivs and difforder.  I.e., it must be a function declared as `function func(; nderivs=0, difforder=0, other_kw_value_pairs)` or as `function func(; nderivs=0, difforder=0, other_kw_value_pairs_dict...)`
    * args    A list of strings indicating names of variables to work with
    * x0      A vector with the value of the variables indicates in args.  **See make_dict() for how to pass both scalars and vectors as variables**

    # IMPORTANT JULIA BUG

    If you modify func, it is possible that keyword_vgh() will still work on the previously defined version. AACK!  
    That's horrible! Alice Yoon's tip on the workaround: instead of func(), use (;params...) -> func(; params...) and then
    everything will be fine. Perhaps this bug will be fixed in Julia 0.6

    # EXAMPLE:

    function tester(;a=10, b=20, c=30, nderivs=0, difforder=0)
        M = ForwardDiffZeros(3, 3; nderivs=nderivs, difforder=difforder)
        M[1,1] = a^2*10
        M[2,2] = b*20
        M[3,3] = a*sqrt(c)*30.1
        return trace(M)
    end

    value, grad, hess = keyword_vgh(tester, ["a", "c"], [10, 3.1])

    value, grad, hess = keyword_vgh((;params...) -> tester(;params...), ["a", "c"], [10, 3.1])

    """
    function keyword_vgh(func, args, x0)

        value, gradient, hessian = vgh(x -> func(;nderivs=length(x), difforder=2, make_dict(args, x)...), x0)

        return value, gradient, hessian    
    end
    
    
else    
    # --------------------------------------------------------------
    #
    #         FOR FORWARDDIFF >= 0.6   (Julia 0.6 and onwards)
    #
    # --------------------------------------------------------------


    """
    function value, gradient, hessian = keyword_vgh(func, args, x0)

    Wrapper for vgh() that computes and returns all three of a function's value, gradient, and hessian, but now
    uses make_dict() to apply it to a function that only takes keyword-value pairs. 

    *Note that if you declare any vectors or matrices inside func() (or inside any function inside func()), 
    you will need to make sure they are ForwardDiff Duals if you want to differentiate w.r.t. their contents.*
    The function get_eltype()  can help you with this.  For example, if you have three variables, a, b, and c,
    and you don't know in advance which one you will differentiate w.r.t., you could declare a matrix using
    
        new_matrix = zeros(get_eltype((a,b,c)), 2, 3)
    
    and that will make sure that it is the right type if you later want derivatives of contents of that 
    matrix.

    # PARAMETERS

    * func    A function that takes keyword-value pairs only, including nderivs and difforder.  I.e., it must be a function declared as `function func(; nderivs=0, difforder=0, other_kw_value_pairs)` or as `function func(; nderivs=0, difforder=0, other_kw_value_pairs_dict...)`
    * args    A list of strings indicating names of variables to work with
    * x0      A vector with the value of the variables indicates in args.  **See make_dict() for how to pass both scalars and vectors as variables**

    # IMPORTANT JULIA BUG

    If you modify func, it is possible that keyword_vgh() will still work on the previously defined version. AACK!  
    That's horrible! Alice Yoon's tip on the workaround: instead of func(), use (;params...) -> func(; params...) and then
    everything will be fine. Not sure yet whether this bug is fixed in Julia 0.6

    # EXAMPLE:

    function tester(;a=10, b=20, c=30)
        M = zeros(get_eltype((a,b,c)), 3, 3)
        M[1,1] = a^2*10
        M[2,2] = b*20
        M[3,3] = a*sqrt(c)*30.1
        return trace(M)
    end

    value, grad, hess = keyword_vgh(tester, ["a", "c"], [10, 3.1])

    value, grad, hess = keyword_vgh((;params...) -> tester(;params...), ["a", "c"], [10, 3.1])

    """
    function keyword_vgh(func, args, x0)

        value, gradient, hessian = vgh(x -> func(;make_dict(args, x)...), x0)

        return value, gradient, hessian    
    end
    
end


# DON'T MODIFY THIS FILE -- the source is in file Cost Function Minimization and Hessian Utilities.ipynb. Look there for further documentation and examples of running the code.



##############################################################################################
#                                                                                            #
#      OLD AND LARGELY UNUSED AND UNMAINTAINED KEYWORD GRADIENTS AND HESSIANS                #
#                                                                                            #
##############################################################################################


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


# DON'T MODIFY THIS FILE -- the source is in file Cost Function Minimization and Hessian Utilities.ipynb. Look there for further documentation and examples of running the code.


"""
function constrained_Hessian_minimization(seed, func; start_eta=10, tol=1e-6, maxiter=400,
    verbose=false)

LARGELY SUPERSEDED BY BBOX_HESSIAN_KEYWORD_MINIMIZATION

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


# DON'T MODIFY THIS FILE -- the source is in file Cost Function Minimization and Hessian Utilities.ipynb. Look there for further documentation and examples of running the code.




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
  


# DON'T MODIFY THIS FILE -- the source is in file Cost Function Minimization and Hessian Utilities.ipynb. Look there for further documentation and examples of running the code.



# Julia 0.5 returns an error if you run Pkg.installed on a package that has not been installed
try
    Pkg.installed("JLD")
catch
    Pkg.add("JLD")
end

# Julia 0.6 does not crash but returns a Void if you run Pkg.installed on a package that has not been installed
if isa(Pkg.installed("JLD"), Void)
    Pkg.add("JLD")
end


using JLD
using MAT


######################################################
#                                                    #
#         BBOX_HESSIAN_KEYWORD_MINIMIZATION          #
#                                                    #
######################################################




"""
function bbox_Hessian_keyword_minimization(seed, args, bbox, func; start_eta=10, tol=1e-6, 
maxiter=400, verbose=false, report_file="")

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
            take nderivs=0, difforder=0  and declare any new matrices using ForwardDiffZeros() instead of zeros().
            func must either return a scalar, or the first output it returns must be a scalar. 
            That scalar is what will be minimized. The trajectory across the minimization of 
            any further outputs that f() returns will be available in ftraj (see RETURNS below)


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
- cpm_traj     A 2-by-nsteps matrix, containing reports from the contrained parabolic minimization at each timestep.
             The first row is niters (how many iterations cpm's 1-d minimization ran for) and the second row is
             Dlambda, the last change in the parameter being minimized in cpm's internal search
- ftraj     Further components for the trajectory, will be an Array{Any}(3, nsteps). First row is gradient,
            second row is Hessian, third row is second-and-further outputs of func, each one at each step of
            the minimization. **NOTE** that if these further outputs contain variables that are being minimized, 
            they'll come out as ForwardDiff Duals, which you might not want!  So, for example, you might want to
            convert vectors and matrices into Float64s before returning them in those extra outputs. E.g.,
            if you want to return sum(err.*err) as the scalar to be minimized, and also return err, in your 
            cost function you would write   " return sum(err.*err), Array{Float64}(err) ".   That way the first,
            scalar output can still be differentiated, for minimization, and the second one comes out in readable form.



# EXAMPLE:  (see also a more complete example in Cost Function Minimization and Hessian Utilities.ipynb)

```
function tester(;x=5, y=10, z=20, nderivs=0, difforder=0)
    return x^2*y + z/tanh(y)
end

params, trajectory = bbox_Hessian_keyword_minimization([0.5, 0.5], ["x", "y"], [1.1 2 ; 1.1 4], tester, 
    verbose=true, tol=1e-12, start_eta=1);
```


"""
function bbox_Hessian_keyword_minimization(seed, args, bbox, func; start_eta=0.1, tol=1e-6, maxiter=400,
    verbose=false, verbose_level=1, verbose_every=1, softbox=true, hardbox=false, report_file="")

    # --- check that saving will be done to a .jld file ---
    if length(report_file)>0 && splitext(report_file)[2] != ".jld"
        if splitext(report_file)[2] == ""
            report_file = resport_file * ".jld"
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
    trajectory = zeros(2+length(params), traj_increment); cpm_traj = zeros(2, traj_increment)
    
    ftraj = Array{Any}(3,0)  # will hold gradient, hessian, and further_out,  per iteration

    further_out =[];  # We define this variable here so it will be available for stashing further outputs from func
    
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
        @printf "%d: eta=%g ps=" 0 eta 
        print_vector(vector_wrap(bbox, args, params))
        @printf "\n"
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
    
    i=0  # here so variable i is available outside the loop
    for i in [1:maxiter;]
        if i > size(trajectory, 2)
            trajectory = [trajectory zeros(2+length(params), traj_increment)]
            cpm_traj   = [cpm_traj   zeros(2, traj_increment)]
        end
        trajectory[1:2, i]   = [eta;cost]
        trajectory[3:end, i] = vector_wrap(bbox, args, params)
        ftraj = [ftraj [grad, hess, further_out]]

        if length(report_file)>0
            save(report_file, Dict("traj"=>trajectory[:,1:i], "cpm_traj"=>cpm_traj[:,1:i], "ftraj"=>ftraj))
        end
        
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
            new_cost, new_grad, new_hess = keyword_vgh(internal_func, args, new_params)   # further_out may mutate
            if verbose && verbose_level >=2
                @printf("bhm: had new_params = : "); print_vector_g(vector_wrap(bbox, args, params)); print("\n");
                @printf("bhm: and my bbox was : "); print(bbox); print("\n")
                @printf("bhm: and my wallwrap output was : "); print(wallwrap(bbox, make_dict(args, new_params))); print("\n")
                @printf("bhm: and this produced new_grad : "); print_vector_g(new_grad); print("\n")
                @printf("bhm:   new_hess :"); print_vector_g(new_hess[:]); print("\n");                                        
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
    if length(report_file)>0
        save(report_file, Dict("traj"=>trajectory, "cpm_traj"=>cpm_traj, "ftraj"=>ftraj))
    end
    
    return vector_wrap(bbox, args, params), trajectory, cost, cpm_traj, ftraj
end


