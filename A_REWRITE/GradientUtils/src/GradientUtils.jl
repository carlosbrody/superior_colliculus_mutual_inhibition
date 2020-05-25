module GradientUtils

using ForwardDiff
using Debugger

export make_dict, get_eltype, get_value, vgh, keyword_vgh


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


@doc """
e = get_eltype(vars)

vars should be a tuple of variables. If any of them is a ForwardDiff Dual, this function returns
the typeof of that one (the first one encountered); otherwise it returns Float64.
""" function get_eltype(vars::Tuple)
    for v in vars
        # for arrays or singletons:
        if eltype(v)<:ForwardDiff.Dual
            return eltype(v)
        end
        # If it's an array of tuples, check each one
        if eltype(v)<:Tuple
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
        if typeof(myv)<:AbstractDict
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

function get_eltype(var)
    return get_eltype(Tuple(var))
end


@doc """
v = get_value(x)

If you're going to print something that might be a ForwardDiff Dual, use this
function. It'll return the value of the number if it is a Dual and just the
number if it was not a Dual, suitable for printing, e.g.

    println(get_value(x))

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
            error("Don't know how to get the value of type ", typeof(x[1]))
        end
    elseif typeof(x)<:Array
        return zeros(size(x))
    elseif typeof(x)<:ForwardDiff.Dual
        if typeof(x.value)<:ForwardDiff.Dual; return x.value.value; else; return x.value; end;
    elseif typeof(x)<:Int64 || typeof(x)<:Float64
        return x
    else
        error("Don't know how to get the value of type ", typeof(x))
    end
end




using DiffResults

@doc """
function value, gradient, hessian = vgh(func, x0)

Wrapper for ForwardDiff.hessian!() that computes and returns all three of a function's value, gradient, and hessian.

EXAMPLE:
========

function tester(x::Vector)

    return sum(x.*x)
end

value, grad, hess = vgh(tester, [10, 3.1])
""" function vgh(func, x0)
    out = DiffResults.HessianResult(x0)
    out = ForwardDiff.hessian!(out, func, x0)
    value    = DiffResults.value(out)
    gradient = DiffResults.gradient(out)
    hessian  = DiffResults.hessian(out)

    return value, gradient, hessian
end


@doc """
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

""" function keyword_vgh(func, args, x0)

    value, gradient, hessian = vgh(x -> func(;make_dict(args, x)...), x0)

    return value, gradient, hessian
end





end # END MODULE
