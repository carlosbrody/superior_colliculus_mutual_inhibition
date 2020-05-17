"""
Takes an R^N->R function, over an unbounded domain, and applies bounds to it.
"""
module TanhWalls

export boundedVersion

"""
   boundedVersion(func, bounds)

   given a function that takes a vector length N as input, and a bounds array
   size N-by-2 where the first column is a lower bound, and the second colum is
   an upper bound, returns a new function, of a new vector, that is a one-to-one
   version of the original, but the new<-->old variable mapping is such that the
   new vars are unbounded and the old vars are bounded.

   = RETURNS:

   -  newfunc     function of new var values, such that
                     newfunc(old2new(x)) == func(x)
                     func(new2old(y))    == newfunc(y)

   -  old2new     a function that takes a vector of the old var values and maps
                  them to new (unbounded) var values

   -  new2old     a function that takes a vector of the new (unbounded) var
                  values and maps them to old (bounded) var values

"""
function boundedVersion(func, bounds)

   # d = 0.5*(upper-lower)
   # m = 0.5*(upper+lower)
   # old = m .+ d.*tanh.((new.-m)./d)
   #
   # new = m .+ d.*atanh.((old .- m)./d)

   M = 0.5 .* sum( bounds, dims=2)[:]
   D = 0.5 .* diff(bounds, dims=2)[:]

   function old2new(x)
      @assert (length(size(x))==1) && (length(x)==size(bounds,1)) "x must be a vector of length $(size(bounds,1))"
      @assert all(x .>= bounds[:,1]) && all(x .<= bounds[:,2])  "$x must be within $bounds"
      return atanh.((x .- M)./D).* D .+ M
   end

   function new2old(x)
      @assert (length(size(x))==1) && (length(x)==size(bounds,1)) "x must be a vector of length $(size(bounds,1))"
      return tanh.((x .- M)./D).* D .+ M
   end

   function newfunc(y)
      return func(new2old(y))
   end

   return newfunc, old2new, new2old

end





end
