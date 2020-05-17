using Random
# using ReverseDiff
using ForwardDiff

function fun3(x)
   Random.seed!(10)
   input = ones(get_eltype(x), size(x))

   input[1] = input[1]*x[1]
   println(eltype(x))
   return sum(input)
end



fun3([1,2,3])

println("F says: ", ForwardDiff.gradient(fun3, [1,2,3]))


# println("R says: ", ReverseDiff.gradient(fun3, [1,2,3.0]))
