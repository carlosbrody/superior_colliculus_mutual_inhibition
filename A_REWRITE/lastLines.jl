using DelimitedFiles

if length(ARGS)>0
   for i=1:length(ARGS)
      println("--- $i: ", ARGS[i])
      A = readdlm(ARGS[i], '\n')
      println(A[end-1])
      println(A[end])
   end
end

##

# fp = fopen("commonSetup.jl", "r")
# A = readdlm("commonSetup.jl", '\n')
