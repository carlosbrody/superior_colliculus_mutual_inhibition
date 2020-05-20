module ResultsAnalysis

using DelimitedFiles

export readCostsFile, writeCostsFile

topline = [
"seedrand"
 "cost"
 "anti_rule_strength"
 "right_light_excitation"
 "vW_AP"
 "sigma"
 "hW_A"
 "sW_P"
 "const_pro_bias"
 "opto_strength"
 "sW_A"
 "dW_PA"
 "hW_P"
 "target_period_excitation"
 "dW_AP"
 "pro_rule_strength"
 "constant_excitation"
 "vW_PA"
]


function  readCostsFile(;fname="negCosts_proanti003.csv", maxCost::Float64=Inf, others...)
   if typeof(fname)<:Array{String}
      G = readCostsFile(fname=fname[1]; maxCost=maxCost, others...)
      for k=2:length(fname)
         G = vcat(G, readCostsFile(fname=fname[k]; maxCost=maxCost, others...)[2:end,:])
      end
      return G
   end

   G = readdlm(fname, ',')

   if G[1,4] == ""
      A = Array{Any}(undef, 1 + Int64(size(G,1)/2), 18)
      A[1,:] .= topline

      rownum = 2;
      for i=1:2:size(G,1)
         A[rownum, 1:2]   .= G[i,1:2]
         A[rownum, 3:end] .= Array{Float64}(G[i+1,1:16])
         rownum += 1
      end
   else
      A = Array{Any}(undef, 1 + size(G,1), 18)
      A[1,:] .= topline
      A[2:end,:] = G
   end

   A = A[[1;findall(A[2:end,2] .<= maxCost) .+ 1],:]
   for k in keys(others)
      mycol = findfirst(A[1,:].==String(k))
      @assert mycol != nothing  "couldn't find column $k"
      mylims = others[k]
      A = A[[1;findall(A[2:end,mycol] .>= mylims[1]) .+ 1],:]
      A = A[[1;findall(A[2:end,mycol] .<= mylims[2]) .+ 1],:]
   end

   return A
end


function writeCostsFile(fname::String, M)
   fp = open(fname, "w")
   # topline = Array{String,2}(undef, 1,18)
   # topline[:] .= M[1,:]
   # writedlm(fp, topline, ',')
   for i=2:size(M,1)
      print(fp, Int64(M[i,1]), ", ")
      writedlm(fp, Array{Float64}(M[i,2:end])', ',')
   end
   close(fp)
end

end # === END MODULE
