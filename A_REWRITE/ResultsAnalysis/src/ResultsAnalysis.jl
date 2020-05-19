module ResultsAnalysis

using DelimitedFiles

export readCostsFile

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


function  readCostsFile(;fname="negCosts_proanti003.csv", maxCost=Inf)

   G = readdlm(fname, ',')

   A = Array{Any}(undef, 1 + Int64(size(G,1)/2), 18)
   A[1,:] .= topline

   rownum = 2;
   for i=1:2:size(G,1)
      A[rownum, 1:2]   .= G[i,1:2]
      A[rownum, 3:end] .= Array{Float64}(G[i+1,1:16])
      rownum += 1
   end

   return A[[1;findall(A[2:end,2] .<= maxCost) .+ 1],:]
end


end # === END MODULE
