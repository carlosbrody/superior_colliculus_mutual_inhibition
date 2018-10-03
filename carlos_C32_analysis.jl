using PyPlot
using MAT
using TSne

tcost_threshold = -0.0001
fignum          = 100   # figure number on which to plot the t-SNE

if ~isdefined(:examples)
    examples = matread("visualization/MiniC32_C32_examples_50.mat")
end


goods = find(examples["results"]["tcost"] .< tcost_threshold)
files = examples["results"]["files"][goods]
params = examples["results"]["params"][goods,:]


"""
    zscore(A; dims=1)

Given a matrix of ntrials-by-nparams, returns the Z-scored
matrix, i.e., subtracts the mean and divides by the standard
deviation across trials.

# OPTIONAL PARAMETERS

- dims    indicates the dimension along which to z-score. The default
          is dims=1, meaning the the first dimension of the matrix indicates
          the trials over which to score
"""
function zscore(A; dims=1)
    return (A .- mean(A, dims)) ./ max.(std(A, dims), eps())
end

figure(fignum)
clf()
teed = tsne(zscore(params), 2, 0, 10000)
scatter(teed[:,1], teed[:,2])


"""
input(prompt::String="")::String

Read a string from STDIN. The trailing newline is stripped.

The prompt string, if given, is printed to standard output without a
trailing newline before reading input.
"""
function myinput(prompt::String="")::String
    print(prompt)
    return chomp(readline())
end


# str = ""; i=4
# while str != "q"
#     plot_farm(files[i], testruns=50)
#     str = myinput(string(i) * ": RET for next, q to quit, b for back> ")
#     if typeof(parse(str))==Int64
#         i = parse(str)
#     elseif str==""
#         i = i+1
#     elseif str=="b"
#         i = i-1
#     end
# end
