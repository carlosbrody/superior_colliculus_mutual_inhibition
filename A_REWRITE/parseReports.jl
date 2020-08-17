##  === Pulling in Reports from NEW-style training and parsing them
#
# Produces a quick report on how many networks are below cost threshold
# and prints out some of their parameters (no figures)

if !in(".", LOAD_PATH)
    push!(LOAD_PATH, ".")
end

using DelimitedFiles
include("sixNodeSetup_C32.jl")
args = argsSeedBounder()[1]   # get the string list of model params



pullReports=false; if pullReports   # if there are new runs on the VMs, get them
    remotes = "spock"

    if remotes == "googleVMs"
        machines = ["proanti006", "proanti007", "proanti008",
                    "proanti009", "proanti010", "proanti011"]

        for i in machines
            run(`./pull_reports.sh $i`)
        end
    elseif remotes == "spock"
        run(`./pull_reports.sh`)
    else
        error("Don't know remotes $remotes")
    end

    # At the of this block, report files will have been put in ../../Reports
end
##
global costs  = fill(0.0, 0)   # the cost for each solution
global fnames = fill("",  0)   # the filename for each solution
global pvals  = fill(0.0, 0, length(args))  # model params for each solution
global srands = fill(0, 0)     # random seeds for each solution

# Now we'll parse each Report file to find where it ended up at.
u = readdir("../../Reports")
for i in filter(x -> startswith(x, "r6_spockRun3"), u)
    # Read the file into string s:
    fp = open("../../Reports/$i", "r")
    s = read(fp, String)
    close(fp)
    println(i)
    # Find the last cost report:
    z = findlast("OVERALL", s)
    if z != nothing  && findlast("random seed", s)!=nothing
        z2 = findfirst("cost=", s[z[end]:end])
        z3 = findfirst(",", s[z[end]:end])
        mycost = parse(Float64, s[z[end]+z2[end]:z[end]+z3[1]-2])
        # println("$i: mycost=$mycost")
        global costs  = vcat(costs, mycost)
        global fnames = vcat(fnames, i)

        z = findlast("ps=[", s)[end].+1
        n = findfirst("]", s[z+1:end])[1] + z-1
        global pvals = vcat(pvals, readdlm(Vector{UInt8}(s[z:n]), ','))

        sz1 = findlast("random seed", s)
        sz2 = findfirst('\n', s[sz1[end]:end])
        global srands = vcat(srands, parse(Int64, s[sz1[end]+2:sz1[end]+sz2-2]))
    end
    # println("Checked $i")
end

println("$(length(costs)) networks went into further training")

nguys = 20   # number of solutions with the best costs that we'll look at
colids = [5, 7, 21, 10, 22, 4, 18, 3] # args[colids] will be the weights collected into "merges"

# we only want to see the best cost solutions
p = sortperm(costs, rev=true)
costs = costs[p]
fnames = fnames[p]
pvals = pvals[p,:]
srands = srands[p]
display([fnames costs])

# merges will hold the chosen weights
merges = [reshape(args[colids],1,length(colids)) ; pvals[end-nguys:end, colids]]
merges = hcat(["cost" ; costs[end-nguys:end]], merges)
display(merges)
