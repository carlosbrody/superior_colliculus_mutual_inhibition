include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")
include("unilateral_analysis.jl")
include("cluster_farms.jl")
include("cluster_example.jl")

#cluster_example_trajectories("C32", "MiniC32"; threshold=-0.0001) # This generates a bunch of examples

# Need to filter
examples,results = load("MiniC32_C32_examples_feb.jld","examples","results");
# N Solutions, 3 opto, pro/anti, 4 nodes, 61 timesteps, 10 trials
# What time window? 50 - 61 is choice period, but lets try 50:55

# n1,n2 = iter_examples(examples, results);
# plot_reviews(n1,n2,line_alpha=0.05,dot_alpha=.5);
# plt[:savefig]

function iter_examples(examples,results;timestep =55,threshold=10,line_alpha=.1,dot_alpha=.2,plot_line=true)
    # Then computes the dprime for each solution
    # Then plots each cell

    badcost = results["cost"] .>=threshold # should already be filtered
    examples = examples[.~badcost,:,:,:,:,:]
    n1 = []
    n2 = []
    
    for i=1:size(examples)[1]
        this_n1,this_n2 = calc_example2(examples,i;timestep=timestep)       
        push!(n1,this_n1)
        push!(n2,this_n2)
    end
    return n1,n2
end

function plot_reviews(all_n1,all_n2; line_alpha=.1,dot_alpha=.2,plot_line=true)
    figure(figsize=(5,5))

    for i=1:size(all_n1)[1]
        n1 = all_n1[i]
        n2 = all_n2[i]
        if plot_line
            plot([n1[1],n2[1]],[n1[2],n2[2]],"k-",alpha=line_alpha)
        end
        if i ==1
            plot(n1[1],n1[2],"ko",alpha=dot_alpha,label="Pro")
            plot(n2[1],n2[2],"ro",alpha=dot_alpha,label="Anti")
        else 
            plot(n1[1],n1[2],"ko",alpha=dot_alpha)
            plot(n2[1],n2[2],"ro",alpha=dot_alpha)
        end
    end
    xlim([-5,5])
    ylim([-5,5])
    plot([-10,10],[0,0],"k--")
    plot([0,0],[-10,10],"k--")
    ylabel("Choice d'",fontsize=16)
    xlabel("Pro/Anti d'",fontsize=16)
    
    ax = gca()
    ax[:set_aspect]("equal")
    plt[:xticks](fontsize=14)
    plt[:yticks](fontsize=14)
    plt[:tight_layout]()
end


function calc_example2(examples, index;timestep=55)
    # This verison merges cells across hemispheres because we only simulate left trials

    # Compute pro/anti dprime on the merged pro units
    p1 = examples[index,1,1,1,timestep,:]
    a1 = examples[index,1,2,1,timestep,:]
    p4 = examples[index,1,1,4,timestep,:]
    a4 = examples[index,1,2,4,timestep,:]
    p14 = vcat(p1,p4) # Merge the two pro units is equivalent to simulating both sides of cues 
    a14 = vcat(a1,a4) # Mering both sides...
    pa14_dprime = calc_dprime(p14,a14)

    # Compute pro/anti dprime on the merge anti units
    p2 = examples[index,1,1,2,timestep,:]
    a2 = examples[index,1,2,2,timestep,:]
    p3 = examples[index,1,1,3,timestep,:]
    a3 = examples[index,1,2,3,timestep,:]
    p23 = vcat(p2,p3)
    a23 = vcat(a2,a3)
    pa23_dprime = calc_dprime(p23,a23)

    # Compute choice dprime on the merged pro units
    ipsi = examples[index,1,:,1,end,:] .>= examples[index,1,:,4,end,:]
    node1 = examples[index,1,:,1,timestep,:]
    node4 = examples[index,1,:,4,timestep,:]
    pro_ipsi = vcat(node1[ipsi],node4[.~ipsi])
    pro_contra = vcat(node1[.~ipsi], node4[ipsi])
    ic14_dprime = calc_dprime(pro_ipsi,pro_contra)
    
    # Compute choice dprime on the merged anti units
    node2 = examples[index,1,:,2,timestep,:]
    node3 = examples[index,1,:,3,timestep,:]
    anti_ipsi = vcat(node2[ipsi],node3[.~ipsi])
    anti_contra = vcat(node2[.~ipsi], node3[ipsi])
    ic23_dprime = calc_dprime(anti_ipsi,anti_contra)
    
    # Return pro/anti and choice dprime for the pro units, and then the anti units
    return [pa14_dprime, ic14_dprime], [pa23_dprime, ic23_dprime]
end
 

function calc_example(examples, index;timestep=55)
    # This version does not merge cells across the hemisphere, and thus suffers from a selection bias of only simulating left trials
    
    #p1 = examples[index,1,1,1,timestep,:]
    #a1 = examples[index,1,2,1,timestep,:]
    #pa1_dprime = calc_dprime(p1,a1)
    #p2 = examples[index,1,1,2,timestep,:]
    #a2 = examples[index,1,2,2,timestep,:]
    #pa2_dprime = calc_dprime(p2,a2)
    #p3 = examples[index,1,1,3,timestep,:]
    #a3 = examples[index,1,2,3,timestep,:]
    #pa3_dprime = calc_dprime(p3,a3)
    #p4 = examples[index,1,1,4,timestep,:]
    #a4 = examples[index,1,2,4,timestep,:]
    #pa4_dprime = calc_dprime(p4,a4)

    #ipsi = examples[index,1,:,1,end,:] .>= examples[index,1,:,4,end,:]
    #temp1 = examples[index,1,:,1,timestep,:]
    #ic1_dprime = calc_dprime(temp1[ipsi],temp1[.~ipsi])
    #temp2 = examples[index,1,:,2,timestep,:]
    #ic2_dprime = calc_dprime(temp2[ipsi],temp2[.~ipsi])
    ## Gotta flip ipsi/contra for this two
    #temp3 = examples[index,1,:,3,timestep,:]
    #ic3_dprime = calc_dprime(temp3[.~ipsi],temp3[ipsi])
    #temp4 = examples[index,1,:,4,timestep,:]
    #ic4_dprime = calc_dprime(temp4[.~ipsi],temp4[ipsi])
    #return [pa1_dprime, ic1_dprime], [pa2_dprime, ic2_dprime], [pa3_dprime, ic3_dprime], [pa4_dprime, ic4_dprime]
end
 
function calc_dprime(vals1,vals2)
    dprime = (mean(vals1) - mean(vals2))./sqrt(.5*(var(vals1) + var(vals2)))
    return dprime
end

