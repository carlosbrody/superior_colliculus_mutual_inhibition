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

if false
    n1,n2 = iter_examples(examples, results;timestep=40);
    plot_reviews(n1,n2,line_alpha=0.05,dot_alpha=.5);
    plt[:savefig]("timestep_40.png")
    n1,n2 = iter_examples(examples, results;timestep=45);
    plot_reviews(n1,n2,line_alpha=0.05,dot_alpha=.5);
    plt[:savefig]("timestep_45.png")
    n1,n2 = iter_examples(examples, results;timestep=50);
    plot_reviews(n1,n2,line_alpha=0.05,dot_alpha=.5);
    plt[:savefig]("timestep_50.png")
    n1,n2 = iter_examples(examples, results;timestep=55);
    plot_reviews(n1,n2,line_alpha=0.05,dot_alpha=.5);
    plt[:savefig]("timestep_55.png")
    n1,n2 = iter_examples(examples, results;timestep=60);
    plot_reviews(n1,n2,line_alpha=0.05,dot_alpha=.5);
    plt[:savefig]("timestep_60.png")
end

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

function plot_dist(n1,n2)
    pro_unit_choice = [x[2] for x in n1]
    ant_unit_choice = [x[2] for x in n2]
    pro_unit_proant = [x[1] for x in n1]
    ant_unit_proant = [x[1] for x in n2]    
    figure(figsize=(5,4))
    plt[:hist](pro_unit_choice, bins=50,color="black",alpha=.5,label="Pro Unit");
    plt[:hist](ant_unit_choice, bins=50, color="red",alpha=.5,label="Anti Unit");
    plt[:xlim](-7,7)
    plt[:xlabel]("Choice d'",fontsize=16)
    plt[:ylabel]("Count", fontsize=16)
    plt[:xticks](fontsize=12)
    plt[:yticks](fontsize=12)
    plt[:legend]()
    plt[:tight_layout]()
    plt[:axvline](0,color="k",linestyle="--")

    figure(figsize=(5,4))
    plt[:hist](pro_unit_proant, bins=50,color="black",alpha=.5,label="Pro Unit");
    plt[:hist](ant_unit_proant, bins=50, color="red",alpha=.5,label="Anti Unit");
    plt[:xlim](-10,10)
    plt[:xlabel]("Pro/Anti d'",fontsize=16)
    plt[:ylabel]("Count", fontsize=16)
    plt[:xticks](fontsize=12)
    plt[:yticks](fontsize=12)
    plt[:legend]()
    plt[:tight_layout]()
    plt[:axvline](0,color="k",linestyle="--")

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
    # N Solutions, 3 opto, pro/anti, 4 nodes, 61 timesteps, 10 trials

    # Compute pro/anti dprime on the merged pro units
    pro_trials_unit1 = examples[index,1,1,1,timestep,:]
    ant_trials_unit1 = examples[index,1,2,1,timestep,:]
    pro_trials_unit4 = examples[index,1,1,4,timestep,:]
    ant_trials_unit4 = examples[index,1,2,4,timestep,:]
    pro_trials_units14 = vcat(pro_trials_unit1,pro_trials_unit4) # Merge the two pro units is equivalent to simulating both sides of cues 
    ant_trials_units14 = vcat(ant_trials_unit1,ant_trials_unit4) # Mering both sides...
    pa_units14_dprime = calc_dprime(pro_trials_units14,ant_trials_units14)

    # Compute pro/anti dprime on the merge anti units
    pro_trials_unit2 = examples[index,1,1,2,timestep,:]
    ant_trials_unit2 = examples[index,1,2,2,timestep,:]
    pro_trials_unit3 = examples[index,1,1,3,timestep,:]
    ant_trials_unit3 = examples[index,1,2,3,timestep,:]
    pro_trials_units23 = vcat(pro_trials_unit2,pro_trials_unit3)
    ant_trials_units23 = vcat(ant_trials_unit2,ant_trials_unit3)
    pa_units23_dprime = calc_dprime(pro_trials_units23,ant_trials_units23)

    # Compute choice dprime on the merged pro units
    ipsi = examples[index,1,:,1,end,:] .>= examples[index,1,:,4,end,:]
    node1 = examples[index,1,:,1,timestep,:]
    node4 = examples[index,1,:,4,timestep,:]
    ipsi_trials_units14 = vcat(node1[ipsi],node4[.~ipsi])
    cont_trials_units14 = vcat(node1[.~ipsi], node4[ipsi])
    ic_units14_dprime = calc_dprime(ipsi_trials_units14, cont_trials_units14)
    
    # Compute choice dprime on the merged anti units
    node2 = examples[index,1,:,2,timestep,:]
    node3 = examples[index,1,:,3,timestep,:]
    ipsi_trials_units23 = vcat(node2[ipsi],node3[.~ipsi])
    cont_trials_units23 = vcat(node2[.~ipsi], node3[ipsi])
    ic_units23_dprime = calc_dprime(ipsi_trials_units23,cont_trials_units23)
    
    # Return pro/anti and choice dprime for the pro units, and then the anti units
    return [pa_units14_dprime, ic_units14_dprime], [pa_units23_dprime, ic_units23_dprime]
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

