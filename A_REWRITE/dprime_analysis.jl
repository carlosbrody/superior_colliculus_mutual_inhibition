include("resultsAnalysis.jl")

# TODO
# Only include hit trials?
# did 4-node include only hit trials?
#
# Do I need to balance the number of pro and anti trials in the choice dprime?
# Did I do that in the 4-node network?

NUM_SOLUTIONS = 36
FILENAME = "node6_solutions.jld"

function save_solutions(;num_trials=10)
    ###
     #   Generates a file which is an array of the outputs of runSolution
     #   Only includes good solutions, so I dont need to filter 
    ###
    all_d = Any[]
    for i=1:NUM_SOLUTIONS
        d = runSolution(i,num_trials)
        push!(all_d,d)
    end
    save(FILENAME, Dict("all_d"=>all_d, "NUM_SOLUTIONS"=>NUM_SOLUTIONS))
end

function plot_dprime(;timestep=55,opto_condition=1,time_condition=1)
    
    # Grab the dprime values
    pro_pa_d, ant1_pa_d, ant2_pa_d, pro_choice_d, ant1_choice_d, ant2_choice_d = get_dprime(timestep=timestep, opto_condition=opto_condition, time_condition=time_condition)

    # Make a histogram of just choice encoding
    figure(figsize=(5,4))     
    plt[:hist](pro_choice_d, bins=10, color="k", alpha=.2, label="Pro")
    plt[:hist](ant1_choice_d, bins=10, color="r", alpha=0.2, label="Anti-1")
    plt[:hist](ant2_choice_d, bins=10, color="m", alpha=0.2, label="Anti-2")
    axvline(0, color="k", linestyle="--")
    plt[:xlabel]("Choice d'")
    plt[:ylabel]("Count")
    plt[:title]("Timestep: "*string(timestep))
    plt[:savefig]("choice_dist_timestep_"*string(timestep)*".png")
    
    # Make a histogram of just rule encoding
    figure(figsize=(5,4))     
    plt[:hist](pro_pa_d, bins=10, color="k", alpha=.2, label="Pro")
    plt[:hist](ant1_pa_d, bins=10, color="r", alpha=0.2, label="Anti-1")
    plt[:hist](ant2_pa_d, bins=10, color="m", alpha=0.2, label="Anti-2")
    axvline(0, color="k", linestyle="--")
    plt[:xlabel]("Pro/Anti d'")
    plt[:ylabel]("Count")
    plt[:title]("Timestep: "*string(timestep))
    plt[:savefig]("proanti_dist_timestep_"*string(timestep)*".png")

    # Make the scatter plot figure
    figure(figsize=(5,4))
    plot(pro_pa_d, pro_choice_d, "ko",label="Pro")
    plot(ant1_pa_d, ant1_choice_d, "ro",label="Anti-1")
    plot(ant2_pa_d, ant2_choice_d, "mo",label="Anti-2")
    axvline(0, color="k", linestyle="--")
    axhline(0, color="k", linestyle="--")
    plt[:xlabel]("Pro/Anti d'")
    plt[:ylabel]("Choice d'")
    plt[:title]("Timestep: "*string(timestep))
    plt[:savefig]("scatter_timestep_"*string(timestep)*".png")
end

function get_dprime(; timestep=55,opto_condition=1,time_condition=1)
    ###
     #   opto_condition should be 1=control, 2=delay_period, 3=target_period
     #   time_condition should be 1:4, for the different combinations of delay and target periods
    ###
    all_d, NUM_SOLUTIONS = load(FILENAME, "all_d", "NUM_SOLUTIONS")
    # Need to iterate over solutions, compute dprime values
    pro_pa_d = []
    ant1_pa_d = []   
    ant2_pa_d = []   
    pro_choice_d = []
    ant1_choice_d = []   
    ant2_choice_d = []   
    for i=1:size(all_d)[1]
        # Extract the pro trial traces, and anti trial traces
        # These will be 6 nodes X 61 timesteps, x 10 trials
        pro_traces = all_d[i]["proValls"][opto_condition, time_condition]
        anti_traces = all_d[i]["antiValls"][opto_condition, time_condition]

        # Compute the Pro/Anti Trial d
        pro_units_pa_d, ant1_units_pa_d, ant2_units_pa_d= compute_pa_dprime(pro_traces, anti_traces, all_d[i],timestep=timestep)
        push!(pro_pa_d, pro_units_pa_d)
        push!(ant1_pa_d, ant1_units_pa_d)
        push!(ant2_pa_d, ant2_units_pa_d)

        # Compute the choice d
        pro_units_choice_d, ant1_units_choice_d, ant2_units_choice_d  = compute_choice_dprime(pro_traces, anti_traces, all_d[i],timestep=timestep)
        push!(pro_choice_d, pro_units_choice_d)
        push!(ant1_choice_d, ant1_units_choice_d)
        push!(ant2_choice_d, ant2_units_choice_d)
     end
    return pro_pa_d, ant1_pa_d, ant2_pa_d, pro_choice_d, ant1_choice_d, ant2_choice_d
end 

function calc_dprime(vals1,vals2)
    dprime = (mean(vals1) - mean(vals2))./sqrt(.5*(var(vals1) + var(vals2)))
    return dprime
end

function compute_pa_dprime(pro_traces, anti_traces,d; timestep=55)
    # Extract just the timestep of interest
    pro_trials_pro_units = pro_traces[d["ProNodeID"],timestep,:]
    pro_trials_ant_units = pro_traces[d["AntiNodeID"],timestep,:]
    ant_trials_pro_units = anti_traces[d["ProNodeID"],timestep,:]
    ant_trials_ant_units = anti_traces[d["AntiNodeID"],timestep,:]

    # concatenate across L/R symmetry
    pro_trials_pro_units_vec  = vcat(pro_trials_pro_units[1,:],pro_trials_pro_units[2,:])
    ant_trials_pro_units_vec  = vcat(ant_trials_pro_units[1,:],ant_trials_pro_units[2,:])
    pro_trials_ant1_units_vec = vcat(pro_trials_ant_units[1,:],pro_trials_ant_units[4,:])
    ant_trials_ant1_units_vec = vcat(ant_trials_ant_units[1,:],ant_trials_ant_units[4,:])
    pro_trials_ant2_units_vec = vcat(pro_trials_ant_units[2,:],pro_trials_ant_units[3,:])
    ant_trials_ant2_units_vec = vcat(ant_trials_ant_units[2,:],ant_trials_ant_units[3,:])

    # Compute dprime
    pro_units_pa_dprime  = calc_dprime(pro_trials_pro_units_vec, ant_trials_pro_units_vec)
    ant1_units_pa_dprime = calc_dprime(pro_trials_ant1_units_vec, ant_trials_ant2_units_vec) 
    ant2_units_pa_dprime = calc_dprime(pro_trials_ant2_units_vec, ant_trials_ant2_units_vec) 

    return pro_units_pa_dprime, ant1_units_pa_dprime, ant2_units_pa_dprime
end


function compute_choice_dprime(pro_traces, anti_traces,d; timestep=55)
    # Extract just the timestep of interest
    pro_trials_pro_units = pro_traces[d["ProNodeID"],timestep,:]
    pro_trials_ant_units = pro_traces[d["AntiNodeID"],timestep,:]
    ant_trials_pro_units = anti_traces[d["ProNodeID"],timestep,:]
    ant_trials_ant_units = anti_traces[d["AntiNodeID"],timestep,:]

    # Extract choice
    pro_trials_ipsi_choice = pro_traces[d["ProNodeID"],end,:][1,:] .> pro_traces[d["ProNodeID"],end,:][2,:]
    ant_trials_ipsi_choice = anti_traces[d["ProNodeID"],end,:][1,:] .> anti_traces[d["ProNodeID"],end,:][2,:]   

    # Merge pro units across pro/anti trials into ipsi/contra choice   
    pro_trials_pro_units_ipsi_choice = vcat(pro_trials_pro_units[1,pro_trials_ipsi_choice],pro_trials_pro_units[2,.~pro_trials_ipsi_choice])
    ant_trials_pro_units_ipsi_choice = vcat(ant_trials_pro_units[1,ant_trials_ipsi_choice],ant_trials_pro_units[2,.~ant_trials_ipsi_choice])
    pro_units_ipsi_choice = vcat(pro_trials_pro_units_ipsi_choice,ant_trials_pro_units_ipsi_choice )

    pro_trials_pro_units_cont_choice = vcat(pro_trials_pro_units[1,.~pro_trials_ipsi_choice],pro_trials_pro_units[2,pro_trials_ipsi_choice])
    ant_trials_pro_units_cont_choice = vcat(ant_trials_pro_units[1,.~ant_trials_ipsi_choice],ant_trials_pro_units[2,ant_trials_ipsi_choice])
    pro_units_cont_choice = vcat(pro_trials_pro_units_cont_choice,ant_trials_pro_units_cont_choice )

    # Merge anti units across pro/anti trials into ipsi/contra choice 
    pro_trials_ant1_units_ipsi_choice = vcat(pro_trials_ant_units[1,pro_trials_ipsi_choice],pro_trials_ant_units[4,.~pro_trials_ipsi_choice])
    ant_trials_ant1_units_ipsi_choice = vcat(ant_trials_ant_units[1,ant_trials_ipsi_choice],ant_trials_ant_units[4,.~ant_trials_ipsi_choice])
    ant1_units_ipsi_choice = vcat(pro_trials_ant1_units_ipsi_choice,ant_trials_ant1_units_ipsi_choice )

    pro_trials_ant1_units_cont_choice = vcat(pro_trials_ant_units[1,.~pro_trials_ipsi_choice],pro_trials_ant_units[4,pro_trials_ipsi_choice])
    ant_trials_ant1_units_cont_choice = vcat(ant_trials_ant_units[1,.~ant_trials_ipsi_choice],ant_trials_ant_units[4,ant_trials_ipsi_choice])
    ant1_units_cont_choice = vcat(pro_trials_ant1_units_cont_choice,ant_trials_ant1_units_cont_choice )

    # Merge anti units across pro/anti trials into ipsi/contra choice 
    pro_trials_ant2_units_ipsi_choice = vcat(pro_trials_ant_units[2,pro_trials_ipsi_choice],pro_trials_ant_units[3,.~pro_trials_ipsi_choice])
    ant_trials_ant2_units_ipsi_choice = vcat(ant_trials_ant_units[2,ant_trials_ipsi_choice],ant_trials_ant_units[3,.~ant_trials_ipsi_choice])
    ant2_units_ipsi_choice = vcat(pro_trials_ant2_units_ipsi_choice,ant_trials_ant2_units_ipsi_choice )

    pro_trials_ant2_units_cont_choice = vcat(pro_trials_ant_units[2,.~pro_trials_ipsi_choice],pro_trials_ant_units[3,pro_trials_ipsi_choice])
    ant_trials_ant2_units_cont_choice = vcat(ant_trials_ant_units[2,.~ant_trials_ipsi_choice],ant_trials_ant_units[3,ant_trials_ipsi_choice])
    ant2_units_cont_choice = vcat(pro_trials_ant2_units_cont_choice,ant_trials_ant2_units_cont_choice )

    # Compute dprime
    pro_units_choice_dprime  = calc_dprime(pro_units_ipsi_choice, pro_units_cont_choice)
    ant1_units_choice_dprime = calc_dprime(ant1_units_ipsi_choice, ant2_units_cont_choice) 
    ant2_units_choice_dprime = calc_dprime(ant2_units_ipsi_choice, ant2_units_cont_choice) 
   
    return pro_units_choice_dprime, ant1_units_choice_dprime, ant2_units_choice_dprime
end
 


