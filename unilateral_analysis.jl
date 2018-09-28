# A set of functions for testing unilateral inactivations
include("pro_anti.jl")
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")

"""
    Iterates through a farm directory and tests all model solutions

INPUTS  
    farm_id, Name of farm, ie "C30"

    farmdir, Directory of farm, ie "MiniC30"

OPTIONAL INPUTS
    testruns, number of trials to simulate for each condition of unilateral inactivation

    threshold, the test cost threshold to use to determine which farms to simulate

OUTPUTS
    Saves a file <farmdir>_<farm_id>_unilateral.jld that includes the unilateral hit% for each condition
    
    returns nothing

"""
function test_farm_unilateral(farm_id, farmdir; testruns=1000, threshold=-0.00025, force_opto=false, opto_strength=0.5,numconditions=4)

    # Get list of good farms to test
    results = load_farm_cost_filter(farm_id, farmdir; threshold = threshold)

    uni = zeros(length(results["cost"]),2,2,numconditions);
    for i=1:length(results["cost"])
        filename = results["files"][i];
        @printf("%d/%d, %s:  \n", i, length(results["cost"]), filename)

        # Ipsi/contra x pro/anti x control/delay/target/full
        uni_hits = test_solution(filename; testruns=testruns, force_opto=force_opto, opto_strength=opto_strength, numconditions=numconditions)
        uni[i,:,:,:]= uni_hits;
    end
    uni_results = merge(copy(results),Dict("uni"=>uni));   

    if force_opto 
    myfilename = farmdir*"_"*farm_id*"_unilateral_opto_"*string(opto_strength)*".jld";
    else
    myfilename = farmdir*"_"*farm_id*"_unilateral.jld";
    end
    save(myfilename, Dict("uni_results"=>uni_results))
end



"""
    Tests a set of model parameters with both ipsilateral and contralateral unilateral inactivations

INPUTS
    filename for the model to test

OPTIONAL INPUTS
    testruns, the number of trials to simulate

OUTPUTS
    returns a matrix 2x2x4, which is the hit% for ipsi/contra inactivations x pro/anti trials x control/delay/target/full trial inactivations 

"""
function test_solution(filename; testruns=1000, force_opto=false, opto_strength=0.5,numconditions=4)
        # load parameters for this solution
        mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")

        # ipsi/contra x pro/anti x control/delay/target/full
        uni_hits = zeros(2,2, numconditions);

        # add another opto_periods condition with full trial inactivation
        if numconditions == 4
            extra_pars[:opto_periods] = [extra_pars[:opto_periods]; "trial_start" "trial_end"]
        elseif numconditions == 1
            extra_pars[:opto_periods] = ["trial_start" "trial_end"]   
        end
   
        if force_opto
            pars3[find(args .=="opto_strength")] = opto_strength;
        end
 
        for period=1:size(extra_pars[:opto_periods],1) 
            # set up inputs, including iterating contra/ipsi and opto-condition
            these_pars = merge(mypars, extra_pars);
            these_pars = merge(these_pars, Dict(:opto_units=>[1,2]));
            these_pars = merge(these_pars, Dict(
            :opto_times=>reshape(extra_pars[:opto_periods][period,:], 1, 2),
            :rule_and_delay_period=>these_pars[:rule_and_delay_periods][1], 
            :target_period=>these_pars[:target_periods][1], 
            :post_target_period=>these_pars[:post_target_periods][1], 
            ));
    
            # Ipsilateral 
            proVs, antiVs, proF, antiF = run_ntrials(testruns, testruns; make_dict(args, pars3, these_pars)...); 
            uni_hits[1,1,period] = sum(proVs[1,:]  .>  proVs[4,:])/testruns;        
            uni_hits[1,2,period] = sum(antiVs[4,:] .> antiVs[1,:])/testruns;        
    
            # Contralateral 
            these_pars = merge(these_pars, Dict(:opto_units=>[3,4]));
            proVs, antiVs, proF, antiF = run_ntrials(testruns, testruns; make_dict(args, pars3, these_pars)...); 
            uni_hits[2,1,period] = sum(proVs[1,:]  .>  proVs[4,:])/testruns;        
            uni_hits[2,2,period] = sum(antiVs[4,:] .> antiVs[1,:])/testruns;        
        end

        return uni_hits
end

"""
    Plots the results of unilateral inactivation

INPUTS
    farm_id, name of farm
    farmdir, location of farm files
    
OPTIONAL INPUTS
    color_clusters, if true, plots each cluster separately as a different color
    inact_type, either "full", "delay", or "choice" determines which trial type to plot

"""
function plot_unilateral_psychometric(farm_id, farmdir; color_clusters="PCA", inact_type="full", uni_filters=[])
    # load unilateral hit data
    unilateral = load(farmdir*"_"*farm_id*"_unilateral.jld","uni_results")
    numfarms = size(unilateral["uni"],1)
    uni = unilateral["uni"].*100;
    # farms x ipsi/contra x pro/anti x control/delay/target/full

    # load cluster info
    if color_clusters == "PCA"
        cluster_info = load(farmdir*"_"*farm_id*"_clusters.jld")
        cluster_ids = cluster_info["idx"];
        all_colors = "bgrcmyk";   
        numclusters = sum(.!isnan.(sort(unique(cluster_ids))));
    elseif color_clusters == "uni"
        uni_results, uni_clusters, unidex, ipsi_unidex, contra_unidex, filters = load_farm_unilateral_filter("C32", "MiniC32"; filters=uni_filters);
        all_colors = "gb";   
        numclusters = 2;
        cluster_ids = uni_clusters;
    else
        cluster_ids = ones(1,numfarms);
        numclusters = 1;
    end

    # parse inactivation type
    if inact_type == "full"
        idex = 4;
    elseif inact_type == "delay"
        idex = 2;
    elseif inact_type =="choice"
        idex = 3;
    else
        error("inact_type must be one of 'full','delay', or 'choice'") 
    end

    # iterate over farms    
    figure();
    for i=1:numfarms
        if !isnan(cluster_ids[i])
        co = (cluster_ids[i]-1)/(numclusters*2);
        # for each farm, plot hit data
        plot(1+co, uni[i,2,1,1],"o",color=string(all_colors[Int64(cluster_ids[i])])) # pro control, ipsi/contra irrelevant
        plot(2+co, uni[i,2,2,1],"x",color=string(all_colors[Int64(cluster_ids[i])])) # anti control, ipsi/contra irrelevant

        plot(4+co, uni[i,2,1,idex],"o",color=string(all_colors[Int64(cluster_ids[i])])) # pro full, ipsi
        plot(5+co, uni[i,2,2,idex],"x",color=string(all_colors[Int64(cluster_ids[i])])) # anti full, ipsi

        plot(7+co, uni[i,1,1,idex],"o",color=string(all_colors[Int64(cluster_ids[i])])) # pro full, contra
        plot(8+co, uni[i,1,2,idex],"x",color=string(all_colors[Int64(cluster_ids[i])])) # anti full, contra
        end
    end

    # plot cluster averages
    for i=1:numclusters
        co = (i-1)/(numclusters*2);
        plot(1+co, mean(uni[vec(cluster_ids .== i),2,1,1]),"ko")
        plot(2+co, mean(uni[vec(cluster_ids .== i),2,2,1]),"ko")
        plot(4+co, mean(uni[vec(cluster_ids .== i),2,1,idex]),"ko")
        plot(5+co, mean(uni[vec(cluster_ids .== i),2,2,idex]),"ko")
        plot(7+co, mean(uni[vec(cluster_ids .== i),1,1,idex]),"ko")
        plot(8+co, mean(uni[vec(cluster_ids .== i),1,2,idex]),"ko")
    end

    xticks([1.25,2.25,4.25,5.25,7.25,8.25],["Pro Control", "Anti Control", "Pro Ipsi", "Anti Ipsi", "Pro Contra", "Anti Contra"])
    ylabel("Accuracy %")
    if color_clusters =="uni"
        plot(vec([4 4.5]), vec([uni_filters[1] uni_filters[1]]), "k--")
        plot(vec([5 5.5]), vec([uni_filters[2] uni_filters[2]]), "k--")
        plot(vec([7 7.5]), vec([uni_filters[3] uni_filters[3]]), "k--")
        plot(vec([8 8.5]), vec([uni_filters[4] uni_filters[4]]), "k--")
    end

end


function plot_unilateral_farm(filename; inact_type="full",fignum=1,force_opto=[],testruns=10,alluni=[])
    if (inact_type != "full" ) 
        error("partial inactivations not implemented, use \"full\"")
    end

    # load file parameters
    mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3");
    extra_pars[:opto_periods][2,:] = ["trial_start" "trial_end"];
    extra_pars[:opto_periods][3,:] = ["trial_start" "trial_end"];

    #override opto parameter
    if !isempty(force_opto)
        pars3[find(args .=="opto_strength")] = force_opto[1];
    end

    # set up plotting
    pygui(true)
    figure(fignum); clf();
    pstrings = ["CONTROL", "IPSI", "CONTRA"]

    for period=1:size(extra_pars[:opto_periods],1) 
        # set up inputs, including iterating contra/ipsi and opto-condition
        these_pars = merge(mypars, extra_pars);
        if period == 2
            these_pars = merge(these_pars, Dict(:opto_units=>[1,2]));
        elseif period==3
            these_pars = merge(these_pars, Dict(:opto_units=>[3,4]));
        end
        these_pars = merge(these_pars, Dict(
        :opto_times=>reshape(extra_pars[:opto_periods][period,:], 1, 2),
        :rule_and_delay_period=>these_pars[:rule_and_delay_periods][1], 
        :target_period=>these_pars[:target_periods][1], 
        :post_target_period=>these_pars[:post_target_periods][1], 
        ));

        # plotting set up
        delete!(these_pars, :plot_list)
        pvax = subplot(4,3,period);   axisHeightChange(0.9, lock="t")
        pdax = subplot(4,3,period+3); axisHeightChange(0.9, lock="c"); 
        avax = subplot(4,3,period+6); axisHeightChange(0.9, lock="c")
        adax = subplot(4,3,period+9); axisHeightChange(0.9, lock="b")
        proVs, antiVs = run_ntrials(testruns, testruns; plot_list=[1:20;], plot_Us=false, 
            ax_set = Dict("pro_Vax"=>pvax, "pro_Dax"=>pdax, "anti_Vax"=>avax, "anti_Dax"=>adax),
        make_dict(args, pars3, these_pars)...);

        hBP = length(find(proVs[1,:]  .> proVs[4,:])) /size(proVs, 2)
        hBA = length(find(antiVs[4,:] .> antiVs[1,:]))/size(antiVs,2)
        safe_axes(pvax); title(@sprintf("%s  PRO hits = %.2f%%", pstrings[period], 100*hBP))
        safe_axes(avax); title(@sprintf("ANTI hits = %.2f%%", 100*hBA))
        safe_axes(pdax); remove_xtick_labels(); xlabel("")
        if period > 1
            remove_ytick_labels([pvax, pdax, avax, adax])
        end
        
        figure(fignum)[:canvas][:draw]()
        pause(0.001)

    end
    println("Opto Strength: "*string(pars3[8]))
        
    if !isempty(alluni)
        dex = find(alluni["files"] .== filename);
    println("Control Pro:  "*string(alluni["uni"][dex,1,1,1]))
    println("Control Anti: "*string(alluni["uni"][dex,1,2,1]))

    println("Ipsi Pro:     "*string(alluni["uni"][dex,1,1,4]))
    println("Ipsi Anti:    "*string(alluni["uni"][dex,1,2,4]))

    println("Contra Pro:   "*string(alluni["uni"][dex,2,1,4]))
    println("Contra Anti:  "*string(alluni["uni"][dex,2,2,4]))
    end
end 

function unilateral_by_strength(filename;testruns=100, plot_stuff=true)

        mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")
        og_opto=pars3[8];
        # ipsi/contra x pro/anti x control/delay/target/full
        uni_hits = zeros(11,2,2);

        # add another opto_periods condition with full trial inactivation
        extra_pars[:opto_periods] = [extra_pars[:opto_periods]; "trial_start" "trial_end"]
        extra_pars[:opto_periods][1,:] = ["trial_start" "trial_end"];  
        extra_pars[:opto_periods][2,:] = ["trial_start" "trial_end"];
        extra_pars[:opto_periods][3,:] = ["trial_start" "trial_end"];  
 
        optos = 0:0.1:1;
        for strength=1:11
            # set up inputs, including iterating contra/ipsi and opto-condition
            these_pars = merge(mypars, extra_pars);
            these_pars = merge(these_pars, Dict(:opto_units=>[1,2]));
            these_pars = merge(these_pars, Dict(
            :opto_times=>reshape(extra_pars[:opto_periods][1,:], 1, 2),
            :rule_and_delay_period=>these_pars[:rule_and_delay_periods][1], 
            :target_period=>these_pars[:target_periods][1], 
            :post_target_period=>these_pars[:post_target_periods][1], 
            ));
            
            # force opto
            pars3[find(args .=="opto_strength")] = optos[strength];
 
            # Ipsilateral 
            proVs, antiVs = run_ntrials(testruns, testruns; make_dict(args, pars3, these_pars)...); 
            uni_hits[strength,1,1] = sum(proVs[1,:]  .>  proVs[4,:])/testruns;        
            uni_hits[strength,1,2] = sum(antiVs[4,:] .> antiVs[1,:])/testruns;        
    
            # Contralateral 
            these_pars = merge(these_pars, Dict(:opto_units=>[3,4]));
            proVs, antiVs = run_ntrials(testruns, testruns; make_dict(args, pars3, these_pars)...); 
            uni_hits[strength,2,1] = sum(proVs[1,:]  .>  proVs[4,:])/testruns;        
            uni_hits[strength,2,2] = sum(antiVs[4,:] .> antiVs[1,:])/testruns;        
        end

        if plot_stuff
        fig, ax=subplots()
        plot(optos, 90.*ones(11,1),"k-",alpha=0.25, label="Pro target")
        plot(optos, 70.*ones(11,1),"k--",alpha=0.25,label="Anti target")
        plot(optos, uni_hits[:,1,1].*100,"b-",label="Pro Go Ipsi")
        plot(optos, uni_hits[:,2,2].*100,"g--",label="Anti Go Ipsi")
        plot(optos, uni_hits[:,2,1].*100,"m-",label="Pro Go Contra")
        plot(optos, uni_hits[:,1,2].*100,"r--",label="Anti Go Contra")
        plot(vec([og_opto og_opto]), vec([0 100]), "k--",alpha=0.25)
        ax[:legend](loc="lower right")
        ylabel("Accuracy %")
        xlabel("Opto Strength (0=full)")
        title(filename)
        end        

        return uni_hits
end

"""
    Plots the results of unilateral inactivation

INPUTS
    farm_id, name of farm
    farmdir, location of farm files
    
OPTIONAL INPUTS
    color_clusters, if true, plots each cluster separately as a different color
    inact_type, either "full", "delay", or "choice" determines which trial type to plot

"""
function plot_forced_unilateral_psychometric(farm_id, farmdir, force_strength; color_clusters=true)
    # load unilateral hit data
    unilateral = load(farmdir*"_"*farm_id*"_unilateral_opto_"*string(force_strength)*".jld","uni_results");
    numfarms = size(unilateral["uni"],1)
    uni = unilateral["uni"].*100;
    # farms x ipsi/contra x pro/anti x control/delay/target/full

    # load cluster info
    if color_clusters
        cluster_info = load(farmdir*"_"*farm_id*"_clusters.jld")
        cluster_ids = cluster_info["idx"];
        all_colors = "bgrcmyk";   
        numclusters = sum(.!isnan.(sort(unique(cluster_ids))));
    else
        cluster_ids = ones(1,numfarms);
        numclusters = 1;
    end

    # iterate over farms    
    figure();
    for i=1:numfarms
        if !isnan(cluster_ids[i])
        co = (cluster_ids[i]-1)/(numclusters*2);
        # for each farm, plot hit data
        plot(1+co, uni[i,1,1,1],"o",color=string(all_colors[Int64(cluster_ids[i])])) # pro full, ipsi
        plot(2+co, uni[i,1,2,1],"x",color=string(all_colors[Int64(cluster_ids[i])])) # anti full, ipsi

        plot(4+co, uni[i,2,1,1],"o",color=string(all_colors[Int64(cluster_ids[i])])) # pro full, contra
        plot(5+co, uni[i,2,2,1],"x",color=string(all_colors[Int64(cluster_ids[i])])) # anti full, contra
        end
    end

    # plot cluster averages
    for i=1:numclusters
        co = (i-1)/(numclusters*2);
        plot(1+co, mean(uni[vec(cluster_ids .== i),1,1,1]),"ko")
        plot(2+co, mean(uni[vec(cluster_ids .== i),1,2,1]),"ko")
        plot(4+co, mean(uni[vec(cluster_ids .== i),2,1,1]),"ko")
        plot(5+co, mean(uni[vec(cluster_ids .== i),2,2,1]),"ko")
    end

    xticks([1.25,2.25,4.25,5.25],["Pro Ipsi", "Anti Ipsi", "Pro Contra", "Anti Contra"])
    ylabel("Accuracy %")

end


function load_farm_unilateral_filter(farmid, farmdir; filters = [90 80 60 40], plot_check=true)

# Get list of all solutions below threshold
# Ipsi/contra x pro/anti x control/delay/target/full
# Ipsi/Contra in MODEL SPACE not RAT SPACE
unilateral  = load(farmdir*"_"*farmid*"_unilateral.jld","uni_results")
numfarms    = size(unilateral["uni"],1)
uni         = unilateral["uni"].*100;

# filter by anti and pro cost
# filters = [ABOVE Pro/Ipsi, ABOVE Anti/Ipsi, BELOW Pro/Contra, BELOW Anti/Contra]
ipsi_unidex     = (uni[:,2,1,4] .>= filters[1]) .& (uni[:,2,2,4] .>= filters[2]);
contra_unidex   = (uni[:,1,1,4] .<= filters[3]) .& (uni[:,1,2,4] .<= filters[4]);
unidex          = ipsi_unidex .& contra_unidex;

# PLOT for checking
if plot_check
    figure()
    subplot(121)
    plot(uni[:,1,1,1], uni[:,1,2,1],"rx")
    plot(uni[:,2,1,4], uni[:,2,2,4],"bo")
    plot(uni[ipsi_unidex,2,1,4], uni[ipsi_unidex,2,2,4],"gx")
    xlabel("pro hit%")
    ylabel("anti hit%")
    title("go ipsi")
    plot(vec([filters[1] filters[1]]),vec([filters[2] 100]),"k--")
    plot(vec([filters[1] 100]),vec([filters[2] filters[2]]),"k--")
    
    # GO contra trials
    subplot(122)
    plot(uni[:,1,1,1], uni[:,1,2,1],"rx")
    plot(uni[:,1,1,4], uni[:,1,2,4],"bo")
    plot(uni[contra_unidex,1,1,4], uni[contra_unidex,1,2,4],"gx")
    xlabel("pro hit%")
    ylabel("anti hit%")
    title("go contra")
    plot(vec([filters[3] filters[3]]),vec([0 filters[4]]),"k--")
    plot(vec([0 filters[3]]),vec([filters[4] filters[4]]),"k--")
end

# return results matrix filtered by unidex
uni_results = Dict()
uni_results["cost"]     = unilateral["cost"][vec(unidex)]
uni_results["tcost"]    = unilateral["tcost"][vec(unidex)]
uni_results["files"]    = unilateral["files"][vec(unidex)]
uni_results["dirs"]     = unilateral["dirs"][vec(unidex)]
uni_results["params"]   = unilateral["params"][vec(unidex),:]
uni_results["uni"]      = unilateral["uni"][vec(unidex),:,:,:]

# make a cluster index vector for plotting purposes
uni_clusters = zeros(numfarms,1)
uni_clusters[unidex] = 1
uni_clusters[.!unidex] = 2;

return uni_results, uni_clusters, unidex, ipsi_unidex, contra_unidex, filters

end
