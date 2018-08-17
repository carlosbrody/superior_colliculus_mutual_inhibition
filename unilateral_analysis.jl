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
function test_farm_unilateral(farm_id, farmdir; testruns=1000, threshold=-0.00025)

    # Get list of good farms to test
    results = load_farm_cost_filter(farm_id, farmdir; threshold = threshold)

    uni = zeros(length(results["cost"]),2,2,4);
    for i=1:length(results["cost"])
        filename = results["files"][i];
        @printf("%d/%d, %s:  \n", i, length(results["cost"]), filename)

        # Ipsi/contra x pro/anti x control/delay/target/full
        uni_hits = test_solution(filename; testruns=testruns)
        uni[i,:,:,:]= uni_hits;
    end
    uni_results = merge(copy(results),Dict("uni"=>uni));   
 
    myfilename = farmdir*"_"*farm_id*"_unilateral.jld";
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
function test_solution(filename; testruns=1000)
        # load parameters for this solution
        mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")

        # ipsi/contra x pro/anti x control/delay/target/full
        uni_hits = zeros(2,2, 4);

        # add another opto_periods condition with full trial inactivation
        extra_pars[:opto_periods] = [extra_pars[:opto_periods]; "trial_start" "trial_end"]
    
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
function plot_unilateral_psychometric(farm_id, farmdir; color_clusters=true, inact_type="full")
    # load unilateral hit data
    unilateral = load(farmdir*"_"*farm_id*"_unilateral.jld","uni_results")
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

        plot(4+co, uni[i,1,1,idex],"o",color=string(all_colors[Int64(cluster_ids[i])])) # pro full, ipsi
        plot(5+co, uni[i,1,2,idex],"x",color=string(all_colors[Int64(cluster_ids[i])])) # anti full, ipsi

        plot(7+co, uni[i,2,1,idex],"o",color=string(all_colors[Int64(cluster_ids[i])])) # pro full, contra
        plot(8+co, uni[i,2,2,idex],"x",color=string(all_colors[Int64(cluster_ids[i])])) # anti full, contra
        end
    end

    # plot cluster averages
    for i=1:numclusters
        co = (i-1)/(numclusters*2);
        plot(1+co, mean(uni[vec(cluster_ids .== i),2,1,1]),"ko")
        plot(2+co, mean(uni[vec(cluster_ids .== i),2,2,1]),"ko")
        plot(4+co, mean(uni[vec(cluster_ids .== i),1,1,idex]),"ko")
        plot(5+co, mean(uni[vec(cluster_ids .== i),1,2,idex]),"ko")
        plot(7+co, mean(uni[vec(cluster_ids .== i),2,1,idex]),"ko")
        plot(8+co, mean(uni[vec(cluster_ids .== i),2,2,idex]),"ko")
    end

    xticks([1.25,2.25,4.25,5.25,7.25,8.25],["Pro Control", "Anti Control", "Pro Ipsi", "Anti Ipsi", "Pro Contra", "Anti Contra"])
    ylabel("Accuracy %")

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







