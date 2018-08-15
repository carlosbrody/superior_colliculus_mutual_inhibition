# A set of functions for testing unilateral inactivations
include("pro_anti.jl")
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")

"""
    Iterates through a farm directory and tests all model solutions
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
"""
function test_solution(filename; testruns=1000)
        # load parameters for this solution
        mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")


        # ipsi/contra x pro/anti x control/delay/target/full
        uni_hits = zeros(2,2, 4);

        # override control condition with full trial inactivation
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



















