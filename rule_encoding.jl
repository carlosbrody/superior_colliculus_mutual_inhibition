
"""
    build_encoding_dataset(farm_id)
        
builds a database of rule encoding indexes for each farm run in with farmdir/farm_id

# PARAMETERS

-farm_id ID of farm to analyze, eq "C17"

# OPTIONAL PARAMETERS 

-farmdir, directory where farm run info is held

-testruns, number of runs to use to compute encoding indicies

- overrideDict, if you want to override any default parameters. Overrides for all farm runs

# RETURNS (and saves)

- encoding, for N farms X Opto Conditions X pro/anti X hit/miss show the rule encoding strength: mean(pro-anti), positive value ==pro encoding, negative value == anti encoding 

- error_types, N farms X opto Conditions X weak/wrong encoding X pro/anti, what percentage of ERROR trials had correct rule encoding at end of delay. Weak coding is |pro-anti| <= 0.2, wrong is the sign is flipped. 

"""
function build_encoding_dataset(farm_id; farmdir="MiniOptimized",testruns=1000, overrideDict=Dict())
    results     = load(farmdir*farm_id*"_results.jld");
    mypars,extra_pars = load(results["files"][1],"mypars","extra_pars")
    sample_point = Int(floor(mypars[:rule_and_delay_periods][1]/mypars[:dt]))

    # for N farms X Opto Conditions X pro/anti X hit/miss show the rule encoding strength: mean(pro-anti)/var(pro-anti)
    # positive value ==pro encoding, negative value == anti encoding 
    encoding    = Array(Float64,length(results["cost"]),size(extra_pars[:opto_periods],1),2,2);
    
    # for N farms X opto Conditions X weak/wrong encoding X pro/anti, what percentage of ERROR trials had correct rule encoding at end of delay ?
    error_types = Array(Float64,length(results["cost"]),size(extra_pars[:opto_periods],1),2,2);
    for i = 1:length(results["cost"])
        filename = results["files"][i];
        @printf("%d/%d, %s:  \n", i, length(results["cost"]), filename)
        mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")

        opto_string = ["Control", "Delay", "Choice"];
        for j=1:size(extra_pars[:opto_periods],1);
        @printf("%s:  \n", opto_string[j])
        # get set of runs
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
        :opto_times=>reshape(extra_pars[:opto_periods][j,:], 1, 2),
        ))
        proVs, antiVs, proF, antiF = run_ntrials(testruns, testruns; merge(make_dict(args, pars3, these_pars), overrideDict)...); 

        # get index of correct and incorrect trials
        hitBP = find(proVs[1,:]  .> proVs[4,:])
        missBP= find(proVs[1,:]  .< proVs[4,:])
        hitBA = find(antiVs[4,:] .> antiVs[1,:])
        missBA= find(antiVs[4,:] .< antiVs[1,:])              

        # for correct pro and anti trials, get pro-anti for each run
        prosh = proF[1,sample_point,hitBP]' + proF[4,sample_point,hitBP]' - proF[2,sample_point,hitBP]' - proF[3,sample_point,hitBP]';
        prosm = proF[1,sample_point,missBP]'+ proF[4,sample_point,missBP]'- proF[2,sample_point,missBP]'- proF[3,sample_point,missBP]';
        antih = antiF[1,sample_point,hitBA]' + antiF[4,sample_point,hitBA]' - antiF[2,sample_point,hitBA]' - antiF[3,sample_point,hitBA]';
        antim = antiF[1,sample_point,missBA]'+ antiF[4,sample_point,missBA]'- antiF[2,sample_point,missBA]'- antiF[3,sample_point,missBA]';
        # compute mean difference for each pro/anti X hit/miss
        encoding[i,j,1,1] = mean(prosh);#/var(prosh);
        encoding[i,j,1,2] = mean(prosm);#/var(prosm);
        encoding[i,j,2,1] = mean(antih);#/var(antih);
        encoding[i,j,2,2] = mean(antim);#/var(antim);

    # for N farms X opto Conditions X weak encoding/wrong encoding X pro/anti, what percentage of ERROR trials had correct rule encoding at end of delay period?
        error_types[i,j,1,1] = sum(prosm .<= .2)/length(prosm);
        error_types[i,j,1,2] = sum(antim .>= -.2)/length(antim);
        error_types[i,j,2,1] = sum(prosm .<= 0)/length(prosm);
        error_types[i,j,2,2] = sum(antim .>= 0)/length(antim);  
        @printf("Pro  hit/miss:  %.2g / %.2g \n", encoding[i,j,1,1], encoding[i,j,1,2])
        @printf("Anti hit/miss:  %.2g / %.2g \n", encoding[i,j,2,1], encoding[i,j,2,2])

        @printf("Weak Encoding pro/anti:  %.2g / %.2g \n", error_types[i,j,1,1], error_types[i,j,1,2])
        @printf("Wrong Encoding pro/anti:  %.2g / %.2g \n", error_types[i,j,2,1], error_types[i,j,2,2])
        end
        @printf("\n")
    end

    myfilename = farmdir*farm_id*"_encoding.jld";
    save(myfilename, Dict("encoding"=>encoding, "error_types"=>error_types))

    return encoding,error_types
end

# displays the encoding information in the database for the farm run at index <fileindex>
function display_encoding(encoding, error_types, fileindex)
    i = fileindex;
    opto_string = ["Control", "Delay", "Choice"];
    for j=1:3;
        @printf("%s:  \n", opto_string[j])
        @printf("Pro  hit/miss:  %.2g / %.2g \n", encoding[i,j,1,1], encoding[i,j,1,2])
        @printf("Anti hit/miss:  %.2g / %.2g \n", encoding[i,j,2,1], encoding[i,j,2,2])
        @printf("Weak Encoding pro/anti:  %2.g%% / %2.g%% \n", error_types[i,j,1,1]*100, error_types[i,j,1,2]*100)
        @printf("Wrong Encoding pro/anti:  %2.g%% / %2.g%% \n", error_types[i,j,2,1]*100, error_types[i,j,2,2]*100)
    end
    @printf("\n")

end



