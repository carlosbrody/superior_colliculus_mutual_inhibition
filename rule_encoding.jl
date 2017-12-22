README="""
    LIST OF ALL FUNCTIONS IN THIS FILE

All functions should work regardless of how many and which parameters were optimized in the model

# BACKEND FUNCTIONS
    build_encoding_dataset(farm_id)

# DATA EXPLORATION FUNCTIONS
    display_encoding(encoding, error_types, fileindex)
    scatter_svd_by_index(data, encoding, error_types, opto_type, encode_index,svd_index)

"""



"""
    build_encoding_dataset(farm_id)
        
builds a database of rule encoding indexes for each farm run in with farmdir/farm_id

# PARAMETERS

-farm_id        ID of farm to analyze, eq "C17"

# OPTIONAL PARAMETERS 

-farmdir        directory where farm run info is held

-testruns       number of runs to use to compute encoding indicies

- overrideDict  if you want to override any default parameters. Overrides for all farm runs

# RETURNS (and saves)

- encoding      A matrix: N farms X Opto Conditions X pro/anti X hit/miss X 3 encoding indicies. (1) Pro-anti, (2) R-L, (3) diagonal1-diagonal2

- error_types   A matrix: N farms X Opto Conditions X pro/anti X hit/miss X 3 encoding indicies. (1) Pro-anti, (2) R-L, (3) diagonal1-diagonal2. Instead of showing the average encoding value, it shows the fraction of trials with a positive index. 


"""
function build_encoding_dataset(farm_id; farmdir="MiniOptimized",testruns=1000, overrideDict=Dict(), update_only=false, old_results=Dict())
    results     = load(farmdir*"_"*farm_id*"_results.jld");
    mypars,extra_pars = load(results["files"][1],"mypars","extra_pars")
    sample_point = Int(floor(mypars[:rule_and_delay_periods][1]/mypars[:dt]))

    # for N farms X Opto Conditions X pro/anti X hit/miss show theencoding strength: Three indexes
    # 1. pro-anti       positive value ==pro encoding, negative value == anti encoding 
    # 2. right-left     positive value ==right choice, negative value == left choice
    # 3. diagonal       (proR+antiL) - (proL+antiR)
    encoding    = Array(Float64,length(results["cost"]),size(extra_pars[:opto_periods],1),2,2,3);
    
    # for N farms X opto Conditions X pro/anti x hit/miss x three indexes, what percentage of trials had each encoding at end of delay?
    error_types = Array(Float64,length(results["cost"]),size(extra_pars[:opto_periods],1),2,2,3);
    
    if update_only
        myfilename = farmdir*"_"*farm_id*"_encoding.jld";
        old_encoding, old_error_types = load(myfilename, "encoding","error_types")
    end

    for i = 1:length(results["cost"])
        filename = results["files"][i];
        @printf("%d/%d, %s:  \n", i, length(results["cost"]), filename)
        mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")

        if update_only & (size(find(filename .== old_results["files"]),1) > 0)
            old_index = find(filename .== old_results["files"]);
            encoding[i,:,:,:,:] = old_encoding[old_index, :,:,:,:];
            error_types[i,:,:,:,:] = old_error_types[old_index, :,:,:,:];
        else
        for j=1:size(extra_pars[:opto_periods],1);
        # get set of runs
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
        :opto_times=>reshape(extra_pars[:opto_periods][j,:], 1, 2),
        :rule_and_delay_period=>these_pars[:rule_and_delay_periods][1], 
        :target_period=>these_pars[:target_periods][1], 
        :post_target_period=>these_pars[:post_target_periods][1], 
        ))

        proVs, antiVs, proF, antiF = run_ntrials(testruns, testruns; merge(make_dict(args, pars3, these_pars), overrideDict)...); 

        # get index of correct and incorrect trials
        hitBP = find(proVs[1,:]  .> proVs[4,:])
        missBP= find(proVs[1,:]  .< proVs[4,:])
        hitBA = find(antiVs[4,:] .> antiVs[1,:])
        missBA= find(antiVs[4,:] .< antiVs[1,:])              

        # INDEX 1
        # for correct pro and anti trials, get pro-anti for each trial
        prosh = proF[1,sample_point,hitBP]' + proF[4,sample_point,hitBP]' - proF[2,sample_point,hitBP]' - proF[3,sample_point,hitBP]';
        prosm = proF[1,sample_point,missBP]'+ proF[4,sample_point,missBP]'- proF[2,sample_point,missBP]'- proF[3,sample_point,missBP]';
        antih = antiF[1,sample_point,hitBA]' + antiF[4,sample_point,hitBA]' - antiF[2,sample_point,hitBA]' - antiF[3,sample_point,hitBA]';
        antim = antiF[1,sample_point,missBA]'+ antiF[4,sample_point,missBA]'- antiF[2,sample_point,missBA]'- antiF[3,sample_point,missBA]';
        # compute mean difference for each pro/anti X hit/miss
        encoding[i,j,1,1,1] = mean(prosh);
        encoding[i,j,1,2,1] = mean(prosm);
        encoding[i,j,2,1,1] = mean(antih);
        encoding[i,j,2,2,1] = mean(antim);
        # for N farms X opto Conditions X pro/anti x hit/miss, percentage of trials with sign of each index?
        error_types[i,j,1,1,1] = sum(prosh .>= 0)/length(prosh);
        error_types[i,j,1,2,1] = sum(prosm .>= 0)/length(prosm);
        error_types[i,j,2,1,1] = sum(antih .>= 0)/length(antih);
        error_types[i,j,2,2,1] = sum(antim .>= 0)/length(antim);  

        # INDEX 2
        # for correct pro and anti trials, get right-left for each run
        prosh = proF[1,sample_point,hitBP]' + proF[2,sample_point,hitBP]' - proF[3,sample_point,hitBP]' - proF[4,sample_point,hitBP]';
        prosm = proF[1,sample_point,missBP]'+ proF[2,sample_point,missBP]'- proF[3,sample_point,missBP]'- proF[4,sample_point,missBP]';
        antih = antiF[1,sample_point,hitBA]' + antiF[2,sample_point,hitBA]' - antiF[3,sample_point,hitBA]' - antiF[4,sample_point,hitBA]';
        antim = antiF[1,sample_point,missBA]'+ antiF[2,sample_point,missBA]'- antiF[3,sample_point,missBA]'- antiF[4,sample_point,missBA]';
        # compute mean difference for each pro/anti X hit/miss
        encoding[i,j,1,1,2] = mean(prosh);
        encoding[i,j,1,2,2] = mean(prosm);
        encoding[i,j,2,1,2] = mean(antih);
        encoding[i,j,2,2,2] = mean(antim);
        # for N farms X opto Conditions X pro/anti x hit/miss, percentage of trials with sign of each index?
        error_types[i,j,1,1,2] = sum(prosh .>= 0)/length(prosh);
        error_types[i,j,1,2,2] = sum(prosm .>= 0)/length(prosm);
        error_types[i,j,2,1,2] = sum(antih .>= 0)/length(antih);
        error_types[i,j,2,2,2] = sum(antim .>= 0)/length(antim);  

        # INDEX 3
        # for correct pro and anti trials, get diagonal1-diagonal2 for each run
        prosh = proF[1,sample_point,hitBP]' + proF[3,sample_point,hitBP]' - proF[2,sample_point,hitBP]' - proF[4,sample_point,hitBP]';
        prosm = proF[1,sample_point,missBP]'+ proF[3,sample_point,missBP]'- proF[2,sample_point,missBP]'- proF[4,sample_point,missBP]';
        antih = antiF[1,sample_point,hitBA]' + antiF[3,sample_point,hitBA]' - antiF[2,sample_point,hitBA]' - antiF[4,sample_point,hitBA]';
        antim = antiF[1,sample_point,missBA]'+ antiF[3,sample_point,missBA]'- antiF[2,sample_point,missBA]'- antiF[4,sample_point,missBA]';
        # compute mean difference for each pro/anti X hit/miss
        encoding[i,j,1,1,3] = mean(prosh);
        encoding[i,j,1,2,3] = mean(prosm);
        encoding[i,j,2,1,3] = mean(antih);
        encoding[i,j,2,2,3] = mean(antim);
        # for N farms X opto Conditions X pro/anti x hit/miss, percentage of trials with sign of each index?
        error_types[i,j,1,1,3] = sum(prosh .>= 0)/length(prosh);
        error_types[i,j,1,2,3] = sum(prosm .>= 0)/length(prosm);
        error_types[i,j,2,1,3] = sum(antih .>= 0)/length(antih);
        error_types[i,j,2,2,3] = sum(antim .>= 0)/length(antim);  

        end
        end
        display_encoding(encoding,error_types,i)
        @printf("\n")
    end

    myfilename = farmdir*"_"*farm_id*"_encoding.jld";
    save(myfilename, Dict("encoding"=>encoding, "error_types"=>error_types))

    return encoding,error_types
end

"""
    display_encoding(encoding, error_types, fileindex)

displays the encoding information in the database for the farm run at index <fileindex>

# PARAMETERS

- encoding      result from build_encoding_dataset()

- error_types   result from build_encoding_dataset()

- fileindex     the index of which farm to display

"""
function display_encoding(encoding, error_types, fileindex)
    i = fileindex;
    opto_string = ["Control", "Delay", "Choice"];
    for j=1:3;
        @printf("%s:  \n", opto_string[j])
        @printf("Pro  HIT  :%10.2g %10.2g %10.2g \n", encoding[i,j,1,1,1], encoding[i,j,1,1,2], encoding[i,j,1,1,3])
        @printf("Pro  MISS :%10.2g %10.2g %10.2g \n", encoding[i,j,1,2,1], encoding[i,j,1,2,2], encoding[i,j,1,2,3])
        @printf("Anti HIT  :%10.2g %10.2g %10.2g \n", encoding[i,j,2,1,1], encoding[i,j,2,1,2], encoding[i,j,2,1,3])
        @printf("Anti MISS :%10.2g %10.2g %10.2g \n", encoding[i,j,2,2,1], encoding[i,j,2,2,2], encoding[i,j,2,2,3])
        @printf("-----------\n")
        @printf("Pro  HIT  :%9d%% %9d%% %9d%% \n",round(error_types[i,j,1,1,1]*100), round(error_types[i,j,1,1,2]*100),round(error_types[i,j,1,1,3]*100))
        @printf("Pro  MISS :%9d%% %9d%% %9d%% \n",round(error_types[i,j,1,2,1]*100), round(error_types[i,j,1,2,2]*100),round(error_types[i,j,1,2,3]*100))
        @printf("Anti HIT  :%9d%% %9d%% %9d%% \n",round(error_types[i,j,2,1,1]*100), round(error_types[i,j,2,1,2]*100),round(error_types[i,j,2,1,3]*100))
        @printf("Anti MISS :%9d%% %9d%% %9d%% \n",round(error_types[i,j,2,2,1]*100), round(error_types[i,j,2,2,2]*100),round(error_types[i,j,2,2,3]*100))
        @printf("\n")
    end
    @printf("\n")

end


"""
    Plots SVD coordinates aginst encoding indicies

# PARAMETERS

- data, output from SVD_interactive(), the coordinates in SVD dimensions

- encoding, output from SVD_interactive(), the encoding indexes

- error_types, output from SVD_interactive(), the encoding fractions

- opto_type, = 1,2,3 for Control Delay, or Choice inactivation

- encode_index, = 1,2,3 for pro/anti index, right/left index, or diagonal1/diagonal2 index

- svd_index, =1 to ndims, where ndims is the number of SVD dimensions requested from SVD_interactive()

# EXAMPLE

data, filenames, encoding, error_types = SVD_interactive("C17";farmdir="MiniOptimized", backend_mode=true);

Plot the Control data, Pro/Anti Index, SVD Dim 1

scatter_svd_by_index(data,encoding, error_types, 1, 1, 1)

"""
function scatter_svd_by_index(data, encoding, error_types, opto_type, encode_index,svd_index)
     index_labels=["Pro/Anti Index", "Right/Left Index", "Diag1/Diag2 Index"];
     svd_label = "SVD Dim ";
     figure()
     subplot(2,2,1)
     scatter(data[:,svd_index], encoding[:,opto_type,1,1,encode_index])
     title("Pro-hit")
     ylabel(index_labels[encode_index])

     
     subplot(2,2,2)
     scatter(data[:,svd_index], encoding[:,opto_type,1,2,encode_index])
     title("Pro-Miss")
     
     subplot(2,2,3)
     scatter(data[:,svd_index], encoding[:,opto_type,2,1,encode_index])
     title("Anti-hit")
     ylabel(index_labels[encode_index])
     xlabel(svd_label*string(svd_index))
     
     subplot(2,2,4)
     scatter(data[:,svd_index], encoding[:,opto_type,2,2,encode_index])
     title("Anti-Miss")
     xlabel(svd_label*string(svd_index))
end

