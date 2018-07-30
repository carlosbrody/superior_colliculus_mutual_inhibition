# this file contains functions for plotting the parameters across an entire farm
#   Alex Piet (07/25/2018)
#
#   List top-level of functions
#   scatter_by_arg()
#   histo_params_by_cluster()    
#   plot_psychometric()
#
#
#
#


"""
    scatter_by_arg(results, args, arg1, arg2)

    Plots a scatter plot of all the parameters in results for parameters arg1 and arg2.

    PARAMETERS
    results     a Dict() that contains entries for "params" a matrix of all parameters for all runs
    args        a list of parameter labels (strings)
    arg1        a string for the first parameter to plot
    arg2        a string for the second parameter to plot

    SETUP and EXAMPLE(if you need it)
    include("svd_cluster.jl")
    f1      = load("MiniC30/farm_C30_Farms_C30_spock_brody01-01_0001.jld")
    args    = f1["args"];
    results = load_farm_cost_filter("C30", "MiniC30")
    scatter_by_arg(results, args, "hW_A","hW_P");
    
    #to include cluster labels
    cluster_info    = load("MiniC30_C30_clusters.jld")
    cluster_ids     = cluster_info["idx"]
    scatter_by_arg(results, args, "hW_A", "hW_P"; cluster_ids = cluster_ids)

"""
function scatter_by_arg(results, args, arg1,arg2; fignum=1, cluster_ids=[])

    # check for clustering
    clustering = !isempty(cluster_ids);

    if !clustering
        scatter_by_arg_backend(results, args, arg1, arg2; fignum=fignum);
    else
        all_colors = "bgrcmyk";
        # iterate over clusters
        ids = sort(unique(cluster_ids));
        for i=1:length(ids)
            # check for NaN
            if !isnan(ids[i])
                # filter results, and call backend
                fresults    = filter_results_by_cluster(results,cluster_ids, ids[i]);
                color       = string(all_colors[i]);
                plot_lines  = i==1;
                scatter_by_arg_backend(fresults, args, arg1, arg2; fignum=fignum, color=color, plot_lines=plot_lines);
            end
        end
    end
end


# helper function
function scatter_by_arg_backend(results, args, arg1, arg2; fignum=1, color="k", plot_lines=true);

    # get index of each parameter from args
    index1 = find(args .== arg1);
    index2 = find(args .== arg2);

    # do some control checks for mis-matched argument labels
    if isempty(index1)
        error("Parameter 1 not found in argument list")
    end
    if isempty(index2)
        error("Parameter 2 not found in argument list")
    end
    if length(index1) > 1
        error("More than 1 match for parameter 1 in argument list")
    end
    if length(index2) > 1
        error("More than 1 match for parameter 2 in argument list")
    end
    
    # open figure
    figure(fignum)

    # get the parameters we want
    x = results["params"][:,index1];
    y = results["params"][:,index2];
    
    # find extrema of parameters for plotting dashed axis lines of appropriate size
    maxx = maximum(x);
    maxy = maximum(y);
    minx = minimum(x);
    miny = minimum(y);
    maxx = 3; maxy = 3;
    minx = -3; miny = -3;    

    # plot axis lines
    if plot_lines
        plot(vec([0 0]), vec([miny maxy]), "k--");   
        plot(vec([minx maxx]), vec([0 0]), "k--");
    end

    # scatter the parameters
    plot(x,y, "o",color=color);
    
    ylabel(arg2)
    xlabel(arg1)

end


"""
    helper function for filtering a dictionary "results" by cluster id

"""
function filter_results_by_cluster(results, cluster_ids, target_cluster)
    
    # make index of good clusters
    gooddex = cluster_ids .== target_cluster

    # filter
    r1 = Dict();
    r1["cost"]  = results["cost"][vec(gooddex)];
    r1["tcost"] = results["tcost"][vec(gooddex)];
    r1["files"] = results["files"][vec(gooddex)];
    r1["dirs"]  = results["dirs"][vec(gooddex)];
    r1["params"]= results["params"][vec(gooddex),:]

    # return filtered matrix
    return r1

end

"""
    Wrapper function that uses the "histo_param()" function to plot histograms for each cluster

"""
function histo_params_by_cluster(results, args, cluster_ids, target_cluster; threshold=-0.0002, cost_choice="cost", further_params...)

    r1 = filter_results_by_cluster(results,cluster_ids, target_cluster);
    histo_params(args, r1["params"], r1["tcost"], r1["cost"], r1["files"]; Dict(further_params)...);


end

"""
    function for plotting the psychometric behavior of all the runs on a farm, can handle coloring each cluster

"""
function plot_psychometric(results; hit_type="standard", color_clusters=false, cluster_ids=[],minidir=true, plot_only = 10)
    println("this function is slow because I have to load many files...")
    println("set plot_only = 10 to plot a smaller number of files")


    if color_clusters && size(results["files"],1) != size(cluster_ids,1)
        error("cluster_ids doesn't match size of results")
    end

    if hit_type == "standard"
        pro_label = "hP"
        anti_label ="hA"
    else
        pro_label = "hBP"
        anti_label ="hBA"
    end

    figure()
    if color_clusters
        all_colors = "bgrcmyk";
        uclusters = sort(unique(cluster_ids));
        numclusters = sum(.!isnan.(uclusters));
    end

    # get list of farms from results
    num2plot = Int64(minimum([length(results["files"]), plot_only]));
    for i=1:num2plot
        if rem(i, 50) == 0
            println(i)
        end
        if minidir
            temp = load(results["files"][i],"files");
            P,A = load(temp,pro_label, anti_label);
        else
            P,A = load(results["files"][i],pro_label, anti_label);
        end
        if i == 1
            extra_pars = load(results["files"][1],"extra_pars");
            targets = extra_pars[:opto_targets].*100;
            plot(vec([1 1.5]), vec([targets[1,1] targets[1,1]]),"k")
            plot(vec([2 2.5]), vec([targets[2,1] targets[2,1]]),"k")
            plot(vec([3 3.5]), vec([targets[3,1] targets[3,1]]),"k")
            plot(vec([1 1.5]), vec([targets[1,2] targets[1,2]]),"r")
            plot(vec([2 2.5]), vec([targets[2,2] targets[2,2]]),"r")
            plot(vec([3 3.5]), vec([targets[3,2] targets[3,2]]),"r")
        end
        if color_clusters && !isnan(cluster_ids[i])
            co = (cluster_ids[i]-1)/(numclusters*2)
            plot(vec([1 2 3]+co), P.*100, "o", color=string(all_colors[Int64(cluster_ids[i])]))
            plot(vec([1 2 3]+co), A.*100, "x", color=string(all_colors[Int64(cluster_ids[i])]))
        else
            plot(vec([1 2 3]), P.*100, "ko")
            plot(vec([1 2 3]), A.*100, "rx")
        end

    end
    ylim([0, 100])
    xlim([.8, 3.7]);
    ylabel("Accuracy %")
    xlabel("Opto Condition")
    title(results["dirs"][1]*"--"*hit_type)
    xticks([1.25, 2.25, 3.25], ["Control", "Delay", "Choice"])
end






