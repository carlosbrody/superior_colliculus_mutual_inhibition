# this file contains functions for plotting the parameters across an entire farm
#   Alex Piet (07/25/2018)
#
#   List of functions
#   scatter_by_arg()
#
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
