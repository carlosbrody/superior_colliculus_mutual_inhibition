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

"""
function scatter_by_arg(results, args, arg1,arg2; fignum=1)

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
    
    # plot axis lines
    plot(vec([0 0]), vec([miny maxy]), "k--");   
    plot(vec([minx maxx]), vec([0 0]), "k--");

    # scatter the parameters
    plot(x,y, "bo");

end



