# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.



if !isdefined(:plot_PA)
    @printf("Loading pro_anti.jl\n")
    include("pro_anti.jl")
end

using HDF5

####################################################################
#
#   Define some helper functions for loading farms,
#   plotting histograms over the parameters, and managing
#   with GUI interactivity.
#
####################################################################


"""
results = farmload(farm_id; farmdir="../NewFarms", verbose=true, verbose_every=10)

Package and return a summary of results from a number of runs.

# PARAMETERS:

- filename     farm_id  (e.g., "C17"); .jld files containing this pattern will be loaded

# OPTIONAL PARAMETERS:

- farmdir      String indicating directory in which the files are found; alternatively, a vector
               of strings; runs from all of those directories will be loaded.

- verbose      If true, print a brief progress report every verbose_every files loaded

- verbose_every  Used if verbose==true (see above)

# RETURNS:

A dictionary with the keys:
    - "tcost" (vector or training costs),
    - "cost" (vector of test costs),
    - "dirs" (vector directories in which files  are found),
    - "files" (vector of filenames, including paths),
    - "qu_out" (deprecated for farm C17 and other non-pre-search farms),
    - "params"" (nruns-by-nparams matrix of paramters),
    - "args" (nparams-long vector of parameter names).
    - "grads" (nruns-by-nparams matrix of training cost gradients). If not available in the file, all elements will be "nothing".
    - "hessians" (nruns vector of nparams-by-nparams hessian matrices). If not available in the file, all elements will be "nothing".

Note that most values are nruns-long vectors but "args" is only nparams-long, because they are
the same names across all the runs.
"""
function farmload(farm_id; farmdir="../NewFarms", verbose=true, verbose_every=50)

    if typeof(farmdir)==String; farmdir=[farmdir]; end

    results = Dict(); dirs=[]; files =[]; qs=[]; tcosts=[]; costs=[]; pars=[]; hBPs=[]; hBAs=[];
    grads = []; hessians = Array{Array{Float64}}(0,1);
    n=0;
    for dd in farmdir
        for f in filter(x -> startswith(x, "farm_" * farm_id * "_"), readdir(dd * "/"))
            n += 1
            myfile = dd * "/" * f;
            fj = jldopen(myfile);
            if exists(fj, "qu_out"); qu_out = load(myfile, "qu_out");
            else;                    qu_out = [-1 -1];
            end
            # if exists(fj, "ftraj3"); ftraj3 = load(myfile, "ftraj3"); mygrad=ftraj3[1,end]'; myhess=ftraj3[2,end]
            # else
                ftraj3 = nothing;            mygrad=nothing;        myhess=nothing
            # end
            close(fj)

            args, params, traj3, cost, hBP, hBA = load(myfile, "args", "pars3", "traj3", "cost", "hBP", "hBA")

            if     !haskey(results, args);       results["args"] = args;
            elseif !all(results["args"].==args); error("Not all files have same args!");
            end

            files = [files ; myfile]; dirs = [dirs ; dd]
            qs = [qs ; qu_out]; tcosts = [tcosts; traj3[2,end]]; costs = [costs; cost]
            if length(pars) ==0; pars =params';  else pars  = [pars  ; params']; end
            if length(hBPs) ==0; hBPs =hBP';     else hBPs  = [hBPs  ; hBP'];    end
            if length(hBAs) ==0; hBAs =hBA';     else hBAs  = [hBAs  ; hBA'];    end

            if mygrad==nothing
                if length(grads)==0; grads=[mygrad];   else grads = [grads ; [mygrad]]; end
            else
                if length(grads)==0; grads=mygrad;     else grads = [grads ; mygrad];   end
            end
            if myhess==nothing
                if length(hessians)==0; hessians = [nothing];
                else hessians = [hessians ; [nothing]];
                end
            else
                hessians = [hessians ; [myhess]]
            end

            if verbose && rem(n, verbose_every)==0
                @printf("%s %g\n", myfile, tcosts[end])
            end
        end
    end

    results["dirs"]     = dirs
    results["files"]    = files
    results["qu_out"]   = qs
    results["tcost"]    = tcosts
    results["cost"]     = costs
    results["params"]   = pars
    results["grads"]    = grads
    results["hessians"] = hessians
    results["hBP"]      = hBPs
    results["hBA"]      = hBAs

    return results
end


# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.



"""
Data structure for interactive scatterplots

`interactive_scatter()` returns one of these objects; used to manage GUI handling.

A mutlidimensional set of data points is plotted in a set of scatterplots, each of which
shows the scatterplot of one dimension against another. Each data point is identified
by a unique string.  When the user clicks on one of the axes, the closest point to it
is identified, and if defined, a callback function is called, with the string ID that
datapoint passed as one of its parameters.

The points can be divided into different subsets of them, with each set plotted
in its own color. In addition, a series of initially invisible dots can also be
added to the plot. If scatter_highlight() is defined as the callback function, then
these initially invidible points will become sequentially visible as the user clicks
on the scatterplots.

The main function that is called to set things up is `interactive_scatter()`. The
usual callback is `scatter_highlight()`.

These functions work with a scatter_data data structure, defined here.

The data structure's fields are:

    Data::Array{Float64}     # npoints-by-ndims Array{Float64} of original data
    stringIDs::Array{String} # npoints-long vector of unique strings, each identifying the corresponding rows of Data. (E.g., the filename from which that point came.)
    I::Array{Array{Int64}}   # Vector of index vectors. Each element, for example, I[i] is a vector, whose elements in turn are integers, indicating the row numbers in Data that correspond to the "set i" points
    # --- stuff about graphics handling of the plots:
    axisHandles::Array{PyCall.PyObject}  # handles to the plotted axes
    axisDims::Array{Int64}               # naxes-by-2. Each row indicates the x- and y-axes dimension for each axis plot. E.g., if row j contains [2 3] that means that the jth scatterplot has column 2 of Data as its x-axis and column 3 of Data as its y-axis
    dotHandles::Array{PyCall.PyObject}   # naxes-by-ndots, handles to the extra, initially invisible, dots plotted
    callback::Any                        # A function that could be called after a button press.
                                         # Function should be defined as taking two parameters, callback(str, SD::scatterdata).
                                         # str will be one of the strings in stringIDs

"""
type scatter_data
    Data::Array{Float64}     # npoints-by-ndims Array{Float64} of data
    stringIDs::Array{String} # npoints-long vector of unique strings, identifying the corresponding rows of Data. (E.g., the filename from which that point came.)
    I::Array{Array{Int64}}   # Vector of index vectors. I[i] is a vector, containing row numbers of "set i" points
    # --- stuff about graphics handling of the plots:
    axisHandles::Array{PyCall.PyObject}  # handles to the plotted axes
    axisDims::Array{Int64}               # naxes-by-2, indicating x- and y-axes dimension for each axis plot
    dotHandles::Array{PyCall.PyObject}   # naxes-by-ndots, extra dots plotted
    callback::Any                        # A function to be called after a button press
end

"""
SD = interactive_scatters(Data, stringIDs; set_indices=nothing,
    plot_set2=false, axisDims = [2 1 ; 3 1], user_callback=nothing,
    n_invisible_dots = 3, invisible_colors = ["c"; "b"; "m"],
    fignum = nothing, axisHandles = nothing, plot_colors = ["r"; "g"; "k"; "y" ; "m"],
    markersize=10, marker=".")

A multidimensional set of data points is plotted in a set of scatterplots, each of which
shows the scatterplot of one dimension against another. Each data point is identified
by a unique string.  When the user clicks on one of the axes, the closest point to it
is identified, and if defined, a callback function is called, with the string ID that
datapoint passed as one of its parameters.

The points can be divided into different subsets of them, with each set plotted
in its own color. In addition, a series of initially invisible dots can also be
added to the plot. If scatter_highlight() is defined as the callback function, then
these initially invidible points will become sequentially visible as the user clicks
on the scatterplots.

The main function that is called to set things up is `interactive_scatter()`. The
usual callback is `scatter_highlight()`.

This function, `interactive_scatter()`: Given a set of multidimensional data points,
puts up at most two scatterplots of different dimensions
against each other. (If data has only two dimensions, only one scatterplot goes up.)
In addition, enables GUI interactivity: users can associate a callback function
with buttonclicks on the plots. Also puts up some invisible points on the plot with can later
be used by `scatter_highlight()`.  Returns a data structure with info as to what is plotted where,
meant to be used by GUI callbacks and functions such as `scatter_highlight()`.

If the user clicks on one of the plots, the plotted data point closest to the clicked point
will be identified, and if the user callback was defined, then user_callback will be called
as user_callback(stringID, SD)

If desired, only a subset of the points in Data can be plotted; and multiple different subsets can be
requested to be plotted, in different colors.

# PARAMETERS:

- Data          An npoints-by-ndims Array{Float64}

- stringIDs     And npoints-long vector of unique strings, each of which will be used to
                identify the corresponding row in Data.

# OPTIONAL PARAMETERS:

- set_indices           Default is a vector with one element, which itself is 1:size(Data,1),
                        i.e., default will plot all points in Data as "set 1" points.
                        If passed, set_indices should be a vector of vectors; each
                        element should be a vectors of integers, each within 1:size(Data,2).
                        The rows of Data in set_indices[i] will be plotted with color plot_colors[i].

- plot_colors           n-row Array, with row i containing the color to be used for set i.

- markersize            Size of each point in the scatterplots

- marker                Type of each point in the scatterplots

- plot_set2             If true, dots for sets 2:end are plotted, otherwise not.

- axisDims              Array of integers, Naxes-by-2 in size. Each row indicates with dimensions
                        to show as horizontal axis (first column of axisDims) and vertical axis (second column)
                        in the corresponding axis.

- user_callback         If set, function that will be called, as `user_callback(stringID, SD)` where
                        stringID is the string for the row that corresponds to the selected point.

- fignum                If optional parameters axisHandles is not passed, then this will indicate the
                        Figure number on which to plot the scatterplots.If not passed, a new figure is created.

- axisHandles           Two element vector containing PyPlot axis handles, indicating where to put up
                        the scatterplots. If not passed, the figure is cleared and two side-by-side
                        subplots will be made.

- n_invisible_dots      Number of auxiliary, initially invisible, dots to put up

- invisible_colors      n_invisible_dots-long vector, indicating the color that each dot, when made
                        visible, should have.  The dots will be plotted with the first color on top.

# RETURNS

- SD::scatter_data      A structure, holindg information about the plot. See documentation for scatter_data


# EXAMPLE

```jldoctest
pygui(true)
npoints    = 10
Data       = randn(npoints,2)
string_IDs = map(x -> @sprintf("%d", x), 1:npoints)
remove_all_BPs()  # delete any previous click handlers, for cleanliness

SD = interactive_scatters(Data, string_IDs, fignum=20, user_callback=scatter_highlight);
```jldoctest

"""
function interactive_scatters(Data, stringIDs; set_indices=nothing,
    plot_set2=false, axisDims = [2 1 ; 3 1], user_callback=nothing,
    n_invisible_dots = 3, invisible_colors = ["c"; "b"; "m"],
    fignum = nothing, axisHandles = nothing, plot_colors = ["r"; "g"; "k"; "y" ; "m"; "r" ; "g" ; "k"],
    markersize=10, marker=".")

    @doc """
    scatter_event_callback(xy, r, linehandle, axhandle, SD::scatter_data)

    Internal function used by `interactive_scatters()` to enable GUI interactivity.
    This function is responsible for turning the position of the
    selected data point into the corresponding filename, and then
    calling the callback that was registered with `interactive_scatters()` (if any was)

    Out of the entries in SD::scatter_data, this function uses axisHandles, axisDims,
    Data, stringIDs, and callback.

    """ function scatter_event_callback(xy, r, linehandle, axhandle, SD::scatter_data)
        # @printf("xy=(%g,%g)\n", xy[1], xy[2])
        idx = nothing
        # Let's go through the axes finding our axis
        for i=1:length(SD.axisHandles)
            if axhandle == SD.axisHandles[i]
                myX = SD.axisDims[i,1]; myY = SD.axisDims[i,2]
                # and now find the index of the point at xy
                idx = find((SD.Data[:, myX].==xy[1]) .& (SD.Data[:, myY].==xy[2]))
                if length(idx)==0;
                    warn(@sprintf("scatter_event_callback: Couldn't find point (%.3f,%.3f), returning\n",
                        xy[1], xy[2]), bt=true);
                    return;
                end
                idx = idx[1]
            end
        end

        @printf("You selected the point with ID %s\n", SD.stringIDs[idx]);
        pause(0.0001)  # just to get the above printed out

        # If there is a user callback, call it:
        if SD.callback != nothing
            SD.callback(SD.stringIDs[idx], SD)
        end
    end

    # ---------------  OK, now the actual interactive_scatters() function -------------

    # Initialize a scatter_data structure
    SD = scatter_data([], [], [], [], [], [], nothing)
    SD.callback  = user_callback
    SD.stringIDs = stringIDs
    SD.Data      = Data

    # Don't try to plot dimensions we don't have:
    axisDims[find(axisDims.>size(Data,2))] = size(Data,2)
    axisDims = unique(axisDims, 1)
    nplots = axisHandles==nothing ? size(axisDims,1) : minimum((size(axisDims, 1), length(axisHandles)))

    # Default is to plot data from all rows as set1
    if set_indices == nothing
        set_indices = [1:size(Data,1)]
    end

    if length(plot_colors) < length(set_indices)
        error(sprintf("Need at least as many plot_colors (%d) as there are different groups of set_indices (%d)",
            length(plot_colors), length(set_indices)))
    end

    # Store indices in the SD structure that will be returned
    SD.I = I = set_indices

    # If we weren't given the axes, make them:
    if axisHandles == nothing
        # If we weren't given a figure, make it:
        if fignum==nothing
            fignum = figure()[:number]
        end
        # If we're making axes, clear the figure for them:
        figure(fignum); clf();
        if nplots==1
            axisHandles = [gca()]
        else
            axisHandles = [subplot(1,2,1), subplot(1,2,2)]
        end
    else
        # We were given axes, get figure number from them:
        if nplots==1; fignum = axisHandles[1][:figure][:number]
        else
            fignum = [axisHandles[1][:figure][:number], axisHandles[2][:figure][:number]]
            if fignum[1]==fignum[2]; fignum=fignum[1]; end
        end
    end
    # Store in return structure:
    SD.axisHandles  = axisHandles
    SD.axisDims     = axisDims

    # Now plot the points:
    for i=1:length(SD.axisHandles)
        safe_axes(SD.axisHandles[i])
        # Find the rows that correspond to these axes:
        myX = SD.axisDims[i,1]; myY = SD.axisDims[i,2]
        for i=length(set_indices):-1:1
            if i==1 || plot_set2
                plot(Data[set_indices[i],myX],  Data[set_indices[i],myY], ".", color=plot_colors[i],
                    markersize=markersize, marker=marker, linestyle="None")
            end
        end
        title(@sprintf("Dim %d vs %d", myY, myX))
    end

    # Add the invisible dots:
    hs =Array{PyCall.PyObject}(0, n_invisible_dots)
    for i=1:length(SD.axisHandles)
        safe_axes(SD.axisHandles[i]); xpos = mean(xlim()); ypos = mean(ylim())
        myh = []
        for j=1:n_invisible_dots
            # Last one to be plotted should be first color (so it'll go on top):
            myh = [myh; plot(xpos, ypos, ".", color=invisible_colors[end-(j-1)], marker=marker)]
        end
        hs = [hs ; reshape(myh[end:-1:1], 1, n_invisible_dots)]
    end
    for h in hs; h[:set_markersize](markersize); h[:set_visible](false); end
    SD.dotHandles = hs

    for f in fignum
        install_nearest_point_callback(figure(f), scatter_event_callback, user_data=SD)
    end

    return SD

end


"""
    scatter_highlight(stringID, SD::scatter_data)

A multidimensional set of data points is plotted in a set of scatterplots, each of which
shows the scatterplot of one dimension against another. Each data point is identified
by a unique string.  When the user clicks on one of the axes, the closest point to it
is identified, and if defined, a callback function is called, with the string ID that
datapoint passed as one of its parameters.

The points can be divided into different subsets of them, with each set plotted
in its own color. In addition, a series of initially invisible dots can also be
added to the plot. If scatter_highlight() is defined as the callback function, then
these initially invidible points will become sequentially visible as the user clicks
on the scatterplots.

The main function that is called to set things up is `interactive_scatter()`. The
usual callback is `scatter_highlight()`.

This function, 'scatter_highlight()`:
Finds the row in SD.stringIDs that equals stringID, and then for each axis in SD.axisHandles,
sets the first SD.dotHandle's x,y data to the positon of the corresponding row of SD.Data;
moves the second SD.dotHandle to where the first used to be; moves the third to where the
second used to be; and so on.

# EXAMPLE:

```jldoctest
SD = interactive_scatters(randn(10,3), map(x -> @sprintf("%d", x), 1:10), fignum=20)

scatter_highlight("1", SD)
scatter_highlight("5", SD)
scatter_highlight("9", SD)
```

"""
function scatter_highlight(stringID, SD::scatter_data)
    idx = find(SD.stringIDs .== stringID)
    if length(idx)==0; @printf("scatter_highlight: Couldn't find stringID %s, returning\n", stringID); return; end

    Data = SD.Data[idx,:]

    # Move our dots along:
    if length(SD.dotHandles) > 0
        # The dot in column X will get the coords of the dot in column X-1:
        for i=1:length(SD.axisHandles)
            for to=size(SD.dotHandles,2):-1:2
                from = to-1
                SD.dotHandles[i,to][:set_xdata](SD.dotHandles[i,from][:get_xdata]())
                SD.dotHandles[i,to][:set_ydata](SD.dotHandles[i,from][:get_ydata]())
                SD.dotHandles[i,to][:set_visible](SD.dotHandles[i,from][:get_visible]())
            end
            # And then the dot in column 1 gets the coords of the red dot closest to the clicked point:
            myX = SD.axisDims[i,1]; myY = SD.axisDims[i,2]
            SD.dotHandles[i,1][:set_xdata](Data[myX])
            SD.dotHandles[i,1][:set_ydata](Data[myY])
            SD.dotHandles[i,1][:set_visible](true)
        end
    end
    safe_axes(SD.axisHandles[end])
    legend(SD.dotHandles[end,1:3], ["current", "1 ago", "2 ago"])
    pause(0.0001)
end





# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.




"""
    Data structure for parameter and cost histograms

`histo_params()` returns one of these objects; used to manage GUI handling.

names::Array{String}                 # name for each parameter histogram
values::Array{Float64}               # matrix of values, nentries-by-length(names)
axisHandles::Array{PyCall.PyObject}  # handle to axis for each parameter, same length as names
LineHandles::Array{PyCall.PyObject}  # handle to the plots on each axis
files::Array{String}                 # nentries long string of files where data came from

"""
type histo_data
    names::Array{String}
    values::Array{Float64}
    axisHandles::Array{PyCall.PyObject}
    LineHandles::Array{PyCall.PyObject}
    files::Array{String}
end




"""
HD = histo_params(args, params, tcosts, costs, files; fignum=1, nbins=10)

Histogram each parameter. params should be nentries-by-length(args) in size.
args should be a vector of strings. tcosts and costs should be nentries in length,
and represent trainig cost, and test cost, respectively. files should be a vector
of filenames.

# OPTIONAL PARAMS

- fignum   The figure in which histograms will be plotted.

- nbins    The number of bins to use in each histogram.

"""
function histo_params(args, params, tcosts, costs, files; fignum=1, nbins=10, linewidth=3)

    pygui(true)
    figure(fignum); clf();

    HD = histo_data([], [], [], [], [])

    nparams = size(params,2)
    nrows = ceil(nparams/3)+1

    for i=1:nparams;
        HD.axisHandles = [HD.axisHandles ; subplot(nrows,3,i)];
        plt[:hist](params[:,i], nbins)
        title(args[i])

        myrow = ceil(i/3)
        if myrow < (nrows+1)/2;     axisHeightChange(0.8, lock="t")
        elseif myrow > (nrows+1)/2; axisHeightChange(0.8, lock="b")
        else                        axisHeightChange(0.8, lock="c")
        end
    end

    HD.axisHandles = [HD.axisHandles ; subplot(nrows, 2, nrows*2-1)]; axisHeightChange(0.8, lock="b");
    axisMove(0, -0.025); plt[:hist](tcosts*1000, nbins); title("training cost*1000")

    HD.axisHandles = [HD.axisHandles ; subplot(nrows, 2, nrows*2)];   axisHeightChange(0.8, lock="b");
    axisMove(0, -0.025); plt[:hist](costs*1000, nbins); title("test cost*1000")

    for ax in HD.axisHandles
        safe_axes(ax)
        h = plot([0 0 0 ; 0 0 0], [ylim()[1] ; ylim()[2]]*ones(1,3), visible=false, linewidth=linewidth)
        h[1][:set_color]("m"); h[2][:set_color]("b"); h[3][:set_color]("c")
        HD.LineHandles = [HD.LineHandles ; reshape(h[end:-1:1], 1, 3)]
    end

    args   = [args ; ["train cost" ; "test cost"]]
    params = [params tcosts*1000 costs*1000]

    HD.names  = args
    HD.values = params
    HD.files  = files

    return HD
end


"""
    histo_params(res; threshold=-0.0001, further_params...)

Wrapper that calls the other histo_params method, after first selecting for
only  runs that have a test cost less than threshold.  further_params are passed
on to the other histo_params method.

- threshold             training costs below this value are considered "successful" (red dots),
                        above it are "unsuccessful" (blue dots)

- cost_choice           String, used to indicate which cost will be used for thresholding. It
                        must be either "cost", indicating the testing cost, or "tcost", the training cost.

- fignum                The figure in which histograms will be plotted.

- nbins                 The number of bins to use in each histogram.


"""
function histo_params(res; threshold=-0.0001, cost_choice="cost", further_params...)
    args   = res["args"]
    params = res["params"]
    tcost  = res["tcost"]
    cost   = res["cost"]
    files  = res["files"]

    if cost_choice=="cost"
        I = find(cost.<threshold)
    elseif cost_choice=="tcost"
        I = find(tcost.<threshold)
    else
        error("cost_choice MUST be one of \"tcost\" or \"cost\"")
    end

    return histo_params(args, params[I,:], tcost[I], cost[I], files[I,:]; Dict(further_params)...)
end


"""
    histo_highlight(filename, HD::histo_data)

Assuming that histo_params was called, and returned the HD that is passed to this function,
this function will put up a vertical bar at the values corresponding to the run indicated
by filename.

Will put up to three bars, blue for the most recent one asked for; green for one call ago
to this function; and magenta for the time the function was called two calls ago. That allows
keeping track of a sequence of requested values.

This function can be used for GUI interactivity-- if a button click selects a run, the
corresponding filename can then be used with this function to highlight that run in the
histograms.

# EXAMPLE:

```jldoctest
res = farmload("C17", farmdir="MiniFarms")

HD = histo_params(res)
I = find(res["cost].<-0.0001)  # to match what histo_params shows

histo_highlight(res["files"][I[1]], HD)
histo_highlight(res["files"][I[2]], HD)
histo_highlight(res["files"][I[3]], HD)
```
"""
function histo_highlight(filename, HD::histo_data)
    idx = find(HD.files .== filename)
    if length(idx)==0; @printf("histo_highlight: Couldn't find filename %s, returning\n", filename); return; end

    for i=1:length(HD.names)
        for to=size(HD.LineHandles,2):-1:2
            from=to-1
            HD.LineHandles[i,to][:set_xdata](HD.LineHandles[i,from][:get_xdata]())
            HD.LineHandles[i,to][:set_ydata](HD.LineHandles[i,from][:get_ydata]())
            HD.LineHandles[i,to][:set_visible](HD.LineHandles[i,from][:get_visible]())
        end
        HD.LineHandles[i,1][:set_xdata]([HD.values[idx,i], HD.values[idx,i]])
        HD.LineHandles[i,1][:set_visible](true)
    end
end


"""
Data structure for PCA scatterplots

`plot_PCA()` returns one of these objects; used to manage GUI handling.
"""
type PCAplot_data
    # --- stuff about the data that is plotted:
    mu::Array{Float64}   # This vector of means is subtracted from a params vector
    sd::Array{Float64}   # Which are then ./ by this sd vector to get the whitened params vector
    V::Array{Float64}    # matrix of principal components, one per column
    C::Array{Float64}    # covariance matrix after whitening
    D::Array{Float64}    # Vector of eigenvalues
    I::Array{Int64}      # index into "successful" runs
    nI::Array{Int64}     # index into "unsuccessful"runs
    Vparams::Array{Float64}  # nruns-by-nparams vector of runs, in PCA space
    files::Array{String} # nruns-long vector of corresponding filenames
    # --- stuff about graphics handling of the plots:
    axisHandles::Array{PyCall.PyObject}  # handles to the plotted axes
    axisPCs::Array{Int64}    # naxes-by-2, indicating x- and y-axes PC # for each axis plot
    dotHandles::Array{PyCall.PyObject} # naxes-by-ndots, extra dots plotted
    BP::PyCall.PyObject # the kbMonitorModule.kb_monitor that will keep track of button presses
    callback::Any   # The function to be called after a button press
end





"""
    PCA_highlight(filename, PC::PCAplot_data)

Assuming that `PCA_plot()` was called, and returned the PC that is passed to this function,
this function will put up a colored dot, across all axes, at the values corresponding to the run
indicated by filename. [Out of the entries in PC::PCAplot_data, this function uses files,
Vparams, axHandles, dothandles, and axisPCs]

Will put up to three dots, blue for the most recent one asked for; green for one call ago
to this function; and magenta for the time the function was called two calls ago. That allows
keeping track of a sequence of requested values.

This function can be used for GUI interactivity-- if a button click selects a run, the
corresponding filename can then be used with this function to highlight that run in the
scatterplots.

# EXAMPLE:

```jldoctest
res = farmload("C17", farmdir="MiniFarms")

PC = plot_PCA(res)
I = find(res["tcost].<-0.0002)  # to match what histo_params shows

PCA_highlight(res["files"][I[1]], PC)
PCA_highlight(res["files"][I[2]], PC)
PCA_highlight(res["files"][I[3]], PC)
```
"""
function PCA_highlight(fname, PC::PCAplot_data)
    idx = find(PC.files .== fname)
    if length(idx)==0; @printf("PCA_highlight: Couldn't find filename %s in the PCAplot_data structure, returning\n", fname); return; end

    Vparams = PC.Vparams[idx,:]

    # Move our non-red dots along:
    if length(PC.dotHandles) > 0
        # The dot in column X will get the coords of the dot in column X-1:
        for i=1:length(PC.axisHandles)
            for to=size(PC.dotHandles,2):-1:2
                from = to-1
                PC.dotHandles[i,to][:set_xdata](PC.dotHandles[i,from][:get_xdata]())
                PC.dotHandles[i,to][:set_ydata](PC.dotHandles[i,from][:get_ydata]())
                PC.dotHandles[i,to][:set_visible](PC.dotHandles[i,from][:get_visible]())
            end
            # And then the dot in column 1 gets the coords of the red dot closest to the clicked point:
            myX = PC.axisPCs[i,1]; myY = PC.axisPCs[i,2]
            PC.dotHandles[i,1][:set_xdata](Vparams[end-(myX-1)])
            PC.dotHandles[i,1][:set_ydata](Vparams[end-(myY-1)])
            PC.dotHandles[i,1][:set_visible](true)
        end
    end
    safe_axes(PC.axisHandles[end])
    legend(PC.dotHandles[end,1:3], ["current", "1 ago", "2 ago"])
    pause(0.001)
end


"""
    PCA_event_callback(xy, r, linehandle, axhandle, PC::PCAplot_data)

Internal function used by `PCA_plot()` to enable GUI interactivity.
This function is responsible for turning the position of the
selected data point into the corresponding filename, and then
calling the callback that was registered with `PCA_plot()` (if any was)

Out of the entries in PC::PCAplot_data, this function uses axisHandles, axisPCs,
Vparams, files, and callback.

"""
function PCA_event_callback(xy, r, linehandle, axhandle, PC::PCAplot_data)
    # @printf("xy=(%g,%g)\n", xy[1], xy[2])
    idx = nothing
    # Let's go through the axes finding our axis
    for i=1:length(PC.axisHandles)
        if axhandle == PC.axisHandles[i]
            myX = PC.axisPCs[i,1]; myY = PC.axisPCs[i,2]
            # and now find the index of the point at xy
            idx = find((PC.Vparams[:, end-(myX-1)].==xy[1]) .& (PC.Vparams[:, end-(myY-1)].==xy[2]))
            if length(idx)==0;
                @printf("PCA_event_callback: Couldn't find point (%.3f,%.3f), returning\n", xy[1], xy[2]);
                return;
            end
            idx = idx[1]
        end
    end

    @printf("You selected file %s\n", PC.files[idx]); pause(0.0001)

    # If there is a user callback, call it:
    if PC.callback != nothing
        PC.callback(PC.files[idx], PC)
    end
end





"""
    PC = plot_PCA(res; threshold=-0.0002, fignum=2, pc_offset=0, plot_unsuccessful=true,
        compute_good_only=true, unsuccessful_threshold=nothing, cost_choice="cost",
        user_callback=nothing)

Computes the principal components for a matrix of parameter values across many runs, and puts
up several scatterplots of the run parameter values projected onto these principal
components. In addition, enables GUI interactivity: users can associate a callback function
with buttonclicks on the figure.

# PARAMETERS:

- res     A dictionary, in the format of the output of `farmload()`

# OPTIONAL PARAMETERS:

- threshold             costs below this are considered "successful" runs

- unsuccessful_threshold   costs above this are considered "UNsuccessful" runs
                        This value defaults to whatever threshold is.  If the two thresholds are
                        different, farms with costs in between, in the "no man's land", are not plotted.

- plot_unsuccessful     If true, dots for unsuccessful runs are shown, otherwise not.

- compute_good_only   If true, all runs, whether successful or not, are used to compute the PCs.
                        If false, only the successful runs.

- cost_choice           String, used to indicate which cost will be used for thresholding. It
                        must be either "cost", indicating the testing cost, or "tcost", the training cost.

- pc_offset             Shows PCs from pc_offset+1 to pc_offset+4

- user_callback         If set, function that will be called, as `usercallback(filename, PC)` where
                        filename is the name of the run that corresponds to the selected point.

- fignum                Figure number on which to plot the PC scatterplots

"""
function plot_PCA(res; threshold=-0.0001, fignum=2, pc_offset=0, plot_unsuccessful=true,
    compute_good_only=true, unsuccessful_threshold=nothing, cost_choice="cost",
    user_callback=nothing)

    if unsuccessful_threshold==nothing
        unsuccessful_threshold = threshold
    end

    # Initialize a PCAplot_data structure
    PC = PCAplot_data([], [], [], [], [], [], [], [], [], [], [], [], nothing, nothing)

    # Get the runs' costs and do selection on them:
    if ~((cost_choice=="tcost")  ||  (cost_choice=="cost"))
        error("cost_choice MUST be one of \"tcost\" or \"cost\"")
    end
    mycost = res[cost_choice];
    PC.I  = I  = find(mycost .< threshold)
    PC.nI = nI = find(mycost .>= unsuccessful_threshold)

    if compute_good_only
        # Use the successful runs (sparams) to define PCA space
        sparams = copy(res["params"][I,:])
        for i=1:size(sparams,2)
            sparams[:,i] -= mean(sparams[:,i])
            sparams[:,i] /= std(sparams[:,i])
        end
        C = sparams'*sparams/size(sparams,1)
        D,V = eig(C)
        pv = 100*D/sum(D)
    end

    # Then z-score everybody
    params = copy(res["params"]); nruns = size(params,1)
    PC.mu = mean(params,1);
    params = params - ones(nruns,1)*PC.mu
    PC.sd = std(params,1);
    params = params ./ (ones(nruns,1)*PC.sd)

    if !compute_good_only
        C = params'*params/nruns
        D,V = eig(C)
        pv = 100*D/sum(D)
    end
    PC.V = V

    # parameters in the eigen-coords:
    PC.Vparams = Vparams = (inv(V)*params')'

    figure(fignum); clf();
    # ax1 = subplot(2,2,1); axisHeightChange(0.9, lock="t")
    # ax2 = subplot(2,2,2); axisMove(0.05, 0); axisHeightChange(0.9, lock="t")
    # ax3 = subplot(2,2,3); axisHeightChange(0.9, lock="b")
    # ax4 = subplot(2,2,4); axisMove(0.05, 0); axisHeightChange(0.9, lock="b")
    # PC.axisHandles = [ax1, ax2, ax3, ax4]
    # PC.axisPCs = [1 2 ; 3 2 ; 1 3 ; 3 4] + pc_offset
    ax1 = subplot(1,2,1); ax2 = subplot(1,2,2)
    PC.axisHandles = [ax1 ; ax2]
    PC.axisPCs = [2 1 ; 3 1]

    for i=1:length(PC.axisHandles)
        safe_axes(PC.axisHandles[i])
        myX = PC.axisPCs[i,1]; myY = PC.axisPCs[i,2]
        if plot_unsuccessful
            plot(Vparams[nI,end-(myX-1)],  Vparams[nI,end-(myY-1)], "g.", markersize=10)
        end
        plot(Vparams[I,end-(myX-1)],  Vparams[I,end-(myY-1)], "r.", markersize=10)
        title(@sprintf("PCA %d (%.2f%%) vs %d (%.2f%%)", myY, pv[end-(myY-1)], myX, pv[end-(myX-1)]))
    end


    # Add the non-red dots:
    hs =Array{PyCall.PyObject}(0, 3)
    for i=1:length(PC.axisHandles)
        safe_axes(PC.axisHandles[i]);
        hs = [hs ; reshape([plot(0, 0, "m.") plot(0, 0, "b.") plot(0, 0, "c.")][end:-1:1], 1, 3)]
        # order of handles gets reversed so the last one plotted -- goes on top -- is first handle
    end
    for h in hs; h[:set_markersize](11); h[:set_visible](false); end
    PC.dotHandles = hs
    legend(hs[end,1:3], ["current", "1 ago", "2 ago"])
    PC.callback = user_callback
    PC.files = res["files"]

    # safe_axes(ax3)
    # if plot_unsuccessful; legend(["unsuccessful", "successful"])
    # else                   legend(["successful"])
    # end

    install_nearest_point_callback(figure(fignum), PCA_event_callback, user_data=PC)

    return PC
end



"""
SV = plot_SVD(;threshold =-0.0002, plot_unsuccessful=false, compute_good_only=false,
        cost_choice="cost", user_callback=nothing, fignum=100)

Uses the pre-computed SVD components for a matrix of parameter values across many runs, and puts
up scatterplots of the run parameter values projected onto these SV components (columns of the U
matrix. In addition, enables GUI interactivity: users can associate a callback function
with buttonclicks on the figure.


# OPTIONAL PARAMETERS:

- threshold             training costs below this value are considered "successful" (red dots),
                        above it are "unsuccessful" (blue dots)

- plot_unsuccessful     If true, dots for unsuccessful runs are shown, otherwise not.

- compute_good_only     If true, only the successful runs are used to compute the SVD space

- cost_choice           String, used to indicate which cost will be used for thresholding. It
                        must be either "cost", indicating the testing cost, or "tcost", the training cost.

- user_callback         If set, function that will be called, as `usercallback(filename, SV)` where
                        filename is the name of the run that corresponds to the selected point.

- fignum                Figure number on which to plot the PC scatterplots


# RETURNS

- SV     A structure of type PCAplot_data. (Only its files, Vparams, axHandles, dotHandles
         axisPCs and callback entries will have assigned values.)


"""
function plot_SVD(;threshold =-0.0001, plot_unsuccessful=false, compute_good_only=false,
    cost_choice="cost", user_callback=nothing, fignum=100)

    # get response matrix
    response, results = load("MiniOptimizedC17_SVD_response_matrix3.jld", "response","results");
    # set up filter by nan
    nanrows = any(isnan.(response),2);

    # Get the runs' costs and do selection on them:
    if ~((cost_choice=="tcost")  ||  (cost_choice=="cost"))
        error("cost_choice MUST be one of \"tcost\" or \"cost\"")
    end
    mycost = results[cost_choice];
    if !compute_good_only
        mycost = mycost[.!vec(nanrows),:];
        disp_cost = copy(mycost);
    end
    badcost = mycost .>= threshold;

    # if we are computing SVD only on the good farms, update nanrows and disp_cost
    if compute_good_only
        nanrows = nanrows .| badcost;
        disp_cost = mycost[.!vec(nanrows),:];
    end

    # Filter response matrix
    r_all = response[.!vec(nanrows),:];
    m = mean(r_all,1);
    r_all = r_all - repmat(m, size(r_all,1),1);
    F = svdfact(r_all);
    u = copy(F[:U]);
    u1 = u[:,1];
    u2 = u[:,2];
    u3 = u[:,3];

    # Make list of just good farms
    if compute_good_only
        # nanrows filter already selects for good farms
        u1good = u1;
        u2good = u2;
        u3good = u3;
    else
        # remove bad cost farms
        u1good = u1[.!vec(badcost),:];
        u2good = u2[.!vec(badcost),:];
        u3good = u3[.!vec(badcost),:];
    end
    files = results["files"];
    files = files[.!vec(nanrows),:];

    SV = PCAplot_data([], [], [], [], [], [], [], [], [], [], [], [], nothing, nothing)
    SV.files = files
    SV.Vparams = u[:,3:-1:1]  # The column vectors are expected in REVERSE order (i.e., as returned by eig(). )

    pygui(true)
    figure(fignum); clf();

    ax1 = subplot(1,2,1)
    if plot_unsuccessful; plot(u[:,3],u[:,1],"bo"); end
    plot(u3good, u1good, "ro")
    title("SVD U columns 3 and 1")
    ylabel("SVD Dim 1")
    xlabel("SVD Dim 3")
    plot(0, 0, "go")[1][:set_visible](false)

    ax2 = subplot(1,2,2)
    if plot_unsuccessful; plot(u[:,2],u[:,1],"bo"); end
    plot(u2good, u1good, "ro")
    title("SVD U columns 2 and 1")
    xlabel("SVD Dim 2")
    remove_ytick_labels()

    SV.axisHandles = [ax1 ; ax2]
    SV.axisPCs   = [3 1 ; 2 1]

    hs =Array{PyCall.PyObject}(0, 3)
    for i=1:length(SV.axisHandles)
        safe_axes(SV.axisHandles[i]);
        hs = [hs ; reshape([plot(0, 0, "m.") plot(0, 0, "b.") plot(0, 0, "c.")][end:-1:1], 1, 3)]
        # order of handles gets reversed so the last one plotted -- goes on top -- is first handle
    end
    for h in hs; h[:set_markersize](11); h[:set_visible](false); end
    SV.dotHandles = hs
    legend(hs[end,1:3], ["current", "1 ago", "2 ago"])

    SV.callback = user_callback
    # [Out of the entries in PC::PCAplot_data, this function uses files,  Vparams, axHandles, dothandles, and axisPCs]

    install_nearest_point_callback(figure(fignum), PCA_event_callback, user_data=SV)

    return SV
end


# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.



plot_farm_trials = 10    #  The number of trials to be plotted per farm run



"""
    params = plot_farm(filename; testruns=400, setup_file=nothing, fignum=3,
        plottables = ["V", "V[1,:]-V[4,:]"], ylabels=["V", "ProR-ProL"],
        ylims = [[-0.02, 1.02], [-1.02, 1.02]], plot_list = [1:20;],
        hit_linestyle="-", err_linestyle="--", xlims=nothing,
        do_plot = true, do_print  = true,
        overrideDict=Dict(), further_params...)

    Plots multiple trials from a single run of a farm. ASSUMES that three different
    types of conditions are run, control, delay opto, and choice opto; and within
    each of these, Pro trials and Anti trials are run.

# PARAMETERS

- filename    Eithe a String, the filename of the .jld file containing the run,
to be loaded; OR an Array{Float64} vector, containing the parameter values.
In the later case, a setup_file must be defined, to set the other parameters
(see optional parameters below).

# OPTIONAL PARAMETERS

- setup_file  if filename is actaually a parameter vector, then this optional parameter
              should be passed as a string pointing to a .jld file that when loaded, will contain
              variables "mypars", "extra_pars", and "search_conditions" -- the names of the arguments
              corresponding to the different entries of the parameter vector will be taken to be
              the string versions of the keys of the dictionary "search_conditions".

- testruns    Number of trials to run. Defaults to value of global variable plot_farm_trials.

- fignum      Figure to put the plot up in.

- overrideDict   A dictionary containing any model parameter values that will
              override any values loaded from the file.  For example
              `overrideDict = Dict(:sigma=>0.001)` will run with that value
              of sigma, no whater what the file said.

- plottables   A vector of strings. Each of these strings indicates something to
            plot; the strings will be evaluated in a context where the variables "V", "U",
            and "t" are instantiated to have the values passed to `plot_PA()`. Thus
            plottables=["V", "V[1,:]-V[4,:]"] indicates that two axes should be used, on
            the first one V will be plotted, on the second V[1,:]-V[4,:].
                The strings can be any arbitrary Julia expression, with the condition that
            their evaluation should produce something with either nrows or ncolumns
            equal to length(t).

- ylabels   A vector of strings, should be equal in length to plottables; Each element
            here will be used as the y label for the corresponding axis.

- ylims     A vector of Any, should be in length to plottables. If an element is
            `nothing`, then allow automatic scaling of the y axes for this axis. Otherwise,
            the element should be a 2-long vector of Float64, indicating the minimum and
            the maximum for the y axis, respectively.

- xlims     Either nothing (autoscale) or a 2-long vector of Float64, to be applied
            to all axes

- hit_linestyle  A string indicating the linestyle for hit trials. E.g., "-" or "--".
               If this is passed as the empty string, "", then hits are not plotted.

- err_linestyle  A string indicating the linestyle for error trials. E.g., "-" or "--".
               If this is passed as the empty string, "", then errors are not plotted.

- further_params    Any further keyword value params are passed on to run_ntrials.
            Note that unilke entries in overrideDict, which take the highest precedence,
            keyword-value pairs here take the lowest-precedence: any kw-val pair that also
            appears in the farm's .jld file, or in overrideDict, will be ignored in favor
            of those higher precedence instances.

- do_plot   If false, no plots of graphics calls are generated

- do_print  If false, no printing of param values to the cosole is done.

- pstrings  Vector of string titles for the different opto conditions in
            extra_pars[:opto_periods]. Must be at least as long as there are rows
            in extra_pars[:opto_periods]


# RETURNS

- params     The parameters of the farm that was plotted.

- hBP_out    A vector with fraction of Pro hits; each element corresponds to one opto_period condition

- hBA_out    A vector with fraction of Pro hits; each element corresponds to one opto_period condition

- PV         nperiods-4-testruns matrix of final unit voltages in Pro trials

- AV         nperiods-4-testruns matrix of final unit voltages in Anti trials

- fPV        nperiods-4-nsteps-testruns matrix of unit voltages as a function of time in Pro trials

- fAV        nperiods-4-nsteps-testruns patrix of unit voltages as a function of time in Anti trials

# EXAMPLE CALL

```jldoctest
plot_farm("MiniOptimized/farm_C17_Farms024_0058.jld", fignum=200, testruns=20, setup_file="Setups/setup_C21.jld", plot_list=[1:5;],
    plottables = ["V", "V[1,:]+V[4,:] - (V[2,:]+V[3,:])"], ylabels=["V", "Rule encoding"],
    ylims = [[-0.02, 1.02], nothing],
    overrideDict=Dict(:dt=>0.0025));

```jldoctest
"""
function plot_farm(filename; testruns=400, setup_file=nothing, fignum=3,
    plottables = ["V", "V[1,:]-V[4,:]"], ylabels=["V", "ProR-ProL"],
    ylims = [[-0.02, 1.02], [-1.02, 1.02]], plot_list = [1:20;],
    hit_linestyle="-", err_linestyle="--", xlims=nothing,
    do_plot=true, do_print=true,
    pstrings = ["CONTROL", "DELAY OPTO", "CHOICE OPTO"],
    overrideDict=Dict(), further_params...)

    mypars=extra_pars=args=pars3=[]  # define to be available outside if block
    if typeof(filename)==String
        mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")
    else
        if setup_file==nothing
            error("If filename is not a string, need a setup file that, when loaded, defines mypars, extra_pars, and search_conditions. E.g., setup_file=\"Setups/setup_C20.jld\"")
        end
        if ! (typeof(filename)<:Array{Float64})
            error("If filename is not a string, it should be a vector of Float64s")
        end
        mypars, extra_pars, search_conditions=load(setup_file, "mypars", "extra_pars", "search_conditions")
        args = map(k -> String(k), keys(search_conditions))
        pars3 = filename
        if length(args) != length(pars3)
            error("If filename is not a string, it should be a vector with length equal to the number of keys in setup_file")
        end
    end

    if do_plot
        pygui(true)
        figure(fignum); clf();
    else
        plot_list = []
    end

    nperiods = size(extra_pars[:opto_periods],1)
    if nperiods > length(pstrings)
        error("need pstrings to be as long as there are rows in extra_pars[:opto_periods]")
    end
    hBP_out = zeros(1, nperiods)
    hBA_out = zeros(1, nperiods)
    PV_out = AV_out = fPV_out = fAV_out = []

    for period = 1:nperiods
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
            :opto_times=>reshape(extra_pars[:opto_periods][period,:], 1, 2),
        ))
        these_pars = merge(Dict(further_params), these_pars)

        # The plot_list should be the one we give it below, not whatever was in the stored parameters
        delete!(these_pars, :plot_list)

        if do_plot
            nplots  = length(plottables)
            pax_set = Array{PyCall.PyObject}(nplots,1)
            aax_set = Array{PyCall.PyObject}(nplots,1)
            for i=1:nplots
                pax_set[i] = subplot(2*nplots, 3, period + 3*(i-1));
                if i<=nplots/2; axisHeightChange(0.9, lock="t")
                else            axisHeightChange(0.9, lock="c")
                end

                aax_set[i] = subplot(2*nplots, 3, period + 3*(nplots+i-1));
                if i<=nplots/2; axisHeightChange(0.9, lock="c")
                else            axisHeightChange(0.9, lock="b")
                end
            end
        else
            pax_set   = []
            aax_set   = []
            plot_list = []
        end

        proVs, antiVs, pro_fullV, anti_fullV = run_ntrials(testruns, testruns; plot_list=plot_list,
            ax_set = [pax_set, aax_set],
            plottables = plottables, ylabels=ylabels, ylims=ylims, xlims=xlims,
            hit_linestyle=hit_linestyle, err_linestyle=err_linestyle,
            merge(make_dict(args, pars3, these_pars), overrideDict)...);
        hBP_out[period] = hBP = length(find(proVs[1,:]  .> proVs[4,:])) /size(proVs, 2)
        hBA_out[period] = hBA = length(find(antiVs[4,:] .> antiVs[1,:]))/size(antiVs,2)
        # @printf("period %d:  hBP=%.2f%%, hBA=%.2f%%\n\n", period, 100*hBP, 100*hBA)

        if length(PV_out) == 0
            PV_out = zeros(nperiods, size(proVs,1), size(proVs,2))
            AV_out = zeros(nperiods, size(proVs,1), size(proVs,2))
            fPV_out = zeros(nperiods, size(pro_fullV,1),  size(pro_fullV, 2),  size(pro_fullV,3))
            fAV_out = zeros(nperiods, size(anti_fullV,1), size(anti_fullV, 2), size(anti_fullV,3))
        end
        PV_out[period,:,:]    = proVs
        AV_out[period,:,:]    = antiVs
        fPV_out[period,:,:,:] = pro_fullV
        fAV_out[period,:,:,:] = anti_fullV

        if do_plot
            safe_axes(pax_set[1]); title(@sprintf("%s  PRO hits = %.2f%%", pstrings[period], 100*hBP))
            safe_axes(aax_set[1]); title(@sprintf("ANTI hits = %.2f%%", 100*hBA))
            for i=1:(nplots-1)
                safe_axes(pax_set[i]); remove_xtick_labels(); xlabel("")
                safe_axes(aax_set[i]); remove_xtick_labels(); xlabel("")
            end
            safe_axes(pax_set[end]); remove_xtick_labels(); xlabel("")
            if period > 1
                remove_ytick_labels(pax_set)
                remove_ytick_labels(aax_set)
            end

            if period==2 && typeof(filename)<:String
                safe_axes(pax_set[1])
                text(mean(xlim()), ylim()[2] + 0.35*(ylim()[2]-ylim()[1]), filename,
                    fontsize=18, horizontalalignment="center")
            end
            figure(fignum)[:canvas][:draw]()
            pause(0.0001)
        end
    end

    if do_print
        for a=1:length(args)
            myarg = args[a]; while length(myarg)<20; myarg=myarg*" "; end
            @printf("%s\t\t%g\n", myarg, pars3[a])
        end
    end

    return pars3, hBP_out, hBA_out, PV_out, AV_out, fPV_out, fAV_out
end


# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.


"""
make_mini_farm(farmid ; fromdirs=["../Farms024", "../Farms025", "../Farms026"],
    todir="MiniFarms")

Takes all the  runs in directories fromdirs, (which might not be on git),
and puts a small-size sumamry of them in todir. That todir directory will
a reasonabale size for git (only MBytes compred to GBytes)

The minifarm directory can be used for browsers in lieu of the original farms, but
note that files in it do not include Hessian information.

# PARAMETERS:

- farmid     A pattern that needs to be matched in a filename for it to be included. E.g., "C17"

# OPTIONAL PARAMETERS:

- fromdirs  Either a String, indicating a directory, or a vector of Strings, indicating multiple
            directories to be treated together.

- todir     A String indicating the directory where the Mini files should go to. If the
            directory did not previously exist, creates it. If the directory was not previously
            empty, does not clear it.

# EXAMPLE:

```jldoctest
make_mini_farm("C17", fromdirs=["../Farms024], todir="MiniNew")
```jldoctest

"""
function make_mini_farm(farmid; fromdirs=["../Farms024", "../Farms025", "../Farms026"],
        todir="MiniFarms")

    res = farmload(farmid, verbose=true, farmdir=fromdirs)

    if ~isdir(todir); mkdir(todir); end;

    sdict = []; dirname=[]; filename=[]

    for i in 1:length(res["tcost"])
        mypars, extra_pars, args, pars3 = load(res["files"][i], "mypars", "extra_pars",
            "args", "pars3")
        sdict = Dict("mypars"=>mypars, "extra_pars"=>extra_pars,
            "args"=>args, "pars3"=>pars3)
        for k in keys(res)
            if k=="tcost"
                sdict["traj3"] = [0 ; res["tcost"][i]; 0]
            elseif !(k=="args" || k=="params")
                sdict[k] = res[k][i]
            end
        end

        dirname, filename = splitdir(res["files"][i]); dirname=splitdir(dirname)[2]

        save(todir * "/" * filename[1:9]*dirname*"_"*filename[10:end], sdict)
        if rem(i, 20)==0
            @printf("Did %d/%d\n", i, length(res["tcost"]))
        end
    end
end




# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.


"""
    make_maxi_farm(farmid, fromdirs, todir)

Takes the runs in a passed set of farm directories and copies each file dir/fname*farmid* into
newdir/fname*farmid*_dir  so that all files live in one directory but don't have names that
step on each other. The files to be copied HAVE to start with "farm_", all others are ignored

"""
function make_maxi_farm(farm_id, fromdirs, todir)

    if ~isdir(todir); mkdir(todir); end;
    if typeof(fromdirs)==String; fromdirs=[fromdirs]; end

    for d in fromdirs
        fromname = split(d, "/")
        while fromname[1]==".." || fromname[1]==""; fromname=fromname[2:end]; end
        while fromname[end]==""; fromname=fromname[1:end-1]; end
        fromname = join(fromname)
        n = length("farm_"*farm_id*"_")
        for f in filter(x -> startswith(x, "farm_" * farm_id * "_"), readdir(d))
            toname = todir*"/farm_"*farm_id*"_"*fromname*"_"*f[n+1:end]
            @printf("cp(%s, %s)\n", d*"/"*f, toname)
            cp(d*"/"*f, toname, remove_destination=true)
        end
    end
end
