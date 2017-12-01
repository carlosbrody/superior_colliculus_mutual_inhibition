# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.



include("pro_anti.jl")

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

- farmdir      Directory in which the files are found

- verbose      If true, print a brief progress report every verbose_every files loaded

- verbose_every  Used if verbose==true (see above)

# RETURNS:

A dictionary with the keys "tcost" (vector or training costs), "cost" (vector of test costs), 
"dirs" (vector directories in which files  are found), "files" (vector of filenames, including 
paths), "qu_out" (deprecated for farm C17 and other non-pre-search farms), "params"" (nruns-by-
nparams matrix of paramters), "args" (nparams-long vector of parameter names).

Note that most values are nruns-long vectors but "args" is only nparams-long, because they are
the same names across all the runs.
"""
function farmload(farm_id; farmdir="../NewFarms", verbose=true, verbose_every=50)

    if typeof(farmdir)==String; farmdir=[farmdir]; end
    
    results = Dict(); dirs=[]; files =[]; qs=[]; tcosts=[]; costs=[]; pars=[]; n=0;
    for dd in farmdir
        for f in filter(x -> startswith(x, "farm_" * farm_id * "_"), readdir(dd * "/"))
            n += 1
            myfile = dd * "/" * f; 
            fj = jldopen(myfile); 
            if exists(fj, "qu_out"); qu_out = load(myfile, "qu_out"); 
            else;                    qu_out = [-1 -1]; 
            end
            close(fj)

            args, params, traj3, cost = load(myfile, "args", "pars3", "traj3", "cost")

            if     !haskey(results, args);       results["args"] = args;
            elseif !all(results["args"].==args); error("Not all files have same args!"); 
            end

            files = [files ; myfile]; dirs = [dirs ; dd]
            qs = [qs ; qu_out]; tcosts = [tcosts; traj3[2,end]]; costs = [costs; cost]
            if length(pars)==0; pars=params'; else pars = [pars ; params']; end
            if verbose && rem(n, verbose_every)==0
                @printf("%s %g\n", myfile, tcosts[end])
            end
        end
    end

    results["dirs"] = dirs
    results["files"]  = files
    results["qu_out"] = qs
    results["tcost"]  = tcosts
    results["cost"]   = costs
    results["params"] = pars

    return results
end




"""
    myres = selectize(res, I)

Given res, the output of `farmload()` and a vector of indices into an nruns-long
vector, returns a version of res in which only those indices are present in the
values -- i.e., vectors will now be length(I)-long, not nruns-long.
"""
function selectize(res, I)
    myres = Dict()
    for k in keys(res)
        if length(res[k])>1 && !(typeof(res[k])<:Array{String})
            myres[k] = res[k][I,:]
        else
            myres[k] = res[k]
        end
    end
    return myres
end


"""
    Data structure for parameter and cost histograms

`histo_params()` returns one of these objects; used to manage GUI handling.
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
        axes(ax)
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
    if length(idx)==0; @printf("PCA_highlight: Couldn't find filename %s, returning\n", fname); return; end

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
    axes(PC.axisHandles[end])    
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
        axes(PC.axisHandles[i])
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
        axes(PC.axisHandles[i]); 
        hs = [hs ; reshape([plot(0, 0, "m.") plot(0, 0, "b.") plot(0, 0, "c.")][end:-1:1], 1, 3)]
        # order of handles gets reversed so the last one plotted -- goes on top -- is first handle
    end
    for h in hs; h[:set_markersize](11); h[:set_visible](false); end
    PC.dotHandles = hs
    legend(hs[end,1:3], ["current", "1 ago", "2 ago"])
    PC.callback = user_callback
    PC.files = res["files"]
        
    # axes(ax3)
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
    response, results = load("SVD_response_matrix.jld", "response","results");
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
        axes(SV.axisHandles[i]); 
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




plot_farm_trials = 10    #  The number of trials to be plotted per farm run

"""
    params = plot_farm(filename; testruns=nothing, fignum=3, overrideDict=Dict())

    Plots multiple trials from a single run of a farm.

# PARAMETERS

- filename    The filename of the .jld file containing the run, to be loaded

# OPTIONAL PARAMETERS

- testruns    Number of trials to run. Defaults to value of global variable plot_farm_trials.

- fignum      Figure to put the plot up in.

- overrideDict   A dictionary containing any model parameter values that will
              override any values loaded from the file.  For example
              `overrideDict = Dict(:sigma=>0.001)` will run with that value
              of sigma, no whater what the file said.

"""
function plot_farm(filename; testruns=400, fignum=3, overrideDict=Dict())

    mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")


    pygui(true)
    figure(fignum); clf();
    
    pstrings = ["CONTROL", "DELAY OPTO", "CHOICE OPTO"]
    for period = 1:3
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
        # :opto_strength=>0.3, 
        :opto_times=>reshape(extra_pars[:opto_periods][period,:], 1, 2),
        # :opto_times=>["target_start-0.4" "target_start"],
        # :opto_times=>["target_start" "target_end"],
        # :post_target_period=>0.3,
        # :rule_and_delay_period=>1.2,
        # :dt=>0.005,
        ))

        # The plot_list should be the one we give it below, not whatever was in the stored parameters
        delete!(these_pars, :plot_list)

        pvax = subplot(4,3,period);   axisHeightChange(0.9, lock="t")
        pdax = subplot(4,3,period+3); axisHeightChange(0.9, lock="c"); 
        avax = subplot(4,3,period+6); axisHeightChange(0.9, lock="c")
        adax = subplot(4,3,period+9); axisHeightChange(0.9, lock="b")

        proVs, antiVs = run_ntrials(testruns, testruns; plot_list=[1:20;], plot_Us=false, 
            ax_set = Dict("pro_Vax"=>pvax, "pro_Dax"=>pdax, "anti_Vax"=>avax, "anti_Dax"=>adax),
        merge(make_dict(args, pars3, these_pars), overrideDict)...);

        hBP = length(find(proVs[1,:]  .> proVs[4,:])) /size(proVs, 2)
        hBA = length(find(antiVs[4,:] .> antiVs[1,:]))/size(antiVs,2)
        # @printf("period %d:  hBP=%.2f%%, hBA=%.2f%%\n\n", period, 100*hBP, 100*hBA)

        axes(pvax); title(@sprintf("%s  PRO hits = %.2f%%", pstrings[period], 100*hBP))
        axes(avax); title(@sprintf("ANTI hits = %.2f%%", 100*hBA))
        axes(pdax); remove_xtick_labels(); xlabel("")
        if period > 1
            remove_ytick_labels([pvax, pdax, avax, adax])
        end
        
        figure(fignum)[:canvas][:draw]()
        pause(0.0001)
    end

    for a=1:length(args)
        myarg = args[a]; while length(myarg)<20; myarg=myarg*" "; end
        @printf("%s\t\t%g\n", myarg, pars3[a])
    end

    return pars3
end



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


#######################################################
#
#     Actual calls to the functions to run everything
#
#######################################################



if !isdefined(:res)
    res = farmload("C17", verbose=true, farmdir="MiniFarms")
end

pygui(true); 
remove_all_BPs(); # Clean up any previous links between clicks on figures and callback functions
plt[:close](1); plt[:close](2); plt[:close](3); plt[:close](4); plt[:close](5)

# Carlos' favored configuration, but adjust to suit -- use capture_current_figure_configuration() 
# to see code that reproduces a configuration you like once you find it

figure(1); set_current_fig_position(1325, 41, 640, 982)   # x, y, width, height
figure(2); set_current_fig_position(645, 785, 680, 408)
figure(3); set_current_fig_position(0, 785, 641, 407)
figure(4); set_current_fig_position(3, 23, 1288, 797)
figure(5); set_current_fig_position(1338, 998, 540, 200)

figure(5)
ax1 = subplot(2,2,1); axisMove(-0.05, 0); axisWidthChange(1.1, lock="l")
ax2 = subplot(2,2,2); axisMove(0.05, 0); axisWidthChange(0.6, lock="r")
ax3 = subplot(2,1,2); axisMove(0.05, 0)
rad = kbMonitorModule.radio_buttons(ax1, ["Don't plot trials", "Plot trials (wait for it)"])
tbx = kbMonitorModule.text_box(ax2, "ntrials to run ", "10")
dbx = kbMonitorModule.text_box(ax3, "Override ", "")


HD = histo_params(res; threshold=-0.0001, cost_choice="cost");
pause(0.001)

# The callback function that will be called after clicking on a data dot:
function highlight_all(fname, PC, SV)
    PCA_highlight(fname, PC);   # Color the selected dots in the PC plot
    PCA_highlight(fname, SV);   # Color the selected dots in the SVD plot
    histo_highlight(fname, HD)  # Color the selected bars in the histograms
    pause(0.001);   # We don't really care about the 1 ms pause; just a convenient way to flush all pending graphics commandsj
    if rad[:value_selected] == "Plot trials (wait for it)"
        if ~isnull(tryparse(Int64, tbx[:text])); 
            ntrials = parse(Int64, tbx[:text])
        else
            @printf("Couldn't parse the ntrials to plot\n")
            ntrials = 10
        end
        try
            plot_farm(fname, testruns=ntrials, fignum=4, overrideDict=eval(parse("Dict("*dbx[:text]*")"))) 
        catch e
            @printf("Couldn't plot farm, error %s", e)
        end
    end
end

PC = plot_PCA(res; threshold=-0.0001, cost_choice="cost", 
    user_callback = (fname, Trash) -> highlight_all(fname, PC, SV),
    plot_unsuccessful=false, unsuccessful_threshold=0.0001, compute_good_only=true, fignum=2);
pause(0.001)

SV = plot_SVD(threshold=-0.0001, cost_choice="cost", 
    user_callback = (fname, Trash) -> highlight_all(fname, PC, SV),
    plot_unsuccessful=false, compute_good_only=true, fignum=3);

docstring = """

Wait for PCA and SVD plots to come up.
Then click on any dot within the PCA plot or SVD plots to
see the corresponding data in the other plots. Click on
'Don't plot trials' in figure 5 if you want to go faster,
without running trials.

The window placement fits a 15-in Macbook Pro, but adjust
at will. Once you find window positions you like, run
    capture_current_figure_configuration()
to get copy-pastable code that reproduces it.

ntrials to run is the number of trials per condition that will
be run in plot_farm() to compute %correct. Up to 20 of them
will be plotted.

Override is a dict that will override any paramater values in
the plotted trials. For example, if you write
    :sigma=>0.2, :rule_and_delay_period=>1.6
then no matter what the run's file said, if you ask for trials
to be plotted, those trials will run
with that sigma and that rule_and_delay_period.

"""

@printf("%s", docstring)



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


"""
    make_mini_farm()

Takes the C17 runs in ../Farms024, ../Farms025, and ../Farms026, which are not on git,
and puts a small-size sumamry of them in "MiniFarms/". That MiniFarms directory is
a reasonabale size for git (only 63MBytes compred to GBytes)

The MiniFarms directory can be used for the C17_browser in lieu of the original farms.

"""
function make_min_farm()
    
    res = farmload("C17", verbose=true, farmdir=["../Farms024", "../Farms025", "../Farms026"])

    if ~isdir("MiniFarms"); mkdir("MiniFarms"); end;

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

        save("MiniFarms/"*filename[1:9]*dirname*"_"*filename[10:end], sdict)
        if rem(i, 20)==0
            @printf("Did %d/%d\n", i, length(res["tcost"]))
        end
    end
end




