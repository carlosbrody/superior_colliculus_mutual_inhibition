# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.


"""
    rad, heb, plb, ylb, slb, dbx, tbx, xlb = setup_control_figure(fignum)

Configures figure fignum to be a control figure for LDA_browser.

# RETURNS:

- rad     Radio buttons controlling whether to plot trials or not

- heb     All/Hits only/Errors only buttons

- plb     Plottables boxes. Non-empty entries here will correspond to an axis
per condition when plotting trials. Text inside these boxes can be arbitrary Julia
expressions, evaluated in a context where the following variables exist:

    - t    time vector

    - V    4-by-ntimesteps matrix of V activations

    - U    4-by-ntimesteps matrix of U activations

    - rule  timesteps-long vector, equal to 0.5*(V[1,:]+V[4,:]-V[2,:]-V[3,:]) (Pro minus Anti)

    - decis timesteps-long vector, equal to V[1,:]-V4,:] (Pro_R minus Pro_L)

    - diag  timesteps-long vector, equal to 0.5*(V[1,:]+V[3,:]-V[2,:]-V[4,:]) (difference of two diagonals)

    So for example you could put in "rule.^2" if for whatever reason you wanted to see 
    the rule squared as a function of time

- ylb     Ylim boxes for the plottables

    Put in "[min, max]" if you want to enforce certain ylimits. Leave empty for auto-scaling.

- slb   selectize box. Arbitrary Julia Boolean expressions, evaluated in a context with
        the same variables as for the plottables boxes, can be put in here. Only trials for
        which the expression is true will be plotted. So, for example, you could write:

        "t0=tbin(t, 1.1); abs.(decis[t0])>0.4"

        to see only trials in which the absolute value of the decision signal at time t=1.1 
        would be greater than 0.4

- dbx     Override dictionary box. Julia Dict-style entries can go here to override
        parameter values in the model. For example, to force running with sigma=0.1 and dt=0.005,
        no matter that the particular parameters in the file to be plotted were, you would
        write:

        :sigma=>0.1, :dt=>0.005

        

- tbx     number of trials to run box

- xlb     xlims box. Write in "[min, max]" if you want particular horizontal axis limits
          for the various plots.

- npla    Number of "plottables" boxes

"""
function setup_control_figure(fignum)
    figure(fignum); clf();
    
    # rad: plot farm versus not
    ax1 = axes([0.01, 0.84, 0.37, 0.15])
    rad = kbMonitorModule.radio_buttons(ax1, ["Plot trials (wait for it)", "Don't plot trials"])

    # heb : plot all or hits only or err only
    ax0 = axes([0.5, 0.8, 0.37, 0.19])
    heb = kbMonitorModule.radio_buttons(ax0, ["all", "hit only", "err only"])

    # plb, ylb :  the content of plottables, and their y limits
    npla = 3; axheight = 0.1; axy = 0.02; 
    plax = Array{PyCall.PyObject}(npla,1); ylax = Array{PyCall.PyObject}(npla,1)
    plb  = Array{PyCall.PyObject}(npla,1); ylb  = Array{PyCall.PyObject}(npla,1)
    for i=1:npla
        plax[i] = axes([0.13, axy, 0.49, axheight])
        plb[i]  = kbMonitorModule.text_box(plax[i], "plottable", "")
        ylax[i] = axes([0.75, axy, 0.24, axheight])
        ylb[i]  = kbMonitorModule.text_box(ylax[i], "ylim", "")

        axy = axy+axheight+0.02
    end
    # Default is to plot V, rule encoding, and Pro_R - Pro_L:
    plb[end][:set_val](" V ");               ylb[end][:set_val]("[-0.02, 1.02] ")
    plb[end-1][:set_val](" rule "); ylb[end-1][:set_val]("[-1.02, 1.02] ")
    plb[end-2][:set_val](" decis "); ylb[end-2][:set_val]("[-1.02, 1.02] ")

    # Julia expression to select individual trials for plotting
    axy += 0.02
    ax6 = axes([0.13, axy, 0.86, axheight])
    slb = kbMonitorModule.text_box(ax6, "Selectize", "")
    # ax7 = axes([0.01, axy, 0.1, axheight])
    # scb = kbMonitorModule.check_buttons(ax7, "", [true])

    axy += axheight + 0.02
    ax3 = axes([0.13, axy, 0.86, axheight])
    dbx = kbMonitorModule.text_box(ax3, "Override ", "")

    axy += axheight + 0.03
    ax4 = axes([0.2, axy, 0.3, axheight])
    tbx = kbMonitorModule.text_box(ax4, "ntrials to run ", "6")

    ax5 = axes([0.68, axy, 0.3, axheight])   # working on adding xlims control
    xlb = kbMonitorModule.text_box(ax5, "xlims", "")

    return rad, heb, plb, ylb, slb, dbx, tbx, xlb, npla
end



# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.


#######################################################
#
#     Actual calls to the functions to run everything
#
#######################################################

using MultivariateStats
using Clustering

if !isdefined(:histo_params)
    @printf("Loading results_analysis.jl\n")
    include("results_analysis.jl")
end

"""
    LDA_browser(farmid, farmdir, threshold=-0.0002, cost_choice="cost", 
        default_nclusters=4)

Runs a browser that puts up histograms, and scatterplots of the different clusters 
in both parameter space (figure 2) and in the original SVD dynamic space (figure 3).

Uses Marino's code to guess the number of clusters, and Alex's code to define 
the dynamics space in which those clusters are found. See top of 
`Results Analysis.ipynb` for a description of how to run Alex and Marino's codes.

Currently (2pm EST 19-Dec-2017), these have been run on three farms of files,
leading to the following three possible calls here:

```jldoctest
LDA_browser("C17", "C17_Optimized_A")
    or
LDA_browser("C17", "Mini_C17_Optimized_B")
    or
LDA_browser("C19", "Mini_C19")
```

All three run on the same parameter sets. The first two were generated identically. 
The last one, C19, was generated by first doing a stringent criterion after optimizing
with only 50 runs/condition, then if the test was passed, further optimizing with 
1600 runs/condition.

# PARAMETERS:

- farmid   The filename pattern common to runs in the desired farm, e.g., "C17"

- farmdir  The directory in which the farm's runs are found, e.g., "Mini_C19"

# OPTIONAL PARAMETERS:

- threshold     Any runs with a cost of type cost_choice (see below) greater 
                than this are dropped from the analysis

- cost_choice    String, one of "cost" (test cost) or "tcost" (training cost)

- default_nclusters   There will be an attempt to read the number of clusters 
                off of a file produced by Marino's code. If that file doesn't exist, this
                will be the number of clusters used.

# FURTHER INFO

See the documentation for `setup_control_figure()` for information
on how to control what is plotted. See (and modify) first few lines of
this function `LDA_browser()` to change the default figure placements.

"""
function LDA_browser(farmid, farmdir, threshold=-0.0002, cost_choice="cost", 
    default_nclusters=4)
# farmid  = "C17"; farmdir = "Mini_C17_Optimized_B"
# farmid  = "C19"; farmdir = "Mini_C19"
# farmid  = "C17"; farmdir = "C17_Optimized_A"


    if !isdefined(:res) || true
        res = farmload(farmid, verbose=true, farmdir=farmdir)
        res["params"][:,4] = abs.(res["params"][:,4])
    end

    pygui(true); 
    remove_all_BPs(); # Clean up any previous links between clicks on figures and callback functions
    plt[:close](1); plt[:close](2); plt[:close](3); plt[:close](4); plt[:close](5)

    # Carlos' favored configuration, but adjust to suit -- 
    # use capture_current_figure_configuration() 
    # to see code that reproduces a configuration you like once you find it

    figure(1); set_current_fig_position(1325, 41, 640, 672)   # x, y, width, height
    figure(2); set_current_fig_position(645, 785, 680, 408)
    figure(3); set_current_fig_position(0, 785, 641, 407)
    figure(4); set_current_fig_position(3, 23, 1288, 797)
    figure(5); set_current_fig_position(1338, 736, 540, 450)

    # -------------- SET UP CONTROL FIGURE --------

    rad, heb, plb, ylb, slb, dbx, tbx, xlb, npla = setup_control_figure(5)

    nsucc = length(find(res[cost_choice].<threshold))
    @printf("\n%d runs being plotted as successful, with cost `%s' less than %g\n\n", nsucc, cost_choice, threshold)

    # Put up the histograms
    HD = histo_params(res; threshold=threshold, cost_choice=cost_choice);
    pause(0.001)

    # The callback function that will be called after clicking on a data dot:
    function highlight_all(fname, SV_space, PARAM_space)
        scatter_highlight(fname, SV_space);      # Color the selected dots in the LDA plot
        scatter_highlight(fname, PARAM_space);   # Color the selected dots in the LDA plot
        histo_highlight(fname, HD)               # Color the selected bars in the histograms
        pause(0.001);   # We don't really care about the 1 ms pause; just a convenient way to flush all pending graphics commandsj

        if rad[:value_selected] == "Plot trials (wait for it)"
            if ~isnull(tryparse(Int64, tbx[:text]));               ntrials = parse(Int64, tbx[:text])
            else; @printf("Couldn't parse the ntrials to plot\n"); ntrials = 10
            end
            plottables = Array{String}(0,1); ylims = Array{Any}(0,1); ylabels = Array{String}(0,1)
            for i=npla:-1:1
                if plb[i][:text] != ""
                    plottables = [plottables ; plb[i][:text]]
                    ylabels    = [ylabels    ; plb[i][:text]]
                    if ylb[i][:text] == ""
                        ylims = [ylims ; nothing]
                    else
                        ylims = [ylims ; [eval(parse(ylb[i][:text]))]]
                    end                
                end
            end
            xlims = xlb[:text]=="" ? nothing : eval(parse(xlb[:text]))
            if heb[:value_selected]     == "all";      hstyle="-"; estyle="--";
            elseif heb[:value_selected] == "hit only"; hstyle="-"; estyle="";
            else                                       hstyle="";  estyle="--";
            end
            selectize = slb[:text]=="" ? "true" : slb[:text]
            try
                plot_farm(fname, testruns=ntrials, fignum=4, 
                    plottables=plottables, ylabels=ylabels, ylims=ylims, xlims=xlims,
                    hit_linestyle = hstyle, err_linestyle = estyle, selectize=selectize,
                    overrideDict=eval(parse("Dict("*dbx[:text]*")"))) 
            catch e
                @printf("Couldn't plot farm, error %s\n", e)
                catch_stacktrace()
            end
        end
    end


    # -----------------  Put up the LDA plot  ---------
    # We'll use the output from Marino's Matlab code to get number of clusters.
    # The rest happens here.-
    nclusters = default_nclusters
    marino_file = "compute_clustering/MarinoCode_"*farmid*"_"*farmdir*".jld"
    if isfile(marino_file)
        nclusters, cost_choice_marino, threshold_marino = 
            load(marino_file, "ngroups", "cost_choice", "threshold")
        nclusters=Int64(nclusters)
        if threshold_marino != threshold || cost_choice_marino != cost_choice
            warn("Marino file not generated with requested threshold and cost choice!\n", 
            "Number of clusters may be off.\n")
        end
    else
        warn(@sprintf("Didn't find file %s with output from Marino's code, guessing nclusters=%d\n",
        marino_file, nclusters))
    end


    alexfile = farmdir*"_"*farmid*"_SVD_response_matrix3_reduced.jld"
    results, response= load(alexfile, "results", "response")   # these include the trial-averaged PSTHS Alex produced
    idx       = find((results[cost_choice].<threshold)  .&  .!any(isnan.(response),2));
    files     = results["files"][idx]
    params    = results["params"][idx,:]
    # Hack params so that sigma is always positive
    args      = load(files[1], "args")
    params[:, find(args.=="sigma")] = abs.(params[:, find(args.=="sigma")])
    # Zero-mean the responses
    response  = response[idx,:]
    response  = response - ones(size(response,1),1)*mean(response,1)

    U, S, V = svd(response)

    # We arbitrarily use the top 4 principal components:
    ncomponents=4
    D = response*V[:,1:ncomponents]
    cluster_ids = assignments(kmeans(D', nclusters, init=:kmcen)) # , display=:iter))
    set_indices = Array{Any}(nclusters)
    for i=1:nclusters; set_indices[i] = find(cluster_ids.==i); end


    # Scatterplot in singular value space
    SV_space = interactive_scatters(D, files, set_indices=set_indices, plot_set2=true, 
        fignum=3, axisDims=[1 2], 
        user_callback=(fname, Trash) -> highlight_all(fname, SV_space, PARAM_space));

    # Now do LDA in parameter space, and put up scatter plot in parameter space
    M = fit(SubspaceLDA, params', nclusters, cluster_ids) # do multiclass LDA in param space with those cluster ids
    ld_data = transform(M, params')'                         # get the parameters in the LDA projection
    PARAM_space = interactive_scatters(ld_data, files, set_indices=set_indices, plot_set2=true,
        fignum=2, axisDims = [2 1; 3 1],
        user_callback=(fname, Trash) -> highlight_all(fname, SV_space, PARAM_space));



    # Print out some instructions for the user
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
end


