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
"""
function setup_control_figure(fignum)
    figure(fignum); clf();
    
    # rad: plot farm versus not
    ax1 = axes([0.01, 0.84, 0.37, 0.15])
    rad = kbMonitorModule.radio_buttons(ax1, ["Don't plot trials", "Plot trials (wait for it)"])

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

    return rad, heb, plb, ylb, slb, dbx, tbx, xlb
end



# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.


#######################################################
#
#     Actual calls to the functions to run everything
#
#######################################################

using MultivariateStats
using Clustering

farmid  = "C17"
farmdir = "MiniOptimized" ; nclusters = 3
# farmdir = "../C17_Optimized_B"; nclusters = 5;

if !isdefined(:histo_params)
    @printf("Loading results_analysis.jl\n")
    include("results_analysis.jl")
end

if !isdefined(:res) || true
    res = farmload(farmid, verbose=true, farmdir=farmdir)
    res["params"][:,4] = abs.(res["params"][:,4])
end

pygui(true); 
remove_all_BPs(); # Clean up any previous links between clicks on figures and callback functions
plt[:close](1); plt[:close](2); plt[:close](3); plt[:close](4); plt[:close](5)

# Carlos' favored configuration, but adjust to suit -- use capture_current_figure_configuration() 
# to see code that reproduces a configuration you like once you find it

figure(1); set_current_fig_position(1325, 41, 640, 672)   # x, y, width, height
figure(2); set_current_fig_position(645, 785, 680, 408)
figure(3); set_current_fig_position(0, 785, 641, 407)
figure(4); set_current_fig_position(3, 23, 1288, 797)
figure(5); set_current_fig_position(1338, 736, 540, 450)

# -------------- SET UP CONTROL FIGURE --------

rad, heb, plb, ylb, slb, dbx, tbx, xlb = setup_control_figure(5)

##########################################################
#
#   next couple of lines define cost type and threshold
#
##########################################################

cost_choice = "cost"   # "cost" is test cost, "tcost" is training cost
threshold   = -0.0002  # Below this is a "successful" run


nsucc = length(find(res[cost_choice].<threshold))
@printf("\n%d runs being plotted as successful, with cost `%s' less than %g\n\n", nsucc, cost_choice, threshold)

# Put up the histograms
HD = histo_params(res; threshold=threshold, cost_choice=cost_choice);
pause(0.001)

# The callback function that will be called after clicking on a data dot:
function highlight_all(fname, LD, SV)
    PCA_highlight(fname, SV);       # Color the selected dots in the SV plot
    scatter_highlight(fname, LD);   # Color the selected dots in the LDA plot
    histo_highlight(fname, HD)      # Color the selected bars in the histograms
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


# Put up the LDA plot
use_marino = false   # if true, loads the cluster labels and linear discriminant projections from the output of Marino's matlab code
if use_marino
    G = matread("compute_clustering/LDA_output.mat")
    idx = find(G[cost_choice].<=threshold)
    G["cluster_ids"] = G["cluster_ids"][idx]
    G["ldparams"]    = G["ldparams"][:,idx]
    G["files"]       = G["files"][idx]
    set_indices = [find(G["cluster_ids"].==1), find(G["cluster_ids"].==2), find(G["cluster_ids"].==3)]
    LD = interactive_scatters(G["ldparams"]', G["files"], set_indices=set_indices, plot_set2=true,
        user_callback=(fname, Trash) -> highlight_all(fname, LD, SV), fignum=2, axisDims = [1 2]);
else # compute cluster labels and linear discriminant projections ourselves
    G1 = load(farmdir*farmid*"_SVD_response_matrix3.jld")   # these include the trial-averaged PSTHS Alex produced
    idx = find(G1["results"][cost_choice].<=threshold)       # choose runs
    dynamics_data = G1["response"][idx,:]                    # get the PSTHS
    fnames = G1["results"]["files"][idx]                     #    and their corresponding filenames
    params = G1["results"]["params"][idx,:]                  #    and their corresponding parameter values

    cluster_ids = assignments(kmeans(dynamics_data', nclusters, init=:kmcen))   # label each run with a cluster id
    
    M = fit(SubspaceLDA, params', nclusters, cluster_ids)            # do multiclass LDA in param space with those cluster ids
    ld_data = transform(M, params')'                         # get the parameters in the LDA projection

    # Now put it on the screen
    set_indices = Array{Any}(nclusters)
    for i=1:nclusters
        set_indices[i] = find(cluster_ids.==i)
    end
    if nclusters==3; axisDims=[1 2]; else axisDims=[1 2 ; 3 4]; end
    LD = interactive_scatters(ld_data, fnames, set_indices=set_indices, plot_set2=true,
        user_callback=(fname, Trash) -> highlight_all(fname, LD, SV), 
        fignum=2, axisDims = axisDims);
    
end
pause(0.001)

# Put up the SVD plot
SV = plot_SVD(threshold=threshold, cost_choice=cost_choice, 
    user_callback = (fname, Trash) -> highlight_all(fname, LD, SV),
    plot_unsuccessful=false, compute_good_only=false, fignum=3);


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



