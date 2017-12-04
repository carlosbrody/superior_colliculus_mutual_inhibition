# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.


#######################################################
#
#     Actual calls to the functions to run everything
#
#######################################################

if !isdefined(histo_params)
    error("You must call include(\"results_analysis.jl\") before calling include(\"C17_browser.jl\"))
end

if !isdefined(:res)
    res = farmload("C17", verbose=true, farmdir="MiniOptimized")
    res["params"][:,4] = abs.(res["params"][:,4])
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
tbx = dbx = []  # define these outside the try/catch so the vars are available outside the try/catch
try 
    tbx = kbMonitorModule.text_box(ax2, "ntrials to run ", "10")
    dbx = kbMonitorModule.text_box(ax3, "Override ", "")
catch
    @printf("\nI couldn't make the ntrials to run and Override text boxes for you.\n")
    tbx = Dict(:text=>"10")
    dbx = Dict(:text=>"")
end


##########################################################
#
#   next couple of lines define cost type and threshold
#
##########################################################

cost_choice = "cost"   # "cost" is test cost, "tcost" is training cost
threshold   = -0.00025  # Below this is a "successful" run


nsucc = length(find(res[cost_choice].<threshold))
@printf("\n%d runs being plotted as successful, with cost `%s' less than %g\n\n", nsucc, cost_choice, threshold)

# Put up the histograms
HD = histo_params(res; threshold=threshold, cost_choice=cost_choice);
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
            @printf("Couldn't plot farm, error %s\n", e)
        end
    end
end

# Put up the PCA plot
PC = plot_PCA(res; threshold=threshold, cost_choice=cost_choice, 
    user_callback = (fname, Trash) -> highlight_all(fname, PC, SV),
    plot_unsuccessful=false, unsuccessful_threshold=0.0001, compute_good_only=true, fignum=2);
pause(0.001)

# Put up the SVD plot
SV = plot_SVD(threshold=threshold, cost_choice=cost_choice, 
    user_callback = (fname, Trash) -> highlight_all(fname, PC, SV),
    plot_unsuccessful=false, compute_good_only=true, fignum=3);

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



# DON'T MODIFY THIS FILE -- the source is in file Results Analysis.ipynb. Look there for further documentation and examples of running the code.


"""
    make_mini_farm()

Takes the C17 runs in ../Farms024, ../Farms025, and ../Farms026, which are not on git,
and puts a small-size sumamry of them in "MiniFarms/". That MiniFarms directory is
a reasonabale size for git (only 63MBytes compred to GBytes)

The MiniFarms directory can be used for the C17_browser in lieu of the original farms.

"""
function make_min_farm(;fromdirs=["../Farms024", "../Farms025", "../Farms026"], todir="MiniFarms")
    
    res = farmload("C17", verbose=true, farmdir=fromdirs)

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

        save("MiniFarms/"*filename[1:9]*dirname*"_"*filename[10:end], sdict)
        if rem(i, 20)==0
            @printf("Did %d/%d\n", i, length(res["tcost"]))
        end
    end
end




