# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


include("pro_anti.jl")

using PyCall
# The following line is PyCall-ese for "add the current directory to the Python path"
unshift!(PyVector(pyimport("sys")["path"]), "")
# We use Python to enable callbacks from the figures.
@pyimport kbMonitorModule


# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


using HDF5

"""
    results = farmload(farm_id; farmdir="../NewFarms", verbose=true, verbose_every=10)

Package and return a summary of results from a number of runs.
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
    histo_params(args, params; fignum=1, nbins=10)

Histogram each parameter. params should be nentries-by-length(args) in size.
args should be a vector of strings
"""
function histo_params(args, params; fignum=1, nbins=10)

    pygui(true)
    figure(fignum); clf();

    nparams = size(params,2)
    nrows = ceil(nparams/3)

    for i=1:nparams;
        ax = subplot(nrows,3,i); 
        plt[:hist](params[:,i], nbins)
        title(args[i])
        
        myrow = ceil(i/3)
        if myrow < (nrows+1)/2;     axisHeightChange(0.8, lock="t")
        elseif myrow > (nrows+1)/2; axisHeightChange(0.8, lock="b")
        else                        axisHeightChange(0.8, lock="c")
        end
    end
end


"""
    histo_params(res; threshold=-0.0002, further_params...)
"""
function histo_params(res; threshold=-0.0002, further_params...)
    args   = res["args"]
    params = res["params"]
    tcost  = res["tcost"]
    
    I = find(tcost.<threshold)
    
    histo_params(args, params[I,:]; further_params...)
end


"""
    plot_PCA(res; threshold=-0.0002, fignum=2)

"""
function plot_PCA(res; threshold=-0.0002, fignum=2, pc_offset=0, plot_unsuccessful=true)

    tcost = res["tcost"]
    I = find(tcost .< threshold)
    nI = find(tcost .>= threshold)

    # Use the successful runs (sparams) to define PCA space
    sparams = copy(res["params"][I,:])
    for i=1:size(sparams,2)
        sparams[:,i] -= mean(sparams[:,i])
        sparams[:,i] /= std(sparams[:,i])
    end
    C = sparams'*sparams/size(sparams,1)
    D,V = eig(C)
    pv = 100*D/sum(D)

    # Then z-score everybody
    params = copy(res["params"])
    for i=1:size(params,2)
        params[:,i] -= mean(params[:,i])
        params[:,i] /= std(params[:,i])
    end

    # parameters in the eigen-coords:
    Vparams = (inv(V)*params')'

    figure(fignum); clf(); 
    subplot(2,2,1); axisHeightChange(0.9, lock="t")
    if plot_unsuccessful
        plot(Vparams[nI,end-pc_offset],  Vparams[nI,end-2-pc_offset], "g.")
    end
    plot(Vparams[I,end-pc_offset],   Vparams[I,end-2-pc_offset], "r.")
    xlabel(@sprintf("PCA %d (%.2f%%)", pc_offset+1, pv[end-pc_offset]))
    ylabel(@sprintf("PCA %d (%.2f%%)", pc_offset+3, pv[end-2-pc_offset]))
    # legend(["unsuccessful", "successful"])
    # ylim(-3.5, 3)

    subplot(2,2,2); axisMove(0.05, 0); axisHeightChange(0.9, lock="t")
    if plot_unsuccessful
        plot(Vparams[nI,end-1-pc_offset], Vparams[nI,end-3-pc_offset], "g.")
    end
    plot(Vparams[I,end-1-pc_offset], Vparams[I,end-3-pc_offset], "r.")
    # ylim(-4, 4)
    # xlim(-4, 3.5)
    xlabel(@sprintf("PCA %d (%.2f%%)", pc_offset+2, pv[end-1-pc_offset]))
    ylabel(@sprintf("PCA %d (%.2f%%)", pc_offset+4, pv[end-3-pc_offset]))

    subplot(2,2,3); axisHeightChange(0.9, lock="b")
    if plot_unsuccessful
        plot(Vparams[nI,end-pc_offset],  Vparams[nI,end-1-pc_offset], "g.")
    end
    plot(Vparams[I,end-pc_offset],   Vparams[I,end-1-pc_offset], "r.")
    xlabel(@sprintf("PCA %d (%.2f%%)", pc_offset+1, pv[end-pc_offset]))
    ylabel(@sprintf("PCA %d (%.2f%%)", pc_offset+2, pv[end-1-pc_offset]))
    legend(["unsuccessful", "successful"])
    # ylim(-3.5, 3)

    subplot(2,2,4); axisMove(0.05, 0); axisHeightChange(0.9, lock="b")
    if plot_unsuccessful
        plot(Vparams[nI,end-2-pc_offset], Vparams[nI,end-3-pc_offset], "g.")
    end
    plot(Vparams[I,end-2-pc_offset], Vparams[I,end-3-pc_offset], "r.")
    # ylim(-4, 4)
    # xlim(-4, 3.5)
    xlabel(@sprintf("PCA %d (%.2f%%)", pc_offset+3, pv[end-2-pc_offset]))
    ylabel(@sprintf("PCA %d (%.2f%%)", pc_offset+4, pv[end-3-pc_offset]))

    return C, V, D, Vparams, I, nI 
end

"""
    myres = selectize(res, I)
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


# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


"""
    params = plot_farm(filename; testruns=400, fignum=3, overrideDict=Dict())
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
    end

    for a=1:length(args)
        myarg = args[a]; while length(myarg)<20; myarg=myarg*" "; end
        @printf("%s\t\t%g\n", myarg, pars3[a])
    end

    return pars3
end



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


# --- Load the data, plot histograms, and show the PCA scatter

# The following line 
# res = farmload("C17", verbose=true, farmdir=["../Farms024" "../Farms025" "../Farms026"])
res = farmload("C17", verbose=true, farmdir="MiniFarms")

threshold = -0.0002

histo_params(res; threshold=threshold, nbins=10)
I  = find(res["tcost"].<threshold)
nI = find(res["tcost"] .>= threshold)

@printf("There were %d/%d (%.2f%%) successful runs\n", length(I), length(res["tcost"]), 
100*length(I)/length(res["tcost"]))


C, V, D, Vparams, I, nI  = plot_PCA(res; threshold=threshold, 
    pc_offset=0, plot_unsuccessful=false);



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


# --- now run the data browser ---

function restart_figure2()
    C, V, D, Vparams, I, nI  = plot_PCA(res; threshold=threshold, 
        pc_offset=0, plot_unsuccessful=false);

    files = res["files"]
    mu = mean(res["params"],1)
    sd = std( res["params"],1)


    # Set up the default blue dot    
    hs = []; 
    figure(2)
    ax1 = subplot(2,2,1); hs = [hs ; [plot(0, 0, "b.") plot(0, 0, "g.") plot(0, 0, "m.")]]
    ax2 = subplot(2,2,2); hs = [hs ; [plot(0, 0, "b.") plot(0, 0, "g.") plot(0, 0, "m.")]]
    ax3 = subplot(2,2,3); hs = [hs ; [plot(0, 0, "b.") plot(0, 0, "g.") plot(0, 0, "m.")]]
    ax4 = subplot(2,2,4); hs = [hs ; [plot(0, 0, "b.") plot(0, 0, "g.") plot(0, 0, "m.")]]


    BP = []

    function event_callback()
        bpe = BP[:buttonlist]()
        if length(bpe)>0 && bpe[1][1] != nothing        

            if bpe[1][1]==ax1
                J = (Vparams[I,end]   - bpe[1][2]).^2 + (Vparams[I,end-2] - bpe[1][3]).^2            
            elseif bpe[1][1]==ax2
                J = (Vparams[I,end-1] - bpe[1][2]).^2 + (Vparams[I,end-3] - bpe[1][3]).^2
            elseif bpe[1][1]==ax3
                J = (Vparams[I,end]   - bpe[1][2]).^2 + (Vparams[I,end-1] - bpe[1][3]).^2
            elseif bpe[1][1]==ax4
                J = (Vparams[I,end-2] - bpe[1][2]).^2 + (Vparams[I,end-3] - bpe[1][3]).^2
            else
                error("Which axes did this buttonpress come from???")
            end
            idx = indmin(J)

            @printf("---- FILE %s: -----\n", files[I[idx]])
            pars = plot_farm(files[I[idx]], testruns=20, fignum=3)
            # pars = load(files[I[idx]], "pars3")
            pars = (pars-mu')./sd'
            vpars = inv(V)*pars
            # print("\n"); print(vpars); print("\n")

            for row=1:4
                for from=2:-1:1
                    to = from+1
                    hs[row,to][:set_xdata](hs[row,from][:get_xdata]())
                    hs[row,to][:set_ydata](hs[row,from][:get_ydata]())
                end
            end

            hs[1,1][:set_xdata](vpars[end  ]); hs[1,1][:set_ydata](vpars[end-2])
            hs[2,1][:set_xdata](vpars[end-1]); hs[2,1][:set_ydata](vpars[end-3])
            hs[3,1][:set_xdata](vpars[end  ]); hs[3,1][:set_ydata](vpars[end-1])
            hs[4,1][:set_xdata](vpars[end-2]); hs[4,1][:set_ydata](vpars[end-3])

            BP[:clear_buttonlist]()
            figure(2)
        end
    end


    BP = kbMonitorModule.kb_monitor(figure(2), callback=event_callback)
    legend(hs[4,1:3], ["current", "1 ago", "2 ago"])
end

@printf("\n\n***Type restart_figure2() into the Julia prompt to start having fun***\n\n")
@printf("Clicking on any red dot in figure 2 will bring up 20 trials of the\n")
@printf("corresponding run on figure 3.\n\n")
@printf("The blue dot is the current farm displayed; green is one click ago;\n")
@printf("and magenta is two clicks ago.\n\n")


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




