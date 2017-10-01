# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


include("rate_networks.jl")  # that will also include genera_utils.jl, constrained_parabolic_minimization.jl, and hessian_utils.jl



# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


"""
    plot_PA(t, U, V; fignum=1, clearfig=true, rule_and_delay_period=1, target_period=1, post_target_period=1,
        other_unused_params...)

Helper function for plotting ProAnti results. Given a time axia and nunits-by-nsteps U and V matrices, will 
plot them (first clearing the figure if clearfig=true but overlaying if clearfig=false) in three vertically-
arranged subplots. The top one has V as a function of time, with the two Pro units in blues and the two Anti 
units in reds, dark colors for the left side of the brain, light colors for the right.  The middle subplot
has Us. And the bottom subplot shows the difference between the two Pro units as a function of time.

"""
function plot_PA(t, U, V; fignum=1, clearfig=true, rule_and_delay_period=1, target_period=1, post_target_period=1,
    other_unused_params...)
    figure(fignum)
    if clearfig; clf(); end
    
    ax1 = subplot(3,1,1)
    h = plot(t, V'); 
    setp(h[1], color=[0, 0, 1])
    setp(h[2], color=[1, 0, 0])
    setp(h[3], color=[1, 0.5, 0.5])
    setp(h[4], color=[0, 1, 1])
    ylabel("V")

    ax = gca()
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([V[:];oldlims[1]])-0.1, maximum([V[:];oldlims[2]])+0.1)
    yl = [ylim()[1], ylim()[2]]
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            -0.05, 1.05, linewidth=2)
    if yl[1]<0.02
        yl[1] = -0.02
    end
    if yl[2]>0.98
        yl[2] = 1.02
    end
    ylim(yl)
    grid(true)
    remove_xtick_labels(ax1)
        
    ax2 = subplot(3,1,2)
    hu = plot(t, U')
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([U[:];oldlims[1]])-0.1, maximum([U[:];oldlims[2]])+0.1)
    setp(hu[1], color=[0, 0, 1])
    setp(hu[2], color=[1, 0, 0])
    setp(hu[3], color=[1, 0.5, 0.5])
    setp(hu[4], color=[0, 1, 1])
    ylabel("U"); 
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            ylim()[1], ylim()[2], linewidth=2)
    remove_xtick_labels(ax2)

    grid(true)
    
    subplot(3,1,3)
    delta = V[1,:] - V[4,:]
    hr = plot(t, delta)
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([delta[:];oldlims[1]])-0.1, maximum([delta[:];oldlims[2]])+0.1)
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            ylim()[1], ylim()[2], linewidth=2)
    xlabel("t"); ylabel("Pro R - Pro L")
    grid(true)
        
end




# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


"""
parsed_times = parse_opto_times(opto_times, model_params)

This function takes period specified with symbols and turns them into times in secs. Each element in the 
opto_times parameter can be a number (time in secs) or a string. This string will be parsed and
evaluated. Four special symbols are available: trial_start (which is really zero), target_start, 
target_end, and trial_end.  To compute these, IT IS ASSUMED that model_params will contain entries
for :rule_and_delay_period, :target period, and :post_target_period

In addition, for backwards compatibility with Alex's code, some numbers are special:
20 means end of trial (but we prefer the string "trial_end"), 100 mens end of rule and delay (but we prefer
"target_start", and 200 means end of target (we prfer "target_end").

# PARAMETERS

* opto_times, an n-by-2 matrix where each row of the matrix specifies af starting time and an ending time. Entries can be numbers or strings representing expressions to be evaluated.

* model_params  a Dict which must contain key-value pairs with the keys :rule_and_delay_period, :target_period, and :post_target_period

# RETURNS

* a matrix the same size as opto_times, but with all strings parsed and evaluated.

# EXAMPLE

> model_params[:target_period] = 1  ; model_params[:rule_and_delay_period] = 1
> opto_times = ["exp(target_end)" 3; 6 8]
> parse_opto_times(opto_times, model_params)

7.38906   3
 6         8

"""
function parse_opto_times(opto_times2, model_params)
    # Make a copy so we don't futz with the original:
    opto_times = copy(opto_times2)
    # we want to be able to stash floats and numbers in by the end:
    if typeof(opto_times)<:Array{String}
        opto_times = convert(Array{Any}, opto_times)
    end

    
    # Let's define the time markers
    trial_start = target_start = target_end = trial_end = 0   # just define them here so vars are available ourside try/catch
    try
        trial_start  = 0
        target_start = model_params[:rule_and_delay_period]
        target_end   = target_start + model_params[:target_period]
        trial_end    = target_end   + model_params[:post_target_period]
    catch y
        @printf("\n\nwhoa, are you sure your model_params had all three of :rule_and_delay_period, :target_period, and :post_target_period?\n\n")
        error(y)
    end
    
    function replacer(P)   # run through an expression tree, replacing known symbols with their values, then evaluate
        if typeof(P)<:Symbol
            if P == :trial_start;  P = trial_start; end;                    
            if P == :target_start; P = target_start; end;                    
            if P == :target_end;   P = target_end; end;                    
            if P == :trial_end;    P = trial_end; end;                    
            return P
        end        
        for i=1:length(P.args)
            if typeof(P.args[i])<:Expr
                P.args[i] = replacer(P.args[i])
            else
                if P.args[i] == :trial_start;  P.args[i] = trial_start; end;                    
                if P.args[i] == :target_start || P.args[i]==100;  P.args[i] = target_start; end;                    
                if P.args[i] == :target_end   || P.args[i]==200;  P.args[i] = target_end; end;                    
                if P.args[i] == :trial_end    || P.args[i]==20;   P.args[i] = trial_end; end;                    
            end
        end
        return eval(P)
    end
    
    for i=1:length(opto_times)
        if typeof(opto_times[i])==String
            # @printf("Got a string, and it is: %s\n", opto_times[i])
            # @printf("Parsing turned it into: ");  print(parse(opto_times[i])); print("\n")
            # @printf("Replacer turned it into: "); print(replacer(parse(opto_times[i]))); print("\n")
            opto_times[i] = replacer(parse(opto_times[i]))
        elseif typeof(opto_times[i])==Int64
            if opto_times[i]==100; opto_times[i] = target_start; end;
            if opto_times[i]==200; opto_times[i] = target_end;   end;
            if opto_times[i]==20;  opto_times[i] = trial_end;    end;
        end            
    end
    return opto_times
end


# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


model_params = Dict(
:dt     =>  0.02,    # timestep, in secs
:tau    =>  0.1,     # tau, in ms
:vW     =>  -1.7,    # vertical weight
:hW     =>  -1.7,    # horizontal weight
:sW     =>  0.2,     # self-connection weight
:dW     =>  0,       # diagonal weight
:nsteps =>  2,       # number of timesteps in the simulation
:noise  =>  [],      # noise added during simulation. Can be empty matrix, or an nunits-by-nsteps matrix
:sigma  =>  0.08,    # standard deviation of Gaussian noise added (will be scaled by sqrt(dt) to be relatively dt-INsensitive)
:input  =>  0,       # input current. Can be scalar, nunits-by-1, or nunits-by-nsteps matrix
:g_leak =>  0.25,    # leak conductance
:U_rest =>  -1,      # resting membrane potential
:theta  =>  1,       # inverse slope of g() function
:beta   =>  1,       # offset to g() function
:constant_excitation      => 0.19,   # constant input, added to all units at all timesteps
:anti_rule_strength       => 0.1,    # input added only to anti units during rule_and_delay_period in Anti trials
:pro_rule_strength        => 0.1,    # input added only to pro units during rule_and_delay_period in Pro trials
:const_pro_bias           => 0,      # input added only to pro units during all times in all trial types
:target_period_excitation => 0.2,      # input added to all units during target_period
:right_light_excitation   => 0.3,    # input added to the Anti and the Pro unit on one side during the target_period
:right_light_pro_extra    => 0,      # input added to the right side Pro unit alone during the target_period
:rule_and_delay_period    => 0.4,    # duration of rule_and_delay_period, in secs
:target_period            => 0.1,    # duration of target_period, in secs
:post_target_period       => 0.5,    # duration of post_target_period, in secs
:const_add => 0,  # from rate_networks.jl, unused here
:init_add  => 0,  # from rate_networks.jl, unused here 
:opto_strength            => 1,      # fraction by which to scale V outputs
:opto_times               => [],     # n-by-2 matrix, indicating start and stop times for opto_strength effect, all within the same trial
:opto_units               => 1:4,    # ids of the units that will be affected by opto_strength effect
)


"""
input, t, nsteps = make_input(trial_type; dt=0.02, nderivs=0, difforder=0, constant_excitation=0.19, anti_rule_strength=0.1, 
    pro_rule_strength=0.1, target_period_excitation=1, right_light_excitation=0.5, right_light_pro_extra=0, 
    rule_and_delay_period=0.4, target_period=0.1, post_target_period=0.4, const_pro_bias=0,
    other_unused_params...)
"""
function make_input(trial_type; dt=0.02, nderivs=0, difforder=0, constant_excitation=0.19, anti_rule_strength=0.1, 
    pro_rule_strength=0.1, target_period_excitation=1, right_light_excitation=0.5, right_light_pro_extra=0, 
    rule_and_delay_period=0.4, target_period=0.1, post_target_period=0.4, const_pro_bias=0,
    other_unused_params...)

    T = rule_and_delay_period + target_period + post_target_period
    t = 0:dt:T
    nsteps = length(t)

    input = constant_excitation + ForwardDiffZeros(4, nsteps, nderivs=nderivs, difforder=difforder)
    if trial_type=="Anti"
        input[2:3, t.<rule_and_delay_period] += anti_rule_strength
    elseif trial_type=="Pro"
        input[[1,4], t.<rule_and_delay_period] += pro_rule_strength
    else
        error("make_input: I don't recognize input type \"" * trial_type * "\"")
    end
    
    input[:,     (rule_and_delay_period.<=t) & (t.<rule_and_delay_period+target_period)] += target_period_excitation
    input[1:2,   (rule_and_delay_period.<=t) & (t.<rule_and_delay_period+target_period)] += right_light_excitation
    input[1,     (rule_and_delay_period.<=t) & (t.<rule_and_delay_period+target_period)] += right_light_pro_extra
    
    input[[1,4],:] += const_pro_bias
    
    return input, t, nsteps
end



"""
proVs, antiVs, pro_fullV, anti_fullV = run_ntrials(nPro, nAnti; plot_list=[], nderivs=0, difforder=0, 
    model_params...)

Runs a set of proAnti model trials.  See model_params above for definition of all the parameters. In addition,

# PARAMETERS

- nPro    number of Pro tials to run

- nAnti   number of Anti trials to run

# OPTIONAL PARAMETERS

- plot_list    A list of trials to plot on the figures. If empty nothing is plotted

- nderivs      For ForwardDiff

- difforder    For ForwardDiff

- model_params   Further optional params, will be passed onto forwardModel() (e.g., opto times)

# RETURNS

- proVs    final V for the four units on Pro trials

- antiVs   final V for the four units on Anti trials

- pro_fullVs   all Vs, across all times, for all Pro trials

- anti_fullVs  all Vs, across all times, for all Anti trials

"""
function run_ntrials(nPro, nAnti; plot_list=[], nderivs=0, difforder=0, model_params...)
    
    pro_input,  t, nsteps = make_input("Pro" ; nderivs=nderivs, difforder=difforder, model_params...)
    anti_input, t, nsteps = make_input("Anti"; nderivs=nderivs, difforder=difforder, model_params...)
    
    model_params = Dict(model_params)
    sW = model_params[:sW]
    hW = model_params[:hW]
    vW = model_params[:vW]
    dW = model_params[:dW]
    model_params = make_dict(["nsteps", "W"], [nsteps, [sW vW dW hW; vW sW hW dW; dW hW sW vW; hW dW vW sW]], 
        model_params)
    model_params = make_dict(["nderivs", "difforder"], [nderivs, difforder], model_params)
    model_params[:opto_times] = parse_opto_times(model_params[:opto_times], model_params)
    
    proVs  = ForwardDiffZeros(4, nPro, nderivs=nderivs, difforder=difforder)
    antiVs = ForwardDiffZeros(4, nAnti, nderivs=nderivs, difforder=difforder)
    pro_fullV  = Array{Any}(nPro,1)
    anti_fullV = Array{Any}(nAnti,1)
    
    # --- PRO ---
    if length(plot_list)>0; figure(1); clf(); end
    model_params = make_dict(["input"], [pro_input], model_params)
    for i=1:nPro
        startU = [-0.3, -0.7, -0.7, -0.3]
        Uend, Vend, pro_fullU, pro_fullV[i] = forwardModel(startU, do_plot=false; model_params...)
        proVs[:,i] = Vend
        if any(plot_list.==i) 
            plot_PA(t, pro_fullU, pro_fullV[i]; fignum=1, clearfig=false, model_params...)
            subplot(3,1,1); title("PRO")
        end
    end

    # --- ANTI ---
    if length(plot_list)>0; figure(2); clf(); end
    model_params = make_dict(["input"], [anti_input], model_params)
    for i=1:nAnti
        startU = [-0.7, -0.3, -0.3, -0.7]
        Uend, Vend, anti_fullU, anti_fullV[i] = forwardModel(startU, do_plot=false; model_params...)
        antiVs[:,i] = Vend
        if any(plot_list.==i) 
            plot_PA(t, anti_fullU, anti_fullV[i]; fignum=2, clearfig=false, model_params...)
            subplot(3,1,1); title("ANTI")
        end
    end
    
    return proVs, antiVs, pro_fullV, anti_fullV
end


# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


function JJ(nPro, nAnti; pro_target=0.9, anti_target=0.7, 
    theta1=0.025, theta2=0.035, cbeta=0.003, verbose=false, 
    pre_string="", zero_last_sigmas=0, seedrand=NaN, 
    rule_and_delay_periods = [0.4], target_periods = [0.1], post_target_periods = [0.5],
    nderivs=0, difforder=0, model_params...)

    nruns = length(rule_and_delay_periods)*length(target_periods)*length(post_target_periods)
    
    cost1s = ForwardDiffZeros(1, nruns, nderivs=nderivs, difforder=difforder)
    cost2s = ForwardDiffZeros(1, nruns, nderivs=nderivs, difforder=difforder)

    if ~isnan(seedrand); srand(seedrand); end
    
    n = totHitsP = totHitsA = totDiffsP = totDiffsA = 0
    for i in rule_and_delay_periods
        for j in target_periods
            for k = post_target_periods
                n += 1
                
                my_params = make_dict(["rule_and_delay_period", "target_period", "post_target_period"],
                [i, j, k], Dict(model_params))
    
                # print("model params is " ); print(model_params); print("\n")
                proVs, antiVs = run_ntrials(nPro, nAnti; nderivs=nderivs, difforder=difforder, my_params...)

                hitsP  = 0.5*(1 + tanh.((proVs[1,:]-proVs[4,:,])/theta1))
                diffsP = tanh.((proVs[1,:,]-proVs[4,:])/theta2).^2
                hitsA  = 0.5*(1 + tanh.((antiVs[4,:]-antiVs[1,:,])/theta1))
                diffsA = tanh.((antiVs[4,:,]-antiVs[1,:])/theta2).^2

                if nPro>0 && nAnti>0
                    cost1s[n] = (nPro*(mean(hitsP) - pro_target).^2  + nAnti*(mean(hitsA) - anti_target).^2)/(nPro+nAnti)
                    cost2s[n] = -cbeta*(nPro*mean(diffsP) + nAnti*mean(diffsA))/(nPro+nAnti)
                elseif nPro>0
                    cost1s[n] = (mean(hitsP) - pro_target).^2
                    cost2s[n] = -cbeta*mean(diffsP)
                else
                    cost1s[n] = (mean(hitsA) - anti_target).^2
                    cost2s[n] = -cbeta*mean(diffsA)
                end

                totHitsP  += mean(hitsP);  totHitsA  += mean(hitsA); 
                totDiffsP += mean(diffsP); totDiffsA += mean(diffsA);
            end
        end
    end
    
    cost1 = mean(cost1s)
    cost2 = mean(cost2s)

    hitsP = totHitsP/n; hitsA = totHitsA/n; diffsP = totDiffsP/n; diffsA = totDiffsA/n
    
    if verbose
        @printf("%s", pre_string)
        @printf("     -- cost=%g,   cost1=%g, cost2=%g\n", 
            convert(Float64, cost1+cost2), convert(Float64, cost1), convert(Float64, cost2))
        if nPro>0 && nAnti>0
            @printf("     -- mean(hitsP)=%g, mean(diffsP)=%g mean(hitsA)=%g, mean(diffsA)=%g\n", 
                convert(Float64, mean(hitsP)), convert(Float64, mean(diffsP)),
                convert(Float64, mean(hitsA)), convert(Float64, mean(diffsA)))
        elseif nPro>0
            @printf("     -- mean(hitsP)=%g, mean(diffsP)=%g (nAnti=0)\n", 
                convert(Float64, mean(hitsP)), convert(Float64, mean(diffsP)))
        else
            @printf("     -- (nPro=0) mean(hitsA)=%g, mean(diffsA)=%g\n", 
                convert(Float64, mean(hitsA)), convert(Float64, mean(diffsA)))
        end        
    end
    
    return cost1 + cost2
end


