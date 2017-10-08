# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


include("rate_networks.jl")  # that will also include general_utils.jl, constrained_parabolic_minimization.jl, and hessian_utils.jl



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
    plot_Us=true, ax_set=Dict(), other_unused_params...)
    figure(fignum)
    if clearfig && ~isempty(ax_set); clf(); end
    
    if ~haskey(ax_set, "Vax")
        if plot_Us; ax1 = subplot(3,1,1); else ax1=subplot(2,1,1); end;
    else
        ax1 = axes(ax_set["Vax"], fontsize=20)
    end
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
        
    if plot_Us
        if ~haskey(ax_set, "Uax")
            ax2 = subplot(3,1,2)
        else
            ax2 = axes(ax_set["Uax"])
        end
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
    end
    
    if ~haskey(ax_set, "Dax")
        if plot_Us; ax3 = subplot(3,1,3); else ax3=subplot(2,1,2); end;
    else
        ax3 = axes(ax_set["Dax"])
    end

    delta = V[1,:] - V[4,:]
    hr = plot(t, delta, "g")
    oldlims = [ylim()[1]+0.1, ylim()[2]-0.1]
    ylim(minimum([delta[:];oldlims[1]])-0.1, maximum([delta[:];oldlims[2]])+0.1)
    vlines([rule_and_delay_period, 
            rule_and_delay_period+target_period,
            rule_and_delay_period+target_period+post_target_period], 
            ylim()[1], ylim()[2], linewidth=2)
    hlines(0, t[1], t[end], linewidth=2)
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
> opto_times = ["exp(target_end)" 3; 6 8]   # NOTE: exp() that means the exponential-- any arbitrary Julia expression is allowed
> parse_opto_times(opto_times, model_params)

7.38906   3
 6         8

"""
function parse_opto_times(opto_times2, model_params)
    # Make a copy so we don't futz with the original:
    opto_times = copy(opto_times2)
    
    # we want to be able to stash floats and numbers in by the end:
    if typeof(opto_times)<:Array{String}  ||  typeof(opto_times)<:Array{Int64}
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

    if FDversion() >= 0.6
        # All the variables that we MIGHT choose to differentiate w.r.t. go into this bag -- further down
        # we'll use get_eltype(varbag) to check for any of them being ForwardDiff.Dual.
        # That is how we'll tell whether new matrices should be regular numbers of ForwardDiff.Dual's.
        # *** if you add a new variable you'll want to differentiate w.r.t., it should be added here too ***
        varbag = (dt, constant_excitation, anti_rule_strength, pro_rule_strength, target_period_excitation, 
        right_light_excitation, right_light_pro_extra, rule_and_delay_period, target_period, post_target_period,
        const_pro_bias)
    end
    
    T = rule_and_delay_period + target_period + post_target_period
    t = 0:dt:T
    nsteps = length(t)
    
    if FDversion() < 0.6
        input = constant_excitation + ForwardDiffZeros(4, nsteps, nderivs=nderivs, difforder=difforder)
    else
        input = constant_excitation + zeros(get_eltype(varbag), 4, nsteps)
    end
    
    if trial_type=="Anti"
        input[2:3, t.<rule_and_delay_period] += anti_rule_strength
    elseif trial_type=="Pro"
        input[[1,4], t.<rule_and_delay_period] += pro_rule_strength
    else
        error("make_input: I don't recognize input type \"" * trial_type * "\"")
    end
    
    if VERSION.major + 0.1*VERSION.minor < 0.6
        input[:,     (rule_and_delay_period.<=t) & (t.<rule_and_delay_period+target_period)] += target_period_excitation
        input[1:2,   (rule_and_delay_period.<=t) & (t.<rule_and_delay_period+target_period)] += right_light_excitation
        input[1,     (rule_and_delay_period.<=t) & (t.<rule_and_delay_period+target_period)] += right_light_pro_extra
    else    
        input[:,     (rule_and_delay_period.<=t) .& (t.<rule_and_delay_period+target_period)] += target_period_excitation
        input[1:2,   (rule_and_delay_period.<=t) .& (t.<rule_and_delay_period+target_period)] += right_light_excitation
        input[1,     (rule_and_delay_period.<=t) .& (t.<rule_and_delay_period+target_period)] += right_light_pro_extra
    end    
    input[[1,4],:] += const_pro_bias
    
    return input, t, nsteps
end



"""
proVs, antiVs, pro_fullV, anti_fullV, opto_fraction, pro_input, anti_input = 
    run_ntrials(nPro, nAnti; plot_list=[], start_pro=[-0.5,-0.5,-0.5,-0.5], start_anti=[-0.5,-0.5,-0.5,-0.5],
        profig=1, antifig=2,
        opto_units = 1:4, nderivs=0, difforder=0, model_params...)

Runs a set of proAnti model trials.  See model_params above for definition of all the parameters. In addition,

# PARAMETERS

- nPro    number of Pro tials to run

- nAnti   number of Anti trials to run

# OPTIONAL PARAMETERS

- plot_list    A list of trials to plot on the figures. If empty nothing is plotted

- profig       The figure number on which Pro traces will be drawn

- antifig      The figure number on which Anti traces will be drawn

- start_pro    A 4-by-1 vector of starting U values on Pro trials for the 4 units

- start_anti   A 4-by-1 vector of starting U values on Anti trials for the 4 units

- nderivs      For ForwardDiff

- difforder    For ForwardDiff

- model_params   Further optional params, will be passed onto forwardModel() (e.g., opto times)


# RETURNS

- proVs    final V for the four units on Pro trials

- antiVs   final V for the four units on Anti trials

- pro_fullVs   all Vs, across all times, for all Pro trials (4-by-nsteps-by-nPro)

- anti_fullVs  all Vs, across all times, for all Anti trials (4-by-nsteps-by-nAnti)

"""
function run_ntrials(nPro, nAnti; plot_list=[], start_pro=[-0.5,-0.5,-0.5,-0.5], start_anti=[-0.5,-0.5,-0.5,-0.5],
    profig=1, antifig=2, plot_Us=true,  clearfig=true, ax_set = Dict(),
    opto_units = 1:4, nderivs=0, difforder=0, model_params...)

    if FDversion() >= 0.6
        # All the variables that we MIGHT choose to differentiate w.r.t. go into this bag -- further down
        # we'll use get_eltype(varbag) to check for any of them being ForwardDiff.Dual.
        # That is how we'll tell whether new matrices should be regular numbers of ForwardDiff.Dual's.
        # *** if you add a new variable you'll want to differentiate w.r.t., it should be added here too ***
        varbag = (start_pro, start_anti, model_params)
        # print("get_eltype(varbag)="); print(get_eltype(varbag)); print("\n")
    end
    
    pro_ax_set = Dict()
    anti_ax_set = Dict()
    if haskey(ax_set, "pro_Vax");  pro_ax_set["Vax"]  = ax_set["pro_Vax"]; end;
    if haskey(ax_set, "pro_Uax");  pro_ax_set["Uax"]  = ax_set["pro_Uax"]; end;
    if haskey(ax_set, "pro_Dax");  pro_ax_set["Dax"]  = ax_set["pro_Dax"]; end;                        
    if haskey(ax_set, "anti_Vax"); anti_ax_set["Vax"] = ax_set["anti_Vax"]; end;
    if haskey(ax_set, "anti_Dax"); anti_ax_set["Dax"] = ax_set["anti_Dax"]; end;
    if haskey(ax_set, "anti_Uax"); anti_ax_set["Uax"] = ax_set["anti_Uax"]; end;
    
    model_params = Dict(model_params)
    pro_input,  t, nsteps = make_input("Pro" ; nderivs=nderivs, difforder=difforder, model_params...)
    anti_input, t, nsteps = make_input("Anti"; nderivs=nderivs, difforder=difforder, model_params...)
    
    sW = model_params[:sW]
    hW = model_params[:hW]
    vW = model_params[:vW]
    dW = model_params[:dW]
    model_params = make_dict(["nsteps", "W"], [nsteps, [sW vW dW hW; vW sW hW dW; dW hW sW vW; hW dW vW sW]], 
        model_params)
    model_params = make_dict(["nderivs", "difforder"], [nderivs, difforder], model_params)
    model_params[:opto_times] = parse_opto_times(model_params[:opto_times], model_params)
    
    if FDversion()<0.6
        proVs  = ForwardDiffZeros(4, nPro, nderivs=nderivs, difforder=difforder)
        antiVs = ForwardDiffZeros(4, nAnti, nderivs=nderivs, difforder=difforder)
    else
        proVs  = zeros(get_eltype(varbag), 4, nPro)
        antiVs = zeros(get_eltype(varbag), 4, nAnti)
    end
    proVall  = zeros(4, nsteps, nPro);
    antiVall = zeros(4, nsteps, nAnti);
    # --- PRO ---
    if length(plot_list)>0; figure(profig); if clearfig; clf(); end; end
    model_params = make_dict(["input"], [pro_input], model_params)
    for i=1:nPro
        Uend, Vend, pro_fullU, temp = forwardModel(start_pro, do_plot=false, opto_units=opto_units; model_params...)
        if FDversion() < 0.6; proVall[:,:,i] = temp; else proVall[:,:,i] = get_value(temp); end
        # print("typeof(proVs) = "); print(typeof(proVs)); print("\n")
        proVs[:,i] = Vend
        if any(plot_list.==i) 
            plot_PA(t, get_value(pro_fullU), get_value(proVall[:,:,i]); clearfig=false,
                fignum=profig, ax_set=pro_ax_set, plot_Us=plot_Us, clearfig=false, model_params...)
            # subplot(3,1,1); title("PRO")
        end
    end

    # --- ANTI ---
    if length(plot_list)>0; figure(antifig); if clearfig; clf(); end; end
    model_params = make_dict(["input"], [anti_input], model_params)
    for i=1:nAnti
        startU = [-0.5, -0.5, -0.5, -0.5]
        Uend, Vend, anti_fullU, temp = forwardModel(start_anti, do_plot=false, opto_units=opto_units; model_params...)
        if FDversion() < 0.6; antiVall[:,:,i] = temp; else antiVall[:,:,i] = get_value(temp); end
        antiVs[:,i] = Vend
        if any(plot_list.==i) 
            plot_PA(t, get_value(anti_fullU), get_value(antiVall[:,:,i]); clearfig=false,
            fignum=antifig, ax_set=anti_ax_set, plot_Us=plot_Us, clearfig=false, model_params...)
            # subplot(3,1,1); title("ANTI")
        end
    end
        
    if haskey(model_params, :opto_fraction)
        opto_fraction = model_params[:opto_fraction]
    else
        opto_fraction = 1
    end
    
    return proVs, antiVs, proVall, antiVall, opto_fraction, pro_input, anti_input
end


# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


"""
cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA, [proValls, antiValls, opto_fraction, pro_input, anti_input] = 
    JJ(nPro, nAnti; pro_target=0.9, anti_target=0.7, 
    opto_targets = Array{Float64}(2,0), opto_periods = Array{Float64}(2,0), 
    model_details = false,
    theta1=0.025, theta2=0.035, cbeta=0.003, verbose=false, 
    pre_string="", zero_last_sigmas=0, seedrand=NaN, 
    rule_and_delay_periods = nothing, target_periods = [0.1], post_target_periods = nothing,
    plot_conditions=false, plot_list = [],
    nderivs=0, difforder=0, model_params...)

Runs a proAnti network, if desired across multiple opto conditions and multiple period durations and returns
resulting cost.

If rule_and_delay_periods and post_target_periods are not passed, it tries to get them from their
singular (not plural) versions in model_params, e.g., model_params[:rule_and_delay_period]. NOTE that
target_period is different, for backwards compatibility it defaults to 0.1.

# PARAMETERS (INCOMPLETE DOCS!!!):

??

# RETURNS (INCOMPLETE DOCS!!!):

- cost   The net cost, composed of squared error cost (promoting performance close to the desired one) plus difference cost (promoting Pro_R and Pro_L differences at the end of the trial)

- cost1s  A vector, for each opto condition, the squared error cost

- cost2s  for each opto condition, the differences costc

- hP     Pro "hits", as computed with the theta1 sigmoid
    
- hA     Anti "hits"

- dP     Pro "diffs", as computed with the theta2 sigmoid
    
- dA     Anti "diffs"

- hBP    Pro binarized hits, as computed by binarizing (equivalent to theta1->0)
    
- hBA    Anti binarized hits

If model_details is set to true, also returns proValls, antiValls, opto_fraction, pro_input, anti_input



"""

function JJ(nPro, nAnti; pro_target=0.9, anti_target=0.7, 
    opto_targets = Array{Float64}(2,0), opto_periods = Array{Float64}(2,0), 
    model_details = false,
    theta1=0.025, theta2=0.035, cbeta=0.003, verbose=false, 
    pre_string="", zero_last_sigmas=0, seedrand=NaN, 
    rule_and_delay_periods = nothing, target_periods = nothing, post_target_periods = nothing,
    plot_conditions=false, plot_list = [],
    nderivs=0, difforder=0, model_params...)

    if FDversion() >= 0.6
        # All the variables that we MIGHT choose to differentiate w.r.t. go into this bag -- further down
        # we'll use get_eltype(varbag) to check for any of them being ForwardDiff.Dual.
        # That is how we'll tell whether new matrices should be regular numbers of ForwardDiff.Dual's.
        # *** if you add a new variable you'll want to differentiate w.r.t., it should be added here too ***
        varbag = (pro_target, anti_target, opto_targets, opto_periods, theta1, theta2, cbeta,model_params)
        # print("get_eltype(varbag)="); print(get_eltype(varbag)); print("\n")
    end
    
    # If the plurals of the periods are not passed in, then use the singular in model_params as the default:
    if rule_and_delay_periods==nothing
        rule_and_delay_periods = model_params[:rule_and_delay_period]
    end
    if target_periods==nothing        
        target_periods = [0.1]
        @printf("\n\nWARNING: JJ() was not given a target_periods parameter, I will *IGNORE* the\n")
        @printf("    model_params[:target_period] entry and will use a default of 0.1\n\n")
    end
    if post_target_periods==nothing
         post_target_periods = model_params[:post_target_period]
    end
    
    if size(opto_targets,1)==0  || size(opto_periods,1)==0 # if there's no opto that is being asked for
        # Then tun with only a single opto_period request, with no opto, and control targets as our targets
        opto_targets = [pro_target anti_target]
        opto_periods = [-1, -1]   # opto effect is before time zero, i.e., is nothing
    end
    
    if ~(size(opto_targets) == size(opto_periods)); 
        error("opto parameters are bad -- need a Pro and Anti performance target for each requested period"); 
    end

    noptos     = size(opto_periods,1)  # of opto conditions
    nruns_each = length(rule_and_delay_periods)*length(target_periods)*length(post_target_periods)    # runs per opto condition
    nruns      = nruns_each*noptos  # total conditions

    if FDversion() < 0.6
        cost1s = ForwardDiffZeros(noptos, nruns, nderivs=nderivs, difforder=difforder)
        cost2s = ForwardDiffZeros(noptos, nruns, nderivs=nderivs, difforder=difforder)
    else
        cost1s = zeros(get_eltype(varbag), noptos, nruns)
        cost2s = zeros(get_eltype(varbag), noptos, nruns)
    end
            
    hP  = zeros(noptos, nruns_each);   # Pro "hits", as computed with the theta1 sigmoid
    hA  = zeros(noptos, nruns_each);   # Anti "hits"
    dP  = zeros(noptos, nruns_each);   # Pro "diffs", as computed with the theta2 sigmoid
    dA  = zeros(noptos, nruns_each);   # Anti "diffs"
    hBP = zeros(noptos, nruns_each);   # Pro binarized hits, as computed by binarizing (equivalent to theta1->0)
    hBA = zeros(noptos, nruns_each);   # Anti binarized hits

    proValls         = [];
    antiValls        = [];
    opto_fraction    = [];
    pro_input        = [];
    anti_input       = [];
    
    for nopto=1:noptos # iterate over each opto inactivation period
        # @printf("size(hBP) is %d, %d\n", size(hBP,1), size(hBP,2))

        # reset random number generator for each opto period, so it cant over fit noise samples
        if ~isnan(seedrand); srand(seedrand); end

        n = 0  # n is a counter over all period duration conditions
        totHitsP = totHitsA = totDiffsP = totDiffsA = 0
        for i in rule_and_delay_periods
            for j in target_periods
                for k = post_target_periods
                    n += 1

                    # include this opto inactivation in the parameters to pass on
                    my_params = make_dict(["rule_and_delay_period","target_period","post_target_period"], [i,j,k])
                    my_params = make_dict(["opto_times"], [reshape(opto_periods[nopto,:], 1, 2)], my_params)
                    my_params = merge(Dict(model_params), my_params)  # my_params takes precedence

                    my_plot_list = [];
                    if typeof(plot_conditions)==Bool && plot_conditions
                        my_plot_list = plot_list;
                    elseif typeof(plot_conditions)<:Array && plot_conditions[nopto]
                        my_plot_list = plot_list;
                    end

                    # print("model params is " ); print(model_params); print("\n")
                    proVs, antiVs, proVall, antiVall, opto_fraction,pro_input,anti_input =
                        run_ntrials(nPro, nAnti; plot_list=my_plot_list, 
                            nderivs=nderivs, difforder=difforder, my_params...)
                        # run_ntrials_opto(nPro, nAnti; nderivs=nderivs, difforder=difforder, my_params...)

                    if length(proValls)==0
                        proValls = zeros(4, size(proVall,2), size(proVall,3), noptos)
                    end
                    if length(antiValls)==0
                        antiValls = zeros(4, size(antiVall,2), size(antiVall,3), noptos)
                    end
                    proValls[:,:,:,nopto]  = get_value(proVall)
                    antiValls[:,:,:,nopto] = get_value(antiVall)
                    # @printf("size of proValls is "); print(size(proValls)); print("\n")
                    
                    hitsP  = 0.5*(1 + tanh.((proVs[1,:]-proVs[4,:,])/theta1))
                    diffsP = tanh.((proVs[1,:,]-proVs[4,:])/theta2).^2
                    hitsA  = 0.5*(1 + tanh.((antiVs[4,:]-antiVs[1,:,])/theta1))
                    diffsA = tanh.((antiVs[4,:,]-antiVs[1,:])/theta2).^2

                    # set up storage  -- we do get_value() to make sure to from ForwardDiff.Dual into Float64 if necessary
                    hP[nopto, n] = mean(get_value(hitsP));
                    hA[nopto, n] = mean(get_value(hitsA));
                    dP[nopto, n] = mean(get_value(diffsP));
                    dA[nopto, n] = mean(get_value(diffsA));
                    hBP[nopto, n] = get_value(sum(proVs[1,:] .>= proVs[4,:,])/nPro);
                    hBA[nopto, n] = get_value(sum(antiVs[4,:] .>  antiVs[1,:,])/nAnti);                    
                    
                    if nPro>0 && nAnti>0
                        # cost1s and cost2s can accept ForwardDiff.Dual, so no get_value() for them
                        cost1s[nopto, n] = (nPro*(mean(hitsP) - opto_targets[nopto,1]).^2  
                                    + nAnti*(mean(hitsA) - opto_targets[nopto,2]).^2)/(nPro+nAnti)
                        cost2s[nopto, n] = -cbeta*(nPro*mean(diffsP) + nAnti*mean(diffsA))/(nPro+nAnti)
                    elseif nPro>0
                        cost1s[nopto, n] = (mean(hitsP) - opto_targets[nopto,1]).^2
                        cost2s[nopto, n] = -cbeta*mean(diffsP)
                    else
                        cost1s[nopto, n] = (mean(hitsA) - opto_targets[nopto,2]).^2
                        cost2s[nopto, n] = -cbeta*mean(diffsA)
                    end
                    totHitsP  += mean(hitsP);  totHitsA  += mean(hitsA); 
                    totDiffsP += mean(diffsP); totDiffsA += mean(diffsA);
                end
            end
        end
        ## @printf("size(hBP) is %d, %d\n", size(hBP,1), size(hBP,2))
        hitsP = totHitsP/n; hitsA = totHitsA/n; diffsP = totDiffsP/n; diffsA = totDiffsA/n
    
        if verbose
            pcost1 = mean(cost1s[nopto,:])   # partial costs
            pcost2 = mean(cost2s[nopto,:])            
            
            # Notice the get_value() calls below, to transform ForwardDiff Duals into Float64s
            @printf("%s", pre_string)
            @printf("Opto condition # %d\n", nopto)
            @printf("     - %d - cost=%g, cost1=%g, cost2=%g\n", nopto,
                get_value(pcost1+pcost2), get_value(pcost1), get_value(pcost2))
            if nPro>0 && nAnti>0
                @printf("     - %d - mean(hitsP)=%g, mean(diffsP)=%g mean(hitsA)=%g, mean(diffsA)=%g\n", nopto,
                    get_value(mean(hitsP)), get_value(mean(diffsP)),
                    get_value(mean(hitsA)), get_value(mean(diffsA)))
            elseif nPro>0
                @printf("     - %d - mean(hitsP)=%g, mean(diffsP)=%g (nAnti=0)\n", nopto,
                    get_value(mean(hitsP)), get_value(mean(diffsP)))
            else
                @printf("     - %d - (nPro=0) mean(hitsA)=%g, mean(diffsA)=%g\n", nopto,
                    get_value(mean(hitsA)), get_value(mean(diffsA)))
            end        
        end
    end
        
    # @printf("size(hBP) is %d, %d\n", size(hBP,1), size(hBP,2))

    cost1 = mean(cost1s)
    cost2 = mean(cost2s)

    if verbose
        @printf("%s", pre_string)
        @printf("OVERALL\n")
        @printf("     -- cost=%g, cost1=%g, cost2=%g\n", 
            get_value(cost1+cost2), get_value(cost1), get_value(cost2))
    end
    
    if model_details
        return cost1 + cost2, cost1s, cost2s, hP,hA,dP,dA,hBP,hBA, 
            proValls, antiValls, opto_fraction, pro_input, anti_input
    else
        return cost1 + cost2, cost1s, cost2s, hP,hA,dP,dA,hBP,hBA
    end
end                    


# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb



"""
model_params, F, nPro, nAnti = load_run(run_name; farmdir="FarmFields")

Loads a run and sets everything into a self-contained model_params so that you could then run it directly:

    JJ(model_params[:nPro], model_params[:nAnti]; model_params...);

If the model_params for the run included :start_pro and :start_anti entries (or "start_pro" and
"start_anti" entries in F), it uses those. Otherwise
sets :start_pro and :start_anti to a default of [-0.5, -0.5, -0.5, -0.5]


# PARAMETERS:

- run_name   A string representing the run. If it doesn't end in .mat, the .mat is added to it

# OPTIONAL PARAMETERS:

- farmdir   A string representing the directory in which the run is found

# RETURNS:

- model_params    The dictionary with all necessary params

- F               A dictionary with the raw matread of the run file

- nPro            equals model_params[:nPro]

- nAnti           equals model_params[:nAnti]



"""
function load_run(run_name; farmdir="FarmFields")

    default_U_start = [-0.5, -0.5, -0.5, -0.5]
    
    if !endswith(run_name, ".mat")
        run_name = run_name * ".mat"
    end
    
    F = matread(farmdir * "/" * run_name)
    model_params = symbol_key_ize(F["model_params"])
    model_params[:rule_and_delay_periods] = F["rule_and_delay_periods"]
    model_params[:rule_and_delay_period]  = model_params[:rule_and_delay_periods][1]
    model_params[:target_period]          = model_params[:target_period]
    model_params[:post_target_periods]    = F["post_target_periods"]
    model_params[:post_target_period]     = model_params[:post_target_periods][1]
    model_params[:seedrand]=F["test_sr"]
    model_params[:cbeta] =F["cb"]
    if ~haskey(model_params, :start_pro)
        if ~haskey(F, "start_pro")
            model_params[:start_pro] = default_U_start
        else
            model_params[:start_pro] = F["start_pro"]
        end
    end
    if ~haskey(model_params, :start_anti)
        if ~haskey(F, "start_anti")
            model_params[:start_anti] = default_U_start
        else
            model_params[:start_anti] = F["start_antia"]
        end
    end
    model_params = make_dict(F["args"], F["pars"], model_params)

    nPro = model_params[:nPro]
    nAnti = model_params[:nAnti]
    
    if ~haskey(model_params, :target_periods)
        model_params[:target_periods] = 0.1  # Use the JJ() default
        model_params[:target_period]  = 0.1  # Use the JJ() default
        @printf("\n\nWARNING: :target_periods was not specified, JJ() would ignore :target_period and\n")
        @printf("would use a default of :target_periods=>[0.1].  Setting that explicitly here, no warning from JJ()\n")
        @printf("will be elicited.\n\n")
    end

    
    return model_params, F, nPro, nAnti
end


