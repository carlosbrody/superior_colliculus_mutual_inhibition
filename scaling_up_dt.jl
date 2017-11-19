# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


# In its own cell so we can run it just once

include("pro_anti.jl")   # Loads all the necessary pre-requisites



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


model_params, F, nPro, nAnti = load_run("farm_LD0003"; farmdir="goodfarms")

# mypars = merge(model_params, Dict(:post_target_periods=>0.05))

cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA, proValls, antiValls, opto_fraction, pro_input, anti_input = 
JJ(model_params[:nPro], model_params[:nAnti]; verbose=true, model_details=true, model_params...);


epochs = ["control";"full";"rule";"delay";"choice"]
titles  = ["epochs" "hitP" "hitA" "diffP" "diffA" "hBP" "hBA"]

@printf("\n\n--- The original run:\n\n")

sleep(0.5)  # A pause just to make sure everything has printed out 
display([titles; epochs hP hA dP dA hBP hBA])



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


# Now scale things to make dt and tau are much larger, and involve 0.25 x the timesteps

mypars = merge(model_params, Dict(:opto_times=>[0 0]))  # default is no opto

mypars = merge(model_params, Dict(
:dt=>0.024,
:tau=>0.09,
:rule_and_delay_period=>1.2,
:target_period=>0.3,
:post_target_period=>0.3,
:anti_rule_strength=>0.054,
# :post_target_period=>0.05, 
# :rule_and_delay_period=>0.8, 
# :opto_times=>[0 0],
# :opto_times=>["target_start/2" "target_start"],
# :opto_times=>["target_start", "target_end"],
:opto_strength=>0.85,
))# , :rule_and_delay_period=>0.8))


mypars[:rule_and_delay_periods] = [mypars[:rule_and_delay_period]]
mypars[:target_periods]         = [mypars[:target_period]]
mypars[:post_target_periods]    = [mypars[:post_target_period]]

mypars[:opto_periods] = [
    "trial_start"     "trial_start"  ;
    "trial_start"     "trial_end"    ; 
    "trial_start"     "target_start/2"  ;
    "target_start/2"  "target_start"  ;
    "target_start"    "target_end"
];



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


# And run with these params.  

# We run 1000 trials of Pro and ANti each, just to get a good estimate of probability correct in each condition

mypars[:seedrand] = Int64(round(time()*1000))

ntrials = 1000

cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA, proValls, antiValls, opto_fraction, pro_input, anti_input = 
JJ(ntrials, ntrials; verbose=true, model_details=true, mypars...);

@printf("\n\n--- The hand-modified run at dt=0.024 at 4x fewer timesteps than the original:\n\n")


sleep(0.5)  # A pause just to make sure everything has printed out 
display([titles; epochs hP hA dP dA hBP hBA])



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


# Finally, display some trajectories:

pygui(true)

# This example is with opto during the delay

mypars = merge(mypars, Dict(
# :opto_times => ["target_start/2"    "target_start"],
:opto_times => [0    0],
:sigma=>0.01
))


proVs, antiVs, pro_fullV, anti_fullV, opto_fraction, pro_input, anti_input = run_ntrials(1000, 1000; 
    plot_list=[1:40;], profig=1, antifig=2, opto_units = 1:4, mypars...);

Np = size(proVs,2); Na = size(antiVs,2)
@printf("hBP = %g, hBA = %g\n", length(find(proVs[1,:] .>= proVs[4,:]))/Np, length(find(antiVs[4,:] .> antiVs[1,:]))/Na)

@printf("\n\n\n---\n\nmypars contains all the parameters for the new run with the scaled up dt, tau, etc.\n\n")



