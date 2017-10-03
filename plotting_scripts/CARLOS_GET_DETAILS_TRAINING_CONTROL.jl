using MAT
using ForwardDiff
using DiffBase

include("pro_anti.jl")
include("pro_anti_opto.jl")

function JJ_opto_plot(nPro, nAnti; opto_targets=[0.9 0.7], theta1=0.025, theta2=0.035, cbeta=0.003, verbose=false, pre_string="", zero_last_sigmas=0, seedrand=NaN, rule_and_delay_periods = [0.4], target_periods = [0.1], post_target_periods = [0.5], opto_periods = [-1 -1],opto_strength=1, nderivs=0, difforder=0,plot_conditions=false,model_details=false, model_params...) #set opto defaults!

    if ~(size(opto_targets) == size(opto_periods)); error("opto parameters are bad"); end

    nruns = length(rule_and_delay_periods)*length(target_periods)*length(post_target_periods)*size(opto_periods)[1]
    nruns_each = length(rule_and_delay_periods)*length(target_periods)*length(post_target_periods)
    
    cost1s = ForwardDiffZeros(size(opto_periods)[1], nruns_each, nderivs=nderivs, difforder=difforder)
    cost2s = ForwardDiffZeros(size(opto_periods)[1], nruns_each, nderivs=nderivs, difforder=difforder)
    hP = zeros(size(opto_periods)[1], nruns_each);
    hA = zeros(size(opto_periods)[1], nruns_each);
    dP = zeros(size(opto_periods)[1], nruns_each);
    dA = zeros(size(opto_periods)[1], nruns_each);
    hBP = zeros(size(opto_periods)[1], nruns_each);
    hBA = zeros(size(opto_periods)[1], nruns_each);

    if model_details
        proVall         = [];
        antiVall        = [];
        opto_fraction   = [];
        pro_input       = [];
        anti_input      = [];
    end

    n = totHitsP = totHitsA = totDiffsP = totDiffsA =nopto= 0
    for kk=1:size(opto_periods)[1] # iterate over each opto inactivation period
    nopto = 0;

    # reset random number generator for each opto period, so it cant over fit noise samples
    if ~isnan(seedrand); srand(seedrand); end

    for i in rule_and_delay_periods
        for j in target_periods
            for k = post_target_periods
                nopto += 1
                
                # include this opto inactivation in the parameters to pass on
                my_params = make_dict(["rule_and_delay_period", "target_period", "post_target_period","opto_period","opto_strength"], [i, j, k, opto_periods[kk,:], opto_strength], Dict(model_params))
    
                # print("model params is " ); print(model_params); print("\n")
                if typeof(plot_conditions)==Bool && ~plot_conditions
                    proVs, antiVs, proVall, antiVall, opto_fraction,pro_input,anti_input  = run_ntrials_opto(nPro, nAnti; nderivs=nderivs, difforder=difforder, my_params...)
                elseif typeof(plot_conditions)==Bool
                    proVs, antiVs, proVall, antiVall, opto_fraction,pro_input,anti_input  = run_ntrials_opto(nPro, nAnti; plot_list=1:10, nderivs=nderivs, difforder=difforder, my_params...)                        
                elseif plot_conditions[kk]
                    proVs, antiVs, proVall, antiVall, opto_fraction,pro_input,anti_input  = run_ntrials_opto(nPro, nAnti; plot_list=1:10, nderivs=nderivs, difforder=difforder, my_params...)                        
                else
                    proVs, antiVs, proVall, antiVall, opto_fraction,pro_input,anti_input  = run_ntrials_opto(nPro, nAnti; nderivs=nderivs, difforder=difforder, my_params...)                                                
                end

                hitsP  = 0.5*(1 + tanh.((proVs[1,:]-proVs[4,:,])/theta1))
                diffsP = tanh.((proVs[1,:,]-proVs[4,:])/theta2).^2
                hitsA  = 0.5*(1 + tanh.((antiVs[4,:]-antiVs[1,:,])/theta1))
                diffsA = tanh.((antiVs[4,:,]-antiVs[1,:])/theta2).^2
               
                # set up storage  
                hP[kk,nopto] = mean(hitsP);
                hA[kk,nopto] = mean(hitsA);
                dP[kk,nopto] = mean(diffsP);
                dA[kk,nopto] = mean(diffsA);
                hBP[kk,nopto] = sum(proVs[1,:] .>= proVs[4,:,])/nPro;
                hBA[kk,nopto] = sum(proVs[4,:] .>  proVs[1,:,])/nAnti;

                if nPro>0 && nAnti>0
                    cost1s[kk,nopto] = (nPro*(mean(hitsP) - opto_targets[kk,1]).^2  + nAnti*(mean(hitsA) - opto_targets[kk,2]).^2)/(nPro+nAnti)
                    cost2s[kk,nopto] = -cbeta*(nPro*mean(diffsP) + nAnti*mean(diffsA))/(nPro+nAnti)
                elseif nPro>0
                    cost1s[kk,nopto] = (mean(hitsP) - opto_targets[kk,1]).^2
                    cost2s[kk,nopto] = -cbeta*mean(diffsP)
                else
                    cost1s[kk,nopto] = (mean(hitsA) - opto_targets[kk,2]).^2
                    cost2s[kk,nopto] = -cbeta*mean(diffsA)
                end

            end
        end
    end
    end
    
    cost1 = mean(cost1s)
    cost2 = mean(cost2s)
    if model_details
        return cost1 + cost2, cost1s, cost2s, hP,hA,dP,dA,hBP,hBA, proVall, antiVall, opto_fraction, pro_input, anti_input
    else
        return cost1 + cost2, cost1s, cost2s, hP,hA,dP,dA,hBP,hBA
    end
end

# find the best run from this farm
farmName = "LA";
farmnum=3;

# load farm and do a run
F = matread("goodfarms/farm_"*farmName*lpad(farmnum[1],4,0)*".mat")
model_params = symbol_key_ize(F["model_params"])

####### WARNING, WARNING, WARNING, WARNING. THIS IS A HACK! THIS IS A HACK!
# because julia is dumb, and its 1:30am, I am doing things this way. I am setting the condition I want to get the details for as the last opto condition to run. Dont change anything unless you understand what the next line of code is doing. 
model_params[:opto_periods][5,:] = model_params[:opto_periods][1,:];
train_func =  (;params...) -> JJ_opto(model_params[:nPro],model_params[:nAnti]; rule_and_delay_periods=F["rule_and_delay_periods"], theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=F["post_target_periods"], seedrand=F["sr"], cbeta=F["cb"], verbose=true,plot_conditions=[false,false,false,false,true],model_details=true,  merge(make_dict(F["args"],F["pars"], merge(model_params, Dict(params))))...)

t_opto_scost, t_opto_scost1, t_opto_scost2, t_opto_hitsP,t_opto_hitsA, t_opto_diffsP, t_opto_diffsA, t_opto_bP, t_opto_bA, proVall, antiVall, opto_fraction, pro_input, anti_input = train_func(;:start_pro=>[-0.5,-0.5,-0.5,-0.5],:start_anti=>[-0.5,-0.5,-0.5,-0.5]);

include("CARLOS_COMPUTE_DELAY_ENCODING.jl")
