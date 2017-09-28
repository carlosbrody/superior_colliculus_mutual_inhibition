# include opto inactivation functionality
#opto_periods: a matrix with start and stop times of each opto inactivation to sample, defaults to [-1 -1], which is no inactivation
#opto_strength: the inactivation fraction (default = 1 is no inactivation), 0 would be complete inactivation

function JJ_opto_nll(nPro, nAnti, data; opto_targets=[0.9 0.7], theta1=0.025, theta2=0.035, cbeta=0.003, verbose=false, pre_string="", zero_last_sigmas=0, seedrand=NaN, rule_and_delay_periods = [0.4], target_periods = [0.1], post_target_periods = [0.5], opto_periods = [-1 -1],opto_strength=1, nderivs=0, difforder=0, model_params...) #set opto defaults!

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

#                    Ppro = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
#                    Panti = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
#                    Vpro = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
#                    Vanti = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
#                    PproData = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
#                    PantiData = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
#


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
                proVs, antiVs = run_ntrials_opto(nPro, nAnti; nderivs=nderivs, difforder=difforder, my_params...)

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
                # things I need, nProData, proData, PproData, nAntiData, antiData, PantiData
                # nProData is a matrix: optoconditions x 1  of total number of trials
                # proData is a matrix: optoconditions x 1   of total hits
                    
                    # compute some values
                    Ppro = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
                    Panti = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
                    Vpro = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
                    Vanti = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
                    PproData = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)
                    PantiData = ForwardDiffZeros(1,1, nderivs=nderivs, difforder=difforder)

                    PproData  = data[:proData]./data[:nProData];
                    PantiData = data[:antiData]./data[:nAntiData];
                    Ppro  = (1/(nPro  + data[:nProData][kk] ))*(nPro*mean(hitsP)  + data[:proData][kk] );
                    Panti = (1/(nAnti + data[:nAntiData][kk]))*(nAnti*mean(hitsA) + data[:antiData][kk]);

                    # compute the cost functions
                    Vpro = (mean(hitsP)-Ppro)^2/(Ppro*(1-Ppro))    + data[:nProData][kk]*(PproData[kk] - Ppro)^2/Ppro      + data[:nProData][kk]*(PproData[kk] - Ppro)^2/(1-Ppro);
                    Vanti= (mean(hitsA)-Panti)^2/(Panti*(1-Panti)) + data[:nAntiData][kk]*(PantiData[kk] - Panti)^2/Panti  + data[:nAntiData][kk]*(PantiData[kk] - Panti)^2/(1-Panti);

                    cost1s[kk,nopto] = Vpro + log(Vpro) + Vanti + log(Vanti);
                    if kk ==1
                    cost2s[kk,nopto] = -cbeta*(nPro*mean(diffsP) + nAnti*mean(diffsA))/(nPro+nAnti)
                    end
                elseif nPro>0
                    cost1s[kk,nopto] = (mean(hitsP) - opto_targets[kk,1]).^2
                    if kk == 1
                    cost2s[kk,nopto] = -cbeta*mean(diffsP)
                    end
                else
                    cost1s[kk,nopto] = (mean(hitsA) - opto_targets[kk,2]).^2
                    if kk == 1
                    cost2s[kk,nopto] = -cbeta*mean(diffsA)
                    end
                end
            end
        end
    end
    end
    
    cost1 = sum(cost1s)
    cost2 = sum(cost2s)
   
    return cost1 + cost2, cost1s, cost2s, hP,hA,dP,dA,hBP,hBA
end

