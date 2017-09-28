# include opto inactivation functionality
#opto_periods: a matrix with start and stop times of each opto inactivation to sample, defaults to [-1 -1], which is no inactivation
#opto_strength: the inactivation fraction (default = 1 is no inactivation), 0 would be complete inactivation

function JJ_opto(nPro, nAnti; opto_targets=[0.9 0.7], theta1=0.025, theta2=0.035, cbeta=0.003, verbose=false, pre_string="", zero_last_sigmas=0, seedrand=NaN, rule_and_delay_periods = [0.4], target_periods = [0.1], post_target_periods = [0.5], opto_periods = [-1 -1],opto_strength=1, nderivs=0, difforder=0, model_params...) #set opto defaults!

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
                    cost1s[kk,nopto] = (nPro*(mean(hitsP) - opto_targets[kk,1]).^2  + nAnti*(mean(hitsA) - opto_targets[kk,2]).^2)/(nPro+nAnti)
                    cost2s[kk,nopto] = -cbeta*(nPro*mean(diffsP) + nAnti*mean(diffsA))/(nPro+nAnti)
                elseif nPro>0
                    cost1s[kk,nopto] = (mean(hitsP) - opto_targets[kk,1]).^2
                    cost2s[kk,nopto] = -cbeta*mean(diffsP)
                else
                    cost1s[kk,nopto] = (mean(hitsA) - opto_targets[kk,2]).^2
                    cost2s[kk,nopto] = -cbeta*mean(diffsA)
                end

              #  totHitsP  += mean(hitsP);  totHitsA  += mean(hitsA); 
              #  totDiffsP += mean(diffsP); totDiffsA += mean(diffsA);
            end
        end
    end
    end
    
    cost1 = mean(cost1s)
    cost2 = mean(cost2s)
#    print(cost1s)
    
#    print(cost2s)
#    hitsP = totHitsP/n; hitsA = totHitsA/n; diffsP = totDiffsP/n; diffsA = totDiffsA/n
    
    
#    if verbose
#        @printf("%s", pre_string)
#        @printf("     -- cost=%g,   cost1=%g, cost2=%g\n", 
#            convert(Float64, cost1+cost2), convert(Float64, cost1), convert(Float64, cost2))
#        if nPro>0 && nAnti>0
#            @printf("     -- mean(hitsP)=%g, mean(diffsP)=%g mean(hitsA)=%g, mean(diffsA)=%g\n", 
#                convert(Float64, mean(hitsP)), convert(Float64, mean(diffsP)),
#                convert(Float64, mean(hitsA)), convert(Float64, mean(diffsA)))
#        elseif nPro>0
#            @printf("     -- mean(hitsP)=%g, mean(diffsP)=%g (nAnti=0)\n", 
#                convert(Float64, mean(hitsP)), convert(Float64, mean(diffsP)))
#        else
#            @printf("     -- (nPro=0) mean(hitsA)=%g, mean(diffsA)=%g\n", 
#                convert(Float64, mean(hitsA)), convert(Float64, mean(diffsA)))
#        end        
#    end 
    
    return cost1 + cost2, cost1s, cost2s, hP,hA,dP,dA,hBP,hBA
end

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
    
    return cost1 + cost2, cost1, cost2, hitsP, hitsA, diffsP, diffsA
end


# returns a vector of the inhibition fraction for each node given the opto period and opto strength
# right now, it always computes bilateral inactivations
function make_opto_input(nsteps; dt = 0.02, opto_strength=1, opto_period=[-1 1],nderivs=0,difforder=0, other_unused_params...);
       
    strength = ForwardDiffZeros(1,1,nderivs=nderivs, difforder=difforder);
    strength[1] = strength[1] + opto_strength;

    if opto_period[2] < 0 # no opto, return all ones
        opto_fraction = 1+ForwardDiffZeros(4,nsteps,nderivs=nderivs, difforder=difforder);
#        print("0 0 out of ")
#        print(nsteps*dt)
#        print(" with command ")
#        print(opto_period)  
#        println()
    else # opto, compute fraction vector
        # check for variable inputs
        other_unused_params = Dict(other_unused_params);
        end_rd =     other_unused_params[:rule_and_delay_period];
        end_target = other_unused_params[:rule_and_delay_period]+other_unused_params[:target_period]; 
        if opto_period[1]==100
            opto_period[1] = end_rd;
        elseif opto_period[1]==200
            opto_period[1] = end_target;
        end        
        if opto_period[2]==100
            opto_period[2] = end_rd;
        elseif opto_period[2]==200
            opto_period[2] = end_target;
        end        

        # check for start and end of trial flags
        start_step = round(opto_period[1]/dt);
        if start_step < 1; start_step = 1; end
        end_step = round(opto_period[2]/dt);
        if end_step > nsteps; end_step = nsteps; end
        if end_step < 1; end_step = 1; end

#        print(start_step*dt)
#        print(" ")
#        print(end_step*dt)
#        print(" out of ")
#        print(nsteps*dt)
#        print(" with command ")
#        print(opto_period)
#        println()       
 
        # compute the fraction
        opto_fraction = 1+ForwardDiffZeros(4,nsteps,nderivs=nderivs, difforder=difforder);
        opto_fraction[:,Int(start_step):Int(end_step)] = repeat(strength,outer=[4,Int(end_step)-Int(start_step)+1]);
    end
    
    return opto_fraction
end


function run_ntrials_opto(nPro, nAnti; plot_list=[], nderivs=0, difforder=0,start_pro=[-0.3, -0.7, -0.7, -0.3], start_anti=[-0.7, -0.3, -0.3, -0.7],  opto_periods=[-1 -1],model_params...)
    pro_input,  t, nsteps = make_input("Pro" ; model_params...);
    anti_input, t, nsteps = make_input("Anti"; model_params...);
    #### make opto fraction vector
    opto_fraction = make_opto_input(nsteps; nderivs=nderivs,difforder=difforder, model_params...);
    model_params = Dict(model_params)
    sW = model_params[:sW]
    hW = model_params[:hW]
    vW = model_params[:vW]
    dW = model_params[:dW]
    model_params = make_dict(["nsteps", "W"], [nsteps, [sW vW dW hW; vW sW hW dW; dW hW sW vW; hW dW vW sW]], model_params)
    model_params = make_dict(["nderivs", "difforder"], [nderivs, difforder], model_params)
    
    ## include opto fraction in model_params dictionary
    model_params = make_dict(["opto_fraction"], [opto_fraction], model_params);

    proVs  = ForwardDiffZeros(4, nPro, nderivs=nderivs, difforder=difforder)
    antiVs = ForwardDiffZeros(4, nAnti, nderivs=nderivs, difforder=difforder)

    # --- PRO ---
    if length(plot_list)>0; figure(1); clf(); end
    model_params = make_dict(["input"], [pro_input], model_params)
    for i=1:nPro
       # startU = [-0.3, -0.7, -0.7, -0.3]
        startU = start_pro;
        Uend, Vend, U, V = forwardModel_opto(startU, opto_fraction, do_plot=false; model_params...)
        proVs[:,i] = Vend
        if any(plot_list.==i) 
            plot_PA(t, U, V; fignum=1, clearfig=false, model_params...)
            subplot(3,1,1); title("PRO")
        end
    end

    # --- ANTI ---
    if length(plot_list)>0; figure(2); clf(); end
    model_params = make_dict(["input"], [anti_input], model_params)
    for i=1:nAnti
        #startU = [-0.7, -0.3, -0.3, -0.7]
        startU = start_anti;

        Uend, Vend, U, V = forwardModel_opto(startU,opto_fraction, do_plot=false; model_params...)
        antiVs[:,i] = Vend
        if any(plot_list.==i) 
            plot_PA(t, U, V; fignum=2, clearfig=false, model_params...)
            subplot(3,1,1); title("ANTI")
        end
    end
    
    return proVs, antiVs
end

function forwardModel_opto(startU, opto_fraction; dt=0.01, tau=0.1, nsteps=100, input=[], noise=[], W=[0 -5;-5 0], 
    init_add=0, start_add=0, const_add=0, do_plot=false, nderivs=0, difforder=0, clearfig=true, fignum=1,
    dUdt_mag_only=false, sigma=0, g_leak=1, U_rest=0, theta=0, beta=1, 
    warn_if_unused_params=false, other_unused_params...)

    if warn_if_unused_params && length(other_unused_params)>0
        @printf("\n\n=== forwardModel warning, had unused params ")
        for k in keys(Dict(other_unused_params))
            @printf("%s, ", k)
        end
    end
    
    my_input = ForwardDiffZeros(size(input,1), size(input,2), nderivs=nderivs, difforder=difforder)
    for i=1:prod(size(input)); my_input[i] = input[i]; end
    input = my_input;
   
    my_opto = ForwardDiffZeros(size(opto_fraction,1), size(opto_fraction,2), nderivs=nderivs, difforder=difforder)
    for i=1:prod(size(opto_fraction)); my_opto[i] = opto_fraction[i]; end
    opto_fraction = my_opto;
    opto_fraction = opto_fraction.*(1+ForwardDiffZeros(size(opto_fraction,1), nsteps, nderivs=nderivs, difforder=difforder))
 
    nunits = length(startU)
    if size(startU,2) > size(startU,1)
        error("startU must be a column vector")
    end
    
    # --- formatting input ---
    if ~(typeof(input)<:Array) || prod(size(input))==1  # was a scalar
        input = input[1]*(1+ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder))
    elseif length(input)==0 # was the empty matrix
        input = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    elseif size(input,2)==1     # was a column vector
        input = input*(1+ForwardDiffZeros(1, nsteps, nderivs=nderivs, difforder=difforder))
    end    
    # --- formatting noise ---
    if ~(typeof(noise)<:Array) || prod(size(noise))==1  # was a scalar
        noise = noise*(1+ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder))
    elseif length(noise)==0 # was the empty matrix
        noise = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    elseif size(noise,2)==1     # was a column vector
        noise = noise*(1+ForwardDiffZeros(1, nsteps, nderivs=nderivs, difforder=difforder))
    end    
    
    U = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    V = ForwardDiffZeros(nunits, nsteps, nderivs=nderivs, difforder=difforder)
    
    if ~(typeof(W)<:Array); W = [W]; end

    W     = reshape(W, nunits, nunits)
    U     = reshape(U, nunits, nsteps)
    V     = reshape(V, nunits, nsteps)
    input = reshape(input, nunits, nsteps)
    noise = reshape(noise, nunits, nsteps)

    input[:,1] += init_add
    input      += const_add

    #@printf("size(U) is (%d,%d), and size(startU) is (%d,%d) and size(noise) is (%d,%d)", 
    #    size(U,1), size(U,2), size(startU,1), size(startU,2), size(noise,1), size(noise,2))
    # @printf("U[1]=%g, noise[1]=%g\n", startU, noise[1])
    U[:,1] = startU + noise[:,1] + start_add; # @printf("Resulting U=%g\n", U[1])
    V[:,1] = g((U[:,1]-theta)/beta); # @printf("Resulting V=%g\n", V[1])
    V[:,1] = V[:,1].*opto_fraction[:,1];
    
    for i=2:nsteps
        dUdt = g_leak*(U_rest -U[:,i-1]) + W*V[:,i-1] + input[:,i-1]
        if dUdt_mag_only; return sum(dUdt.*dUdt); end;
        # @printf("dUdt=%g\n", dUdt[1])
        # @printf("i=%g\n", i)
        # @printf("noise[2]=%g\n", noise[2])
        U[:,i] = U[:,i-1] + (dt/tau)*dUdt + noise[:,i] + sigma*sqrt(dt)*randn(size(U,1),1)
        # @printf("Resulting U[2]=%g\n", U[2])
        V[:,i] = g((U[:,i]-theta)/beta)
        V[:,i] = V[:,i].*opto_fraction[:,i]
        # @printf("Resulting V[2]=%g\n", V[2])
    end

    if do_plot
        figure(fignum)
        if length(startU)==1
            if clearfig; clf(); end;
            t = (0:nsteps-1)*dt
            plot(t, V[1,:], "b-")
            plot(t[1], V[1,1], "g.")
            plot(t[end], V[1,end], "r.")
            xlabel("t"); ylabel("V1"); ylim([-0.01, 1.01])
        elseif length(startU)>=2
            if clearfig; clf(); end;
            plot(V[1,:], V[2,:], "b-")
            plot(V[1,1], V[2,1], "g.")
            plot(V[1,end], V[2,end], "r.")
            xlabel("V1"); ylabel("V2"); 
            xlim([-0.01, 1.01]); ylim([-0.01, 1.01])
        end
    end

    return U[:,end], V[:,end], U, V
end


