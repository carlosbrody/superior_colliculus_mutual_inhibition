#make index plots
using Distributions
examples,results = load("MiniC32_C32_examples_50_long.jld","examples","results");
examplesS = get_synthetic_LR_trials(examples);
dex = results["cost"] .<= -0.0001;
ex = examples[dex,:,:,:,:,:];
exS = examplesS[dex,:,:,:,:,:];
dt = 0.024;
tvec = vec(collect(0:74).*dt);


# AVERAGE OVER POPULATION
# make choice index
choice_dex = zeros(size(ex,1)*size(ex,6)*size(ex,3), size(ex,5));
count = 1;
for i=1:size(ex,1)
    for j=1:size(ex,6)
        choice_dex[count,:] = (ex[i,1,1,1,:,j] - ex[i,1,1,4,:,j]) ;
        count += 1;
        choice_dex[count,:] = (ex[i,1,2,1,:,j] - ex[i,1,2,4,:,j]) ;
        count += 1;
    end
end
mean_choice_dex = mean(abs.(choice_dex),1);
std_choice_dex = std(abs.(choice_dex),1);
#figure();
#plot(tvec, vec(mean_choice_dex),"k")
#plot(tvec, vec(mean_choice_dex + std_choice_dex),"r")
#plot(tvec, vec(mean_choice_dex - std_choice_dex),"r")

rchoice = choice_dex[:,end] .> 0;
muR = mean(choice_dex[rchoice,:],1);
muL = mean(choice_dex[.!rchoice,:],1);
varR = var(choice_dex[rchoice,:],1);
varL = var(choice_dex[.!rchoice,:],1);
choice_dprime = (muR-muL)./(sqrt.(0.5.*( varR+varL  )  ) );

figure(figsize=(4,2));
plot(tvec, vec(choice_dprime),"k")
ylabel("Average Choice d'",fontsize=12)
xlabel("Time",fontsize=12)
plot(vec([1.2 1.2]), ylim(),color=(179/255,179/255,179/255));
xlim(tvec[1],tvec[end])
ylim(0, ylim()[2])
xticks([0, 0.5, 1, 1.5],("0","0.5","1","1.5"),fontsize=12)
yticks(fontsize=12)

# make rule index
rule_dex = zeros(size(ex,1)*size(ex,6)*size(ex,3), size(ex,5));
count = 1;
for i=1:size(ex,1)
    for j=1:size(ex,6)
        rule_dex[count,:] = 0.5.*(ex[i,1,1,1,:,j] + ex[i,1,1,4,:,j]) - 0.5.*(ex[i,1,1,2,:,j] + ex[i,1,1,3,:,j]);
        count += 1;
        rule_dex[count,:] = 0.5.*(ex[i,1,2,1,:,j] + ex[i,1,2,4,:,j]) - 0.5.*(ex[i,1,2,2,:,j] + ex[i,1,2,3,:,j]);
        count += 1;
    end
end
mean_rule_dex = mean(abs.(rule_dex),1);
std_rule_dex = std(abs.(rule_dex),1);
#figure();
#plot(tvec, vec(mean_rule_dex),"k")
#plot(tvec, vec(mean_rule_dex + std_rule_dex),"r")
#plot(tvec, vec(mean_rule_dex - std_rule_dex),"r")

prodex  = collect(1:2:37300);
antidex  = collect(2:2:37300);
muP = mean(abs.(rule_dex[prodex,:]),1);
muA = mean(abs.(rule_dex[antidex,:]),1);
varP = var(abs.(rule_dex[prodex,:]),1);
varA = var(abs.(rule_dex[antidex,:]),1);
rule_dprime = (muP-muA)./(sqrt.(0.5.*( varP+varA  )  ) );

figure(figsize=(4,2));
plot(tvec, vec(-rule_dprime),"k")
ylabel("Average Task d'",fontsize=12)
xlabel("Time",fontsize=12)
plot(vec([1.2 1.2]), ylim(),color=(179/255,179/255,179/255));
xlim(tvec[1],tvec[end])
#ylim(0, ylim()[2])
xticks([0, 0.5, 1, 1.5],("0","0.5","1","1.5"),fontsize=12)
yticks(fontsize=12)

########
########
# Average over solutions
########
########
# make choice index
if false
    choice_dex = zeros(size(ex,1),size(ex,6)*size(ex,3), size(ex,5));
    choice_dprime = zeros(size(ex,1),size(ex,5));
    for i=1:size(ex,1)
        hit = trues(200,1);
        for j=1:size(ex,6)
            choice_dex[i,j,:] = (ex[i,1,1,1,:,j] - ex[i,1,1,4,:,j]) ;
            choice_dex[i,j+100,:] = (ex[i,1,2,1,:,j] - ex[i,1,2,4,:,j]) ;
            hit[j] = ex[i,1,1,1,end,j] > ex[i,1,1,4,end,j];
            hit[j+100] = ex[i,1,2,1,end,j] < ex[i,1,2,4,end,j];
        end
    #    rchoice = vec((choice_dex[i,:,end] .> 0) .& hit );
    #    lchoice = vec((choice_dex[i,:,end] .< 0) .& hit );
        rchoice = vec((choice_dex[i,:,end] .> 0) );
        lchoice = vec((choice_dex[i,:,end] .< 0) );
    
        muR = mean(choice_dex[i,rchoice,:],1);
        muL = mean(choice_dex[i,lchoice,:],1);
        varR = var(choice_dex[i,rchoice,:],1);
        varL = var(choice_dex[i,lchoice,:],1);
        choice_dprime[i,:] = (muR-muL)./(sqrt.(0.5.*( varR+varL  )  ) );    
    end
else
    choice_dex = zeros(size(ex,1),size(ex,6)*size(ex,3), size(ex,5));
    choice_dprime = zeros(size(ex,1),size(ex,5));
    for i=1:size(ex,1)
 #       hit = trues(200,1);
        for j=1:size(ex,6)
            choice_dex[i,j,:] = (exS[i,1,1,1,:,j] - exS[i,1,1,4,:,j]) ;
            choice_dex[i,j+100,:] = (exS[i,1,2,1,:,j] - exS[i,1,2,4,:,j]) ;
#            hit[j] = ex[i,1,1,1,end,j] > ex[i,1,1,4,end,j];
#            hit[j+100] = ex[i,1,2,1,end,j] < ex[i,1,2,4,end,j];
        end
        rchoice = vec((choice_dex[i,:,end] .> 0) );
        lchoice = vec((choice_dex[i,:,end] .< 0) );
    
        muR = mean(choice_dex[i,rchoice,:],1);
        muL = mean(choice_dex[i,lchoice,:],1);
        varR = var(choice_dex[i,rchoice,:],1);
        varL = var(choice_dex[i,lchoice,:],1);
        choice_dprime[i,:] = (muR-muL)./(sqrt.(0.5.*( varR+varL  )  ) );    
    end
end
# make rule index
rule_dex = zeros(size(ex,1),size(ex,6)*size(ex,3), size(ex,5));
rule_dprime = zeros(size(ex,1), size(ex,5));
for i=1:size(ex,1)
    for j=1:size(ex,6)
        rule_dex[i,j,:] = 0.5.*(ex[i,1,1,1,:,j] + ex[i,1,1,4,:,j]) - 0.5.*(ex[i,1,1,2,:,j] + ex[i,1,1,3,:,j]);
        rule_dex[i,j+100,:] = 0.5.*(ex[i,1,2,1,:,j] + ex[i,1,2,4,:,j]) - 0.5.*(ex[i,1,2,2,:,j] + ex[i,1,2,3,:,j]);
    end
    prodex  = collect(1:1:size(ex,6));
    antidex  = collect(101:1:size(ex,6)*2);
    muP = mean(rule_dex[i,prodex,:],1);
    muA = mean(rule_dex[i,antidex,:],1);
    varP = var(rule_dex[i,prodex,:],1);
    varA = var(choice_dex[i,antidex,:],1);
    rule_dprime[i,:] = (muP-muA)./(sqrt.(0.5.*( varP+varA  )  ) );
end



figure(figsize=(4,1.5));
plot(tvec, vec(mean(choice_dprime,1)),"k")
plot(tvec, vec(mean(choice_dprime,1)-std(choice_dprime,1)./sqrt(373)),color=(179/255,179/255,179/255))
plot(tvec, vec(mean(choice_dprime,1)+std(choice_dprime,1)./sqrt(373)),color=(179/255,179/255,179/255))
ylabel("Avg. Choice d'",fontsize=12)
xlabel("Time",fontsize=12)
ylim(0, ylim()[2])
plot(vec([1.2 1.2]), ylim(),color=(179/255,179/255,179/255));
xlim(tvec[1],tvec[end])
xticks([0, 0.5, 1, 1.5],("0","0.5","1","1.5"),fontsize=12)
yticks(fontsize=12)

figure(figsize=(4,1.5));
meanD = vec(mean(choice_dprime,1));
lbD = vec(mean(choice_dprime,1)-std(choice_dprime,1)./sqrt(373));
ubD = vec(mean(choice_dprime,1)+std(choice_dprime,1)./sqrt(373));
meanD_p = cdf(Normal(0,1), meanD./sqrt(2));
lbD_p = cdf(Normal(0,1), lbD./sqrt(2));
ubD_p = cdf(Normal(0,1), ubD./sqrt(2));
plot(tvec, meanD_p.*100,"k")
plot(tvec, lbD_p.*100,color=(179/255,179/255,179/255))
plot(tvec, ubD_p.*100,color=(179/255,179/255,179/255))
ylabel("Avg. Choice decoding",fontsize=12)
xlabel("Time",fontsize=12)
ylim(50, 100)
plot(vec([1.2 1.2]), ylim(),color=(179/255,179/255,179/255));
xlim(tvec[1],tvec[end])
xticks([0, 0.5, 1, 1.5],("0","0.5","1","1.5"),fontsize=12)
yticks(fontsize=12)

figure(figsize=(4,1.5));
plot(tvec, vec(mean(rule_dprime,1)),"k")
plot(tvec, vec(mean(rule_dprime,1) - std(rule_dprime,1)./sqrt(373) ),color=(179/255,179/255,179/255))
plot(tvec, vec(mean(rule_dprime,1) + std(rule_dprime,1)./sqrt(373) ),color=(179/255,179/255,179/255))
ylabel("Avg. Task d'",fontsize=12)
xlabel("Time",fontsize=12)
ylim(0, ylim()[2])
plot(vec([1.2 1.2]), ylim(),color=(179/255,179/255,179/255));
xlim(tvec[1],tvec[end])
xticks([0, 0.5, 1, 1.5],("0","0.5","1","1.5"),fontsize=12)
yticks(fontsize=12)

figure(figsize=(4,1.5));
meanR =vec(mean(rule_dprime,1));
lbR =  vec(mean(rule_dprime,1) - std(rule_dprime,1)./sqrt(373) );
ubR =  vec(mean(rule_dprime,1) + std(rule_dprime,1)./sqrt(373) );
meanR_p = cdf(Normal(0,1), meanR./sqrt(2));
lbR_p = cdf(Normal(0,1), lbR./sqrt(2));
ubR_p = cdf(Normal(0,1), ubR./sqrt(2));
plot(tvec, meanR_p.*100,"k")
plot(tvec, lbR_p.*100,color=(179/255,179/255,179/255))
plot(tvec, ubR_p.*100,color=(179/255,179/255,179/255))
ylabel("Avg. Task decoding",fontsize=12)
xlabel("Time",fontsize=12)
ylim(50, 100)
plot(vec([1.2 1.2]), ylim(),color=(179/255,179/255,179/255));
xlim(tvec[1],tvec[end])
xticks([0, 0.5, 1, 1.5],("0","0.5","1","1.5"),fontsize=12)
yticks(fontsize=12)


