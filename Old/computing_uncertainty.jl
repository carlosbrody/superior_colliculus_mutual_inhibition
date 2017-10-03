n = 200;
opto_hits_all = zeros(2,4,100);
for i=1:100
println(i)
new_sr = convert(Int64, round(time()))
test_func =  (;params...) -> JJ_opto_plot(n,n; rule_and_delay_periods=F["rule_and_delay_periods"], theta1=model_params[:theta1], theta2=model_params[:theta2], post_target_periods=F["post_target_periods"], seedrand=new_sr, cbeta=F["cb"], verbose=true,plot_conditions=[false, false, false,true,false],merge(make_dict(F["args"],F["pars"], merge(model_params, Dict(params))))...)


opto_scost, opto_scost1, opto_scost2, opto_hitsP,opto_hitsA, opto_diffsP, opto_diffsA, opto_bP, opto_bA = test_func();
 
opto_hits_all[:,:,i] = opto_hitsP;
end

diffs = zeros(1,8);
c = 0
for j=1:2
for k=1:4
c +=1
x = opto_hits_all[j,k,:];
std_sample = std(x);
errbars = zeros(1,100);
for i=1:100
    errbars[i] = sqrt((1/n)*x[i]*(1-x[i]));
end
[mean(errbars) std_sample]
diffs[c] = mean(errbars)/std_sample
end
end

