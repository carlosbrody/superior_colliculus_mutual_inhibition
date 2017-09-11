
using PyCall
using PyPlot
using ForwardDiff
using DiffBase

pygui(true)


import Base.convert
convert(::Type{Float64}, x::ForwardDiff.Dual) = Float64(x.value)
function convert(::Array{Float64}, x::Array{ForwardDiff.Dual}) 
    y = zeros(size(x)); 
    for i in 1:prod(size(x)) 
        y[i] = convert(Float64, x[i]) 
    end
    return y
end


include("hessian_utils.jl")




#RUN_DYNAMICS FUNCTION; RUNS 1 TRIAL
#PARAMS: vert w; horiz w; pro bias; delay input
function run_dynamics(trial_type, params::Vector ; opto_cue=1, opto_delay=1, opto_choice=1,
    light_input=12, noisefr=0.1, threshold=0.18,
    tau=4.4, dt=0.05, start_U = [-25, -25, -25, -25],
    g_leak = 1, U_rest = 0, theta1 = 5, beta1 = 50, theta2=0.15, theta3=0.15, do_plot = false, fignum=1,
    cue_period = 200, delay_period = 200, choice_period = 50)

vwi = params[1]; hwi = params[2]; pro_bias = params[3]; delay_input = params[4]

t = [0 : dt : cue_period + delay_period + choice_period;]

V = zeros(eltype(params), 4, length(t))   # the eltype(params) is for ForwardDiff
U = zeros(eltype(params), 4, length(t))

U[:,1] = start_U

W = [0 -vwi -hwi 0; -vwi 0 0 -hwi;
-hwi 0 0 -vwi; 0 -hwi -vwi 0]


for i in [2:length(t);]  # the funny semicolon appears to be necessary in Julia

    dUdt = W * V[:,i-1] + g_leak*(U_rest - U[:,i-1])/tau

    if t[i] < cue_period + delay_period
        if trial_type=="anti"
            dUdt[[2,4]] += delay_input
        elseif trial_type == "pro"
            dUdt[[1,3]] += delay_input
        else
            error("invalid trial type")
        end

    elseif t[i] < cue_period + delay_period + choice_period
        dUdt[[1,2]] += light_input
    end
    
    dUdt[[1,3]] += pro_bias



    U[:,i] = U[:,i-1] +  dt*dUdt

    V[:,i] = 0.5*tanh((U[:,i]-theta1)/beta1) + 0.5

    if t[i] < cue_period
        V[:,i]=V[:,i]*opto_cue
    elseif t[i] < cue_period+delay_period
        V[:,i]=V[:,i]*opto_delay
    elseif t[i] < cue_period+delay_period+choice_period
        V[:,i]=V[:,i]*opto_choice
    end

    V[:,i] = V[:,i] + noisefr*randn(4)*sqrt(dt)


end


if do_plot
    figure(fignum);
    h = plot(t, V');
    setp(h[1], color=[0, 0, 1])
    setp(h[2], color=[1, 0, 0])
    setp(h[3], color=[1, 0.5, 0.5])
    setp(h[4], color=[0, 1, 1])

    ax = gca()
    yl = [ylim()[1], ylim()[2]]
    vlines([cue_period, cue_period+delay_period,
        cue_period+delay_period+choice_period],
        0.05, 1.05, linewidth=2)
    if yl[1]<0.02
     yl[1] = -0.02
 end
 if yl[2]>0.98
     yl[2] = 1.02
 end
 ylim(yl)
 grid(true)
end



if trial_type=="anti"
    answer_out  = V[3,end]  - V[1,end]
    # answer_out  = 0.5*(1 + tanh.((V[3,end]  - V[1,end])/theta2))
    answer_dif  = tanh.((V[3,end]  - V[1,end])/theta3).^2

elseif trial_type == "pro"
    answer_out  = V[1,end]  - V[3,end]
    # answer_out  = 0.5*(1 + tanh.((V[1,end]  - V[3,end])/theta2))
    answer_dif  = tanh.((V[1,end]  - V[3,end])/theta3).^2
else
    error("invalid trial type")
end



#compute reaction time

reac=NaN;
for i in [8001:length(t)-15;]
 val1=mean(V[1,i-15:i+15]);
 val3=mean(V[3,i-15:i+15]);
 if(abs(val1-val3)>threshold)
     reac=t[i];
     break
 end
end


return answer_out,answer_dif, reac, t, U, V
end






# # #### LET'S TRY IT!

# params = [36, 1, 0.854, 1]

# answer_out, answer_dif, reac = run_dynamics("pro", params,do_plot=true,fignum=1)

# println(answer_out)
# println(answer_dif)
# println(reac)

# show()








#THIS FUNCTION RUNS RUN_DYNAMICS N TIMES


function run_n_times(trial_type, params::Vector, ntrials=20; opto_cue=1, opto_delay=1, opto_choice=1,
  random_seed=321, light_input=12, noisefr=0.1, threshold=0.18)

srand(random_seed)

#initialize these guys or it doesn't work
vec_ans_out=1;
vec_ans_dif=1;
vec_reac=1;


for i in [1:ntrials;]
 answer_out, answer_dif, reac = run_dynamics(trial_type, params, noisefr=noisefr, light_input=light_input, opto_cue=opto_cue, opto_delay=opto_delay, opto_choice=opto_choice)

 if(i==1)
     vec_ans_out=answer_out;
     vec_ans_dif=answer_dif;
     vec_reac=reac;
 else
     vec_ans_out=[vec_ans_out answer_out];
     vec_ans_dif=[vec_ans_dif answer_dif];
     vec_reac=[vec_reac reac];
 end

end

return vec_ans_out, vec_ans_dif, vec_reac

end








# #### LET'S TRY IT!

# params = [36, 1, 0.854, 1]

# vec_ans_out, vec_ans_dif, vec_reac = run_n_times("anti",params,20)

# println(vec_ans_out)
# println(vec_ans_dif)
# println(vec_reac)






### COST FUNCTION WITH OPTO


function J_opto(params, targets; ntrials=20, noisefr=0.1, random_seed=321, beta2=0, verbose=false)

    srand(random_seed)


    pro_ans_out, pro_ans_dif, pro_reac = run_n_times("pro",params,ntrials)

    anti_ans_out, anti_ans_dif, anti_reac = run_n_times("anti",params,ntrials)

    pro_ans_delay_out, pro_ans_delay_dif, pro_reac_delay = run_n_times("pro",params,ntrials,opto_delay=0.95)

    anti_ans_delay_out, anti_ans_delay_dif, anti_reac_delay = run_n_times("anti",params,ntrials,opto_delay=0.95)


    pro_ans_choice_out, pro_ans_choice_dif, pro_reac_choice = run_n_times("pro",params,ntrials,opto_choice=0.95)

    anti_ans_choice_out, anti_ans_choice_dif, anti_reac_choice = run_n_times("anti",params,ntrials,opto_choice=0.95)



    cost1 = (mean(pro_ans_out) - targets[1])^2 + (mean(anti_ans_out) - targets[2])^2 + 
    (mean(pro_ans_delay_out) - targets[3])^2 + (mean(anti_ans_delay_out) - targets[4])^2 + 
    (mean(pro_ans_choice_out) - targets[5])^2 + (mean(anti_ans_choice_out) - targets[6])^2 



    cost2 = -mean(pro_ans_dif) - mean(anti_ans_dif) - mean(pro_ans_delay_dif) - mean(anti_ans_delay_dif) - mean(pro_ans_choice_dif) - mean(anti_ans_choice_dif)

    cost=cost1 + beta2*cost2



    results=[mean(pro_ans_out) mean(anti_ans_out) mean(pro_ans_delay_out) mean(anti_ans_delay_out) mean(pro_ans_choice_out) mean(anti_ans_choice_out)]


    if verbose
        @printf("cost1=%.3f, cost2=%.3f, total cost=%.3f\n", convert(Float64, cost1), beta2*convert(Float64, cost2), convert(Float64, cost))
        println("RESULTS")
        println(results)
    end

    return cost, results, cost1, cost2


end






# #### LET'S TEST IT!

# params = [36, 1, 0.854, 1]
# targets=[0.8,0.7,0.8,0.5,0.8,0.7]
# ntrials=20;


# out = J_opto(params,targets,noisefr=5,verbose=true)
# cost=out[1]
# results=out[2]
# cost1=out[3]
# cost2=out[4]
# println(cost)
# println(results)
# println(cost1)
# println(cost2)



# out = DiffBase.HessianResult(params)
# ForwardDiff.hessian!(out, x -> J_opto(x, targets, ntrials=ntrials)[1], params)
# cost = DiffBase.value(out)
# grad = DiffBase.gradient(out)
# hess = DiffBase.hessian(out)

# println(cost)
# println(grad)
# println(hess)















params = [15, 15, 0.5, 0.5]


targets = [0.8, 0.7, 0.8, 0.5, 0.8, 0.7]   # Fraction correct in Pro and Anti


ntrials = 20

eta = 10;





J_out=J_opto(params,targets,ntrials=ntrials)
cost=J_out[1]
results=J_out[2]



out = DiffBase.HessianResult(params)
ForwardDiff.hessian!(out, x -> J_opto(x, targets, ntrials=ntrials)[1], params)
cost = DiffBase.value(out)
grad = DiffBase.gradient(out)
hess = DiffBase.hessian(out)

hessgrad=inv(hess)*grad



i=0; 
while eta > 1e-6


    i=i+1
    new_params = params - eta*grad


    J_out=J_opto(new_params,targets,ntrials=ntrials)
    new_cost=J_out[1]
    new_results=J_out[2]



    out = DiffBase.HessianResult(new_params)
    ForwardDiff.hessian!(out, x -> J_opto(x, targets, ntrials=ntrials)[1], new_params)
    new_cost2 = DiffBase.value(out)
    grad = DiffBase.gradient(out)
    hess = DiffBase.hessian(out)

    if abs(new_cost-new_cost2)>0.0001
        println((new_cost-new_cost2)/new_cost)
        error("yyy")
    end

    new_hessgrad=inv(hess)*grad

    if new_cost < cost

       params = new_params
       cost   = new_cost
       hessgrad   = new_hessgrad
       results = new_results
       eta = eta*1.1
   else

       eta = eta/2
   end



   if rem(i, 1)==0
       @printf "%d: eta=%f, cost=%.5f, params=[%.3f, %.3f, %.3f, %.3f]\n" i eta cost params[1] params[2] params[3] params[4]
       println("HESSIAN * GRADIENT")
       println(hessgrad)
       println("RESULTS")
       println(results)
       println("*********************")
   end




end





