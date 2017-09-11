
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
function run_dynamics( ;target1=0.75,target2=0.75,target3=0.75,target4=0.7,target5=0.55,target6=0.7,
    vwi=9, hwi=0.25, pro_bias=0.2135, opto_effect=0.9,
    delay_input=0.25, light_input=3, noisefr=0.005, threshold=0.18, ntrials=50, random_seed=321,
    tau=17.6, dt=10, start_U = [-25, -25, -25, -25],
    g_leak = 1, U_rest = 0, theta1 = 5, beta1 = 50, theta2=0.15, theta3=0.15,
    cue_period = 200, delay_period = 200, choice_period = 200, nderivs=0, difforder=0)




vec_ans_out = ForwardDiffZeros(ntrials, 6; nderivs=nderivs, difforder=difforder)
vec_reac = ForwardDiffZeros(ntrials,6; nderivs=nderivs, difforder=difforder)




for jjj in [1:6;] #trial types: pro, pro delay, pro chioce, anti, anti delay, anti choice


    if (jjj==1)||(jjj==2)||(jjj==3) 
        trial_type="pro"
    else     
        trial_type="anti"
    end



    if (jjj==1)||(jjj==4)
        opto_delay=1;
        opto_choice=1;
    elseif (jjj==2)||(jjj==5)
        opto_delay=opto_effect;
        opto_choice=1;
    elseif (jjj==3)||(jjj==6)
        opto_delay=1;
        opto_choice=opto_effect;
    end



    srand(random_seed)



    for iii in [1:ntrials;]



        t = [0 : dt : cue_period + delay_period + choice_period;]


        V = ForwardDiffZeros(4, length(t); nderivs=nderivs, difforder=difforder)
        U = ForwardDiffZeros(4, length(t); nderivs=nderivs, difforder=difforder)

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


            if t[i] < cue_period+delay_period
                V[:,i]=V[:,i]*opto_delay
            elseif t[i] < cue_period+delay_period+choice_period
                V[:,i]=V[:,i]*opto_choice
            end

            V[:,i] = V[:,i] + noisefr*randn(4)


        end



        if trial_type=="anti"    
            answer_out  = 0.5*(1 + tanh.((V[3,end]  - V[1,end])/theta2))

        elseif trial_type == "pro"
            answer_out  = 0.5*(1 + tanh.((V[1,end]  - V[3,end])/theta2))

        else
            error("invalid trial type")
        end



        #compute reaction time

        reac=NaN;
        for i in [161:length(t)-15;]
           val1=mean(V[1,i-15:i+15]);
           val3=mean(V[3,i-15:i+15]);
           if(abs(val1-val3)>threshold)
               reac=t[i];
               break
           end
       end



       vec_ans_out[iii,jjj]=answer_out;
       vec_reac[iii,jjj]=reac;

   end

end

x=pro_bias
prowin=(1-1*tanh(x*1000)).*(-x+0.007)

cost = (mean(vec_ans_out[:,1]) - target1)^2 + (mean(vec_ans_out[:,2]) - target2)^2 + 
(mean(vec_ans_out[:,3]) - target3)^2 + (mean(vec_ans_out[:,4]) - target4)^2 + 
(mean(vec_ans_out[:,5]) - target5)^2 + (mean(vec_ans_out[:,6]) - target6)^2 + prowin



vec_ans_out=mean(vec_ans_out,1)

return cost, vec_ans_out, vec_reac



end









# # #### LET'S TRY IT!

# cost, answer_out, reac = run_dynamics(ntrials=100)

# println(cost)

# println(answer_out)

# error("ciao")




# params1=["vwi", "hwi", "pro_bias", "delay_input"];
# # params2=[9.0,0.25,0.2135,0.25];
# params2=[4.0,4.0,0.2,0.2];

# out = DiffBase.GradientResult(params2)  # out must be same length as whatever we will differentiate w.r.t.
# keyword_gradient!(out, (;pars...) -> run_dynamics(;pars...)[1], params1, params2)  # note initial values must be floats
# grad = DiffBase.gradient(out)
# cost    = DiffBase.value(out)

# println(grad)
# println(cost)


# error("ciao")









eta = 0.5;




params1=["vwi", "hwi", "pro_bias", "delay_input"];
params2=[1.0,1.0,0.2,0.2];
# params2=[9.0,0.25,0.2135,0.25];
# params2=[4.0,4.0,0.2,0.2];



out = DiffBase.GradientResult(params2)  # out must be same length as whatever we will differentiate w.r.t.
keyword_gradient!(out, (;pars...) -> run_dynamics(;pars...)[1], params1, params2)  # note initial values must be floats
grad = DiffBase.gradient(out)
cost    = DiffBase.value(out)

badstuff,results=run_dynamics(vwi=params2[1],hwi=params2[2],pro_bias=params2[3],delay_input=params2[4])



i=0; 
while eta > 1e-6


    i=i+1
    new_params2 = params2 - eta*grad


    out = DiffBase.GradientResult(new_params2)  # out must be same length as whatever we will differentiate w.r.t.
    keyword_gradient!(out, (;pars...) -> run_dynamics(;pars...)[1], params1, new_params2)  # note initial values must be floats
    grad = DiffBase.gradient(out)
    new_cost    = DiffBase.value(out)


    new_cost2,new_results=run_dynamics(vwi=new_params2[1],hwi=new_params2[2],pro_bias=new_params2[3],delay_input=new_params2[4])    
    if abs(new_cost-new_cost2)>0.0001
        println((new_cost-new_cost2)/new_cost)
        error("yyy")
    end


    new_grad=grad;

    if new_cost < cost

     params2 = new_params2
     cost   = new_cost
     grad   = new_grad
     results = new_results
     eta = eta*1.1
 else    
     eta = eta/2
 end



 if rem(i, 1)==0
   @printf "%d: eta=%f, cost=%.5f, params=[%.3f, %.3f, %.3f, %.3f]\n" i eta cost params2[1] params2[2] params2[3] params2[4]
   # println("eta")
   # println(eta)
   # println("cost")
   # println(cost)
   # println("params")
   # println(params2)
   println("GRADIENT")
   println(grad)
   println("RESULTS")
   println(results)
   println("*********************")
   println("*********************")
end




end





