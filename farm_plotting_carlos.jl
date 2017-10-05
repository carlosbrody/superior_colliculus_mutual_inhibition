# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


# If running this in the notebook, you'll have to run the cell below defining load_run() 
# before you can run this one

# for fname in filter(x -> startswith(x, "farm_"), readdir("goodfarms"))

function plot_farm(fname)
    model_params, F, nPro, nAnti = load_run(fname, farmdir="goodfarms");

    # proFullVs = Array{Array{Float64}}(3,1)

    proVs, antiVs, proFullVs, antiFullVs = 
          run_ntrials(10, 10; plot_list=[1:10;], opto_times=model_params[:opto_periods][1,:], model_params...);
    figure(1); subplot(3,1,1); title("PRO -- control")
    figure(2); subplot(3,1,1); title("ANTI -- control")
    proVs, antiVs, proFullVs, antiFullVs = 
          run_ntrials(10, 10; plot_list=[1:10;], opto_times=model_params[:opto_periods][4,:], 
            profig=3, antifig=4, model_params...);
    figure(3); subplot(3,1,1); title("PRO -- delay")
    figure(4); subplot(3,1,1); title("ANTI -- delay")
    proVs, antiVs, proFullVs, antiFullVs = 
          run_ntrials(10, 10; plot_list=[1:10;], opto_times=model_params[:opto_periods][5,:], 
            profig=5, antifig=6, model_params...);
    figure(5); subplot(3,1,1); title("PRO -- choice")
    figure(6); subplot(3,1,1); title("ANTI -- choice")
end

function plot_farm2(fname)
    model_params, F, nPro, nAnti = load_run(fname, farmdir="goodfarms");

    # proFullVs = Array{Array{Float64}}(3,1)

    figure(1); clf();
    ax_set = Dict("pro_Vax"=>subplot(5,3,1), "pro_Dax"=>subplot(5,3,4), 
        "anti_Vax"=>subplot(5,3,10), "anti_Dax"=>subplot(5,3,13))
    proVs = antiVs, proFullVs, antiFullVs = 
        run_ntrials(10, 10; plot_list=[1:10;], opto_times=model_params[:opto_periods][1,:], 
        ax_set=ax_set, profig=1, antifig=1, clearfig=false, plot_Us=false, model_params...);
    
    ax_set = Dict("pro_Vax"=>subplot(5,3,2), "pro_Dax"=>subplot(5,3,5), 
        "anti_Vax"=>subplot(5,3,11), "anti_Dax"=>subplot(5,3,14))
    proVs, antiVs, proFullVs, antiFullVs = 
        run_ntrials(10, 10; plot_list=[1:10;], opto_times=model_params[:opto_periods][4,:], 
        ax_set=ax_set, profig=1, antifig=1, clearfig=false, plot_Us=false, model_params...);

    ax_set = Dict("pro_Vax"=>subplot(5,3,3), "pro_Dax"=>subplot(5,3,6), 
        "anti_Vax"=>subplot(5,3,12), "anti_Dax"=>subplot(5,3,15))
    proVs, antiVs, proFullVs, antiFullVs = 
          run_ntrials(10, 10; plot_list=[1:10;], opto_times=model_params[:opto_periods][5,:], 
        ax_set=ax_set, profig=1, antifig=1, clearfig=false, plot_Us=false, model_params...);
    
end

function histoit(F)
    hP = zeros(5,1)
    hA = zeros(5,1)

    for i=1:5,
        VR = F["runs_pro"][:,end,:,i,1][1,:]
        VL = F["runs_pro"][:,end,:,i,1][4,:]

        hP[i] = mean(sign(VR-VL))

        VR = F["runs_anti"][:,end,:,i,1][1,:]
        VL = F["runs_anti"][:,end,:,i,1][4,:]

        hA[i] = mean(sign(VL-VR))
    end
    figure(7); clf()
    bar(1:3, hP[[1,4,5]], 0.25)
    bar((1:3)+0.25, hA[[1,4,5]], 0.25)
end


# DON'T MODIFY THIS FILE -- the source is in file ProAnti.ipynb


pygui(true)

for resname in filter(x -> startswith(x, "res_"), readdir("../for_carlos_without_runs/goodfarms"))
    F = matread("../for_carlos_without_runs/goodfarms/" * resname)
    histoit(F)
    fname = F["orig_file"]
    plot_farm(fname)
    
    @printf("==> This was %s\n\n", fname)
    @printf("Type anything to go on to next run, type q to quit\n")
    q = readline()
    if startswith(q, "q") || startswith(q, "Q")
        break
    end
end


