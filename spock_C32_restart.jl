include("pro_anti.jl")
# Script for restarting fits of C32 that stopped due to time limit


# Figure out which file to load
if length(ARGS)>0  &&  ~isnull(tryparse(Int64, ARGS[1])); 
    my_run_number = parse(Int64, ARGS[1]); 
else                                                      
#    my_run_number = 1; # I am process my_run_number
   error("ARGS was weird")
end
moddex = mod(my_run_number,8) + 1;
if moddex <= 4
    dex1 = "01";
    dex2 = moddex;
else
    dex1 = "02";
    dex2 = moddex - 4;
end
ndex = fld(my_run_number-1,8)+1;
loadname = "../Farms_C32/farm_C32_spock-brody"*dex1*"-0"*string(dex2)*"_"*lpad(string(ndex),4,0)*".jld";


# Figure out report filename
reports_dir = "../Reports"
if !isdir(reports_dir); mkdir(reports_dir); end
report_file = reports_dir * "/" * "refit_"*dex1*"-0"*string(dex2)*"_"*lpad(string(ndex),4,0)

   
# load file
f = load(loadname)


# check file
vals, vecs = eig(f["ftraj3"][2,end])
if all(vals .> 0) && isreal(vals)
    # no need to refit
    append_to_file(report_file, @sprintf("\n\n**** No need to refit **** %s ---\n\n", Dates.format(now(), "e, dd u yyyy HH:MM:SS")))
else
    # need to refit
    append_to_file(report_file, @sprintf("\n\n**** Need to refit **** %s ---\n\n", Dates.format(now(), "e, dd u yyyy HH:MM:SS")))

    extra_pars  = f["extra_pars"]
    mypars      = f["mypars"]
    func_chatty =  (;params...) -> JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=true, verbose_file = report_file, merge(merge(mypars, extra_pars), Dict(params))...)[1]

    
    ## Enter loop of optimizing, and then saving the result
    pars3       = f["pars3"]
    args        = f["args"]
    bbox        = f["bbox"]
    cost3       = f["cost3"];
    old_cost3   = cost3+10;
    while cost3 + 1.0e-12 < old_cost3

         append_to_file(report_file, @sprintf("\n\n**** New training iteration **** %s ---\n\n", Dates.format(now(), "e, dd u yyyy HH:MM:SS")))    
        old_cost3 = cost3;
    
        pars3, traj3, cost3, cpm_traj3, ftraj3 = bbox_Hessian_keyword_minimization(pars3, args, bbox, func_chatty,  verbose_timestamp=true, start_eta = 0.01, tol=1e-12, verbose_file=report_file, verbose=true, verbose_every=10, maxiter=50)
    
        # Save Intermediate result
        f["pars3"] = pars3
        f["traj3"] = traj3
        f["cost3"] = cost3
        f["cpm_traj3"] =  cpm_traj3
        f["ftraj3"] = ftraj3
        save(loadname,f)
    end

    # do final eval
    append_to_file(report_file, @sprintf("\n\n**** Finished Fitting, evaluating now **** %s ---\n\n", Dates.format(now(), "e, dd u yyyy HH:MM:SS")))

    # evaluate the result with many trials, for accuracy
    cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA = JJ(extra_pars[:testruns], extra_pars[:testruns]; verbose=false, make_dict(args, pars3, merge(merge(mypars, extra_pars)))...)

    append_to_file(report_file, @sprintf("\n\n ****** writing to file %s *******\n\n", loadname))
    f["cost"] = cost
    f["cost1s"] = cost1s
    f["cost2s"] = cost2s
    f["hP"] = hP 
    f["hA"] = hA
    f["dP"] = dP
    f["dA"] = dA
    f["hBP"] = hBP
    f["hBA"] = hBA
    save(loadname, f)
end



