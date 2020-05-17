if length(ARGS)>0  &&  ~isnull(tryparse(Int64, ARGS[1])); my_run_number = parse(Int64, ARGS[1]);
else                                                      my_run_number = 1; # I am process my_run_number
end
if length(ARGS)>0  &&  ~isnull(tryparse(Int64, ARGS[2])); tot_n_runs    = parse(Int64, ARGS[2]);
else                                                      tot_n_runs = 1;   # I'm being run as part of tot_n_run processes

println("I am $my_run_number out of $tot_n_runs")
