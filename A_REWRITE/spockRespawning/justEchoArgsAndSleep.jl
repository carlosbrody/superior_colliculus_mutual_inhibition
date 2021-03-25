using ArgParse

s = ArgParseSettings()
@add_arg_table! s begin
    "--logfileID"
        help = "the string representing where the output will be logged"
        required = true
    "--bspockID", "-b"
        help = "the job ID that the respawner is keeping track of for iterations of this process"
        required = true
    "--respawn", "-r"
        help = "we should respawn, using the logfile for initial conditions. "*
            "Otherwise start up from scratch."
        action = :store_true
    "my_run_number"
        help = "if not on spock, will be first positional argument"
        arg_type = Int64
    "tot_n_runs"
        help = "if not on spock, will be second positional argument"
        arg_type = Int64
    "--hostname"
        help = "string used within farm and report filenames"
        arg_type = String
        default = String(chomp(read(`hostname`, String)))
end

parsed_args = parse_args(s)

println(parsed_args)

for k in keys(parsed_args)
    eval(Meta.parse("$k = parsed_args[\"$k\"]"))
end
