using ArgParse

s = ArgParseSettings()
@add_arg_table! s begin
    "--sleepFor"
       help = "how many seconds to sleep for, then end"
       arg_type = Int64
       default = 10
    "--bspockID", "-b"
        help = "the job ID that the spock respawner is keeping track of for iterations of this process"
        arg_type = String
    "--onSpock", "-s"
        help = "flag used to indicate we're running on spock"
        action = :store_true
    "--respawn", "-r"
        help = "we should respawn, using the logfile for initial conditions. "*
            "Otherwise start up from scratch."
        action = :store_true
    "--hostname"
        help = "string used within farm and report filenames"
        arg_type = String
        default = String(chomp(read(`hostname`, String)))
end

parsed_args = parse_args(s)

println(parsed_args)

sleep(parsed_args["sleepFor"])



