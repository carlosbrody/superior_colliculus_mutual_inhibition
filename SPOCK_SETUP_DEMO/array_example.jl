# farm function

using Random

function array_example(farm_id, job_id, task_id)
    # print out some info
    println(farm_id)
    println(job_id)
    println(task_id)

    # set up random number
    random_seed = parse(Int64,job_id)*parse(Int64,task_id)
    Random.seed!(random_seed)
    random_num = rand(2)

    # save a data file with inputs and random number
    if !isdir(String(farm_id))
        mkdir(String(farm_id))
    end
    myfilename = String(farm_id)*"/"*string(farm_id)*"_"*string(task_id)*".mat"
    matwrite(myfilename, Dict("farm_id"=>farm_id, "job_id"=>job_id, "task_id"=>task_id,"random_seed"=>random_seed,"random_num"=>random_num))

    println("All done!")
end

# set up packages to be used
using MAT

# Call the function using the inputs from the bash script
array_example(ARGS[1], ARGS[2], ARGS[3])
