POTENTIAL ISSUES
    Julia 0.6.3? Does it matter vs 0.6.0?

SETUP NEW USER
    (1) Checking out the GITHUB Repository
    >> git clone "http://alexpiet@github.com/carlosbrody/superior_colliculus_mutual_inhibition.git"
    >> # will be prompted for your github password

    (2) Adding Julia packages
    >> sbatch ./setup.sh
    If this script crashes, open julia and manually load.
    You only need to do this once
    >> module load julia/0.6.3
    >> julia
    >> Pkg.add("MAT") #load each package manually

    (3) Learn how to submit jobs
    Example script              >> sbatch ./example.sh
        This demonstrates the basics of submitting a job
    Example job array           >> sbatch --array=0-9 ./array_example.sh "C1"
        This will start ten farm runs of "C1", each identified 0-9.
        Each saves a file in directory C1/C1_0.mat
        Each farm run gets a unique job ID, as well as their number 0-9.
        These two numbers can be used to seed a random number generator uniquely and simultaneously

    While jobs are running use 'ls -al' to look at name of logfile
    Then use 'tail -f <logfilename>' to see updates while it works, lags behind realtime

    (4) Grow Forward Diff Chunks    >> ./grow_ForwardDiffChunks.sh


