### 2020-05-20 : old 4-node circuit was on C30 settings, not C32

Unbelievable. Was running on C30 settings, not C32. The difference is critical "Farm C32: Just like C30, but with variable delay periods, variable target periods, and no post-target periods".  The results for the 4-node circuit were not matching previous results, this was a possible reason. Now rerunning.

Output goes into `neg50Costs_$hostname.csv` after the first pass with 50 trials, and then passes with 1600 trials go into `neg1600Costs_hostname.csv`, set to do several passes of 50 iters each, replacing output as it goes, and stopping if the cost falls below -0.00028.  Won't do more than 150 iters on the 1600 trials.

Old C30 results are in C30/ directory, and backups of the neg1600 files are in BACKUPS/ directory

```julia
# Enough trials, iters for a real run:
extra_pars[:few_trials]                = 50       # number of trials to use in first pass
extra_pars[:firstPassNIter]            = 500      # maximum iterations in first pass
extra_pars[:many_trials]               = 1600     # of trials to use in further pass
extra_pars[:secondPassNIter]           = 50       # maximum iterations in further pass
extra_pars[:first_pass_cost_threshold] = -0.0002  # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:stoppingCostThreshold]     = -0.00028  # if below this cost, stop the minimization
extra_pars[:nFurtherPasses]            = 2        # after one further pass at many_trials and secondPassNIter, how many more of those to do before giving up
```

### 2020-05-19 : old 4-node circuit

#### multipass farming

On `proanti002`, running `softwallHessianFarming.jl` which is currently set to do a first pass of 500 iterations for 50 trials, and then another pass of 40 iterations with 1600.  

The script first includes `commonSetup.jl`, which in turn contains

```julia
extra_pars[:few_trials]                = 50       # number of trials to use in first pass
extra_pars[:firstPassNIter]            = 500      # maximum iterations in first pass
extra_pars[:many_trials]               = 1600     # of trials to use in second pass
extra_pars[:secondPassNIter]           = 40       # maximum iterations in second pass
extra_pars[:first_pass_cost_threshold] = -0.0002  # maximum cost threshold for a first pass run to seed a second pass run
```

Output goes into `neg50Cost_proanti002.csv` (successful first pass) and `neg1600Cost_proanti002.csv` (after second pass)

**Note:** readCostsFile will need upgrading from a file where lines alternated in number of items versus a file where every line has the same number of items.

#### refining previous farming

On `proanti003`, running `refineFistPassFarm.jl`, with 1600 trials for 40 iters (hardcoded), which picks up seeds from a previous run, in file `negCost_proanti003.csv`. Output will go into `negRefinedCost_proanti003.csv`. Expecting about 60 solutions. Between that and multipass farming we should be good. 

### 2020-05-17 : with Hessian is much better

`Optim.jl` has a bounded box constrained fit, but it doesn't take Hessians. Compared 50-trials fits, those with soft tanh walls and `NewtonTrustRegion()` finished much faster that the bounding box with only gradient. Using that from now on.
