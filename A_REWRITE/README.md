### 2020-05-23

#### 11:58pm

Latest approach seemed an improvement, will now turn over both `proanti002` and `proanti003` to it.
```julia
extra_pars[:nPasses]                   =  8       # of pass blocks below

extra_pars[:pass1NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass1NIter]                = 400      # maximum iterations in first pass
extra_pars[:pass1CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass1RnD]                  = [1.1]    # rule and delay period range

extra_pars[:pass2NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass2NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass2CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass2RnD]                  = [1.075  1.125]    # rule and delay period range

extra_pars[:pass3NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass3NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass3CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass3RnD]                  = [1.05   1.15]    # rule and delay period range

extra_pars[:pass4NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass4NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass4CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass4RnD]                  = [1.025   1.175]    # rule and delay period range

extra_pars[:pass5NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass5NIter]                = 400      # maximum iterations in first pass
extra_pars[:pass5CostThreshold]        = -0.0001  # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass5RnD]                  = [1.0   1.2]    # rule and delay period range

extra_pars[:pass6NTrials]              = 100      # number of trials to use in first pass
extra_pars[:pass6NIter]                = 400      # maximum iterations in first pass
extra_pars[:pass6CostThreshold]        = -0.0001  # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass6RnD]                  = [1.0   1.2]    # rule and delay period range

extra_pars[:pass7NTrials]              = 400       # number of trials to use in first pass
extra_pars[:pass7NIter]                = 400       # maximum iterations in first pass
extra_pars[:pass7CostThreshold]        = -0.00015  # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass7RnD]                  = [1.0   1.2]    # rule and delay period range

extra_pars[:pass8NTrials]              = 1600       # number of trials to use in first pass
extra_pars[:pass8NIter]                = 400        # maximum iterations in first pass
extra_pars[:pass8CostThreshold]        = -0.00028   # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass8RnD]                  = [1.0   1.2]    # rule and delay period range
```

#### 7:00pm

The initial [1.075 1.125] rule and delay period variability still seemed too much. Ran a trial on `proanti003`, doing an initial *no* variability pass [1.1], and then adding variability. 

#### 4:24pm

Stopping `proanti005`, will focus on stepwise broadening of rule-and-delay variability. First pass, on [1.075 1.125], will have 400 iters, it doesn't make it to 0 otherwise. After killing processes on `proanti002`, `proanti003`, and `proanti005`, `proanti005` will be stopped (save some money) and both of `proanti002` and `proanti003` will run the following, with the focus at first on whether anything gets past the first pass:
```julia
extra_pars[:pass1NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass1NIter]                = 400      # maximum iterations in first pass
extra_pars[:pass1CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass1RnD]                  = [1.075 1.125]    # rule and delay period range

extra_pars[:pass2NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass2NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass2CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass2RnD]                  = [1.05  1.1]    # rule and delay period range

extra_pars[:pass3NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass3NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass3CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass3RnD]                  = [1.0   1.2]    # rule and delay period range
```

#### 2:52pm

`proanti003` restarted, variability in delay rule now halved to
```julia
:rule_and_delay_period  =>          1.1,
:rule_and_delay_periods =>          [1.075 1.125],
```
and old results saved in `failed003.tz`

#### Current status:

* `proanti001`
First word count is number of runs that passed first pass; second word count is number of attempted first passes.
```julia
date; gcloud compute ssh --command "cd A_REWRITE; grep BELOW ju* | wc ; grep loop ju* | wc" carlosbrody@proanti001
Sat May 23 14:22:55 EDT 2020
      1       7      74
    275    1916   18331

```
* `proanti002`
```julia
date; gcloud compute ssh --command "cd A_REWRITE; grep BELOW ju* | wc ; grep loop ju* | wc" carlosbrody@proanti002
Sat May 23 14:23:13 EDT 2020
      4      28     295
    603    4221   40314
```
* `proanti003`
```julia
date; gcloud compute ssh --command "cd A_REWRITE; grep BELOW ju* | wc ; grep loop ju* | wc" carlosbrody@proanti003
Sat May 23 14:23:55 EDT 2020
      0       0       0
    271    1897   18116
```
* `proanti005`
```julia
date; gcloud compute ssh --command "cd A_REWRITE; grep BELOW ju* | wc ; grep loop ju* | wc" carlosbrody@proanti005
Sat May 23 14:24:32 EDT 2020
      5      35     370
    852    5964   56955
```


#### Current summary close to C32: [githash](https://github.com/carlosbrody/superior_colliculus_mutual_inhibition/commit/01481cda6b20902f8a3a5f929d4c0326f9099c8e)

- `proanti003` to see whether starting with smaller variability can let us get to a good regime in full variability, is on 200iters/25trials/0threshold/sigma[-2 2] first pass has smaller variability rule and delay:
```julia
:rule_and_delay_period  =>          1.1,
:rule_and_delay_periods =>          [1.05 1.15],
```
second pass changes threshold, keeps trials and iters: 200iters/25trials/-0.0001threshold/sigma[-2 2] and changes to full C32 variability:
```julia
   mypars[:rule_and_delay_periods] = [1.0 1.2]
```

### 2020-05-22

#### Current summary C32 settings:

- `proanti001` on 1500iters/50trials/-0.0001threshold/sigma[0 2] first pass, out of 257 attempts 1 made it to 2nd pass.
- `proanti002` on 200iters/25trials/0threshold/sigma[-2 2] first pass, out of 342 attempts 1 made it to 2nd pass.
- `proanti005` on 100iters/50trials/0threshold/sigma[0 2] first pass, out of 576 attempts 3 made it to 2nd pass.

All of these use full C32 settings
```julia
:rule_and_delay_period  =>          1.2,
:rule_and_delay_periods =>          [1.0 1.2],
:target_period          =>          0.6,
:target_periods         =>          [0.45 0.6],
:post_target_period     =>          0,
:post_target_periods    =>          [0],
```




### 2020-05-21 : C32 settings make finding solutions much harder

#### Update 2

Using `proanti002` to step from C30 settings to C32 settings. Everything up to but not including variable delay periods seems to work fine in terms of finding solutions within a reasonable time. Using only 25 trials in first pass for speed.
```julia
:rule_and_delay_period  =>          1.2,
:rule_and_delay_periods =>          [1.2],
:target_period          =>          0.6,
:target_periods         =>          [0.45 0.6],
:post_target_period     =>          0,
:post_target_periods    =>          [0],

extra_pars[:few_trials]                = 25       # number of trials to use in first pass
extra_pars[:firstPassNIter]            = 200      # maximum iterations in first pass
extra_pars[:many_trials]               = 400     # of trials to use in further pass
extra_pars[:secondPassNIter]           = 200       # maximum iterations in further pass
extra_pars[:first_pass_cost_threshold] = 0         # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:stoppingCostThreshold]     = -0.00028  # if below this cost, stop the minimization
```

#### Update 1

After 282 attempts on `proanti003`, not a single one went into second pass.  Also, on `proanti005`, where we now look at parameter values as they evolve, some suggestion that we might be hitting a boundary at `sigma=0`. Old C32 settings had bounds on `sigma` of `[-2 2]`. Killing the useless `proanti003` attempts, and now using that VM to try the `[-2 2]` bounds, with only 25 trials, and 200 nIters on first pass so we can see evolution better, on `proanti003`.
```julia
extra_pars[:few_trials]                = 25       # number of trials to use in first pass
extra_pars[:firstPassNIter]            = 200      # maximum iterations in first pass
extra_pars[:many_trials]               = 1600     # of trials to use in further pass
extra_pars[:secondPassNIter]           = 200       # maximum iterations in further pass
extra_pars[:first_pass_cost_threshold] = 0         # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:stoppingCostThreshold]     = -0.00028  # if below this cost, stop the minimization
```

#### Morning

Ran on `proanti002` and `proanti003` doing 500 first pass iterations and 150 second pass; on `proanti001` and `proanti005` doing 1500 first pass iters; found only one potential solution, still running (`proanti001:julia_outs_proanti001_13`). In that one, cost had gone negative by 100th initial iter. Another potential solution on `proanti005` did not get to negative until more than 1000 first pass iters.

Now running with 250 first pass maximum iters (`proanti002`, cost threshold for second pass -0.0001) and only 100 first pass maximum iters (`proanti003`, cost threshold 0), trying to quickly eliminate fruitless attempts on first pass. Still using Hessian and NewtonTrustRegion(). Now using Optim's callback function, and using `constantFarm.sh` to start processes; this latter script runs in the background, respawning processes when they are missing. Using this because there is some bug that causes OpenBLAS to crash if to many calls to `Optim.optimize` are made. (Each one can run for a long time, it seems to be the number of outer calls.)

Output goes to `neg50Costs_$hostname` (successful first pass) and `neg1600Costs_$hostname` (result of second pass).

**TO-DO:** Considering 
- (a) look at minimizer values as they evolve, for a diagnostic (are we hitting walls?); *Done: running on proanti005* with params:
```julia
extra_pars[:few_trials]                = 50       # number of trials to use in first pass
extra_pars[:firstPassNIter]            = 100      # maximum iterations in first pass
extra_pars[:many_trials]               = 1600     # of trials to use in further pass
extra_pars[:secondPassNIter]           = 200       # maximum iterations in further pass
extra_pars[:first_pass_cost_threshold] = 0         # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:stoppingCostThreshold]     = -0.00028  # if below this cost, stop the minimization
```
- (b) doing a two-tiered stopping callback, e.g., must reach certain cost by 60 iters but continues with few trials even if that first threshold is reached, so as not to jump into slow 1600 trial version too soon.
- (c) do a sanity check on running on C30 settings to see whether those are still easy
- (d) stepwise go from C30 settings to C32 settings to see where the hard part comes in

### 2020-05-20 : old 4-node circuit was on C30 settings, not C32

#### extending nIter on proanti001 and proanti002

Looks like the C32 settings are much harder to train -- indeed, now I remember that!  So in addition to running as below on `proanti002` and `proanti003` VMs, have now revived `proanti001` and `proanti004`, and am running on them with the following parameters, bigger nIter and easier threshold from going from few trials to many trials.

```julia
# Enough trials, iters for a real run:
extra_pars[:few_trials]                = 50       # number of trials to use in first pass
extra_pars[:firstPassNIter]            = 1500      # maximum iterations in first pass
extra_pars[:many_trials]               = 1600     # of trials to use in further pass
extra_pars[:secondPassNIter]           = 150       # maximum iterations in further pass
extra_pars[:first_pass_cost_threshold] = -0.0001  # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:stoppingCostThreshold]     = -0.00028  # if below this cost, stop the minimization
extra_pars[:nFurtherPasses]            = 2        # after one further pass at many_trials and secondPassNIter, how many more of those to do before giving up
```

#### fits on proanti002 and proanti003

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
