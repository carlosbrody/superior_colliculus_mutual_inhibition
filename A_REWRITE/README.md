### 2020-06-14. Stopping 6-node run -- enough results in for now.

Have now stopped the last VMs running, `proanti006`, `proanti007`, `proanti008`, `proanti009`, `proanti010`, `proanti011`. Pulling reports from them finds 66 solutions with a cost less than -0.0001. Run the cell at the top of `sandbox.jl`. Not looking at final Farms outputs, but at Reports as it goes with 1000 trials. 

### 2020-06-03. Stopping old 4-node run to save money.

At this point it seems clear that the old 4-node code rerpoduces the old results, even in Julia 1.0.5, even if it hasn't produced final farms. Reading from the Reports files, filtering for 2nd-pass training, and sorting by cost on 1000 trials, we get the following out of the first cell in `sandbox.jl`:
```julia
  "cost"        "hW_P"     "vW_PA"    "dW_PA" 
 -4.43512e-5   -0.31714   -2.77002   -0.326953 
 -4.51673e-5    0.821592  -2.37589    1.79628  
 -4.8334e-5    -0.520949  -2.03382   -0.0968841
 -5.73134e-5    1.12911   -1.63477    0.799051 
 -5.95355e-5    0.443505  -1.52204    1.44696  
 -6.45954e-5    0.298148  -1.455      0.796008 
 -8.82276e-5    0.24758   -1.36222    0.671435 
 -0.000107178  -0.424304  -1.62051    0.455756 
 -0.00010807    0.311738  -2.18619    1.47656  
 -0.000116492   1.00027   -2.24554    0.0822045
 -0.000117643   0.454465  -1.23349    0.376089 
 -0.000141724  -1.03584   -1.14418    0.0712079
 -0.000161144   0.705797  -2.05363    0.117806 
 -0.000162982  -0.75624   -2.08015   -0.0950036
 -0.000178077  -0.893418  -2.20795    0.300835 
 -0.000216462  -0.794968  -1.54534   -0.51226  
 ```
 One can clearly see that `hW_P` has mixed signs, `vW_PA` is always negative, and `vW_PA < dW_PA`, even when `dW_PA` is negative.  I.e., confirming old results.
 
Am therefore stopping `proanti002`, `proanti004`, and `proanti005` to save money. No need to keep them going. They can be restarted if need be.

### 2020-06-01. Starting 6-node on old optimization code

None of the runs on `proanti006` through `proanti011` were going past pass 5, and only two made it past pass 4. The strategy has clearly failed. The old code would train to optimum on 50 trials (with maxIter=1000); there do an evaluation, and if the model passed the evaluation, train to optimum on 1000 trials. The evaluation has now been put into `evaluateModel()` in the `ProAnti.jl` module. 

We now have the question of whether new optimization procedures, using `Optim.jl` would run that old strategy faster than the old optimization procedures. Not answering that yet. For now, am taking what seems like a more guaranteed path. We are simply starting a 6-node run on 6 VMs (`proanti006` through `proanti011`), using the old optimization code. We are using files `sixNode_reduced_farm_C32.jl` and `sixNodeSetup_C32.jl`, which are almost identical to `opto_reduced_farm_C32.jl` and `setup_C32.jl`. Report files start with `r6` and farms go into `../../Farms_6N`.

6-node code has 26 parameters, old code had 16; 26^2/16^2 ~= 2.65. In practice, 6-node code looks like it runs about 3 times slower. **So we should expect results in about 3 weeks**. :(

#### old code is slow, but is successfully reproducing old results, now in Julia 1.0

A straight re-run of the old code (merely upgraded to Julia 1.0) seems to be reproducing the old results -- as it should!! Running on `proanti002`, `proanti004`, and `proanti005` after 7 days it has produced no farms yet, but looking at report files shows they are on their way. 158 CPUs are running, and at this point 63 of them have gone into pass 2 further training. 

Of those 63, a few have gone beyond the -0.0001 threshold and display the expected results (a) `hW_P` is varied in sign; (b) `vW_PA` is almost always negative; (c) `dW_PA` is usually positive, and when it is negative, it is still greater than `vW_PA`. Fron `sandbox.jl`:
```julia
 "cost"        "hW_P"     "vW_PA"    "dW_PA" 
 -3.66935e-5    0.84933   -1.50096    0.815532 
 -4.06454e-5   -0.298154  -2.66791   -0.331266 
 -4.247e-5     -0.507152  -2.00155   -0.0776131
 -8.65474e-5    0.244958  -1.29238    0.602195 
 -9.73053e-5   -0.423981  -1.6231     0.457797 
 -0.000101394   0.330891  -2.17932    1.49553  
 -0.000109802   0.992769  -2.21862    0.0740309
 -0.000119361  -0.871576  -2.20647    0.315527 
 -0.000121367  -1.02623   -1.13486    0.0724739
 -0.000151369   0.692988  -2.03862    0.104209 
 -0.00017751   -0.728513  -1.54116   -0.603525 
 ```


### 2020-05-30 -- Noise vs trials numbers, Higher Iters, Recruiting Spock

#### Costs change with number of trials used to assess them (see [Costs_vs_ntrials.ipynb](Costs_vs_ntrials.ipynb) )

You might think that different  numbers of trials just affect variability in the cost function, not the mean value.  **This is wrong.**

Take a well-trained network and a random seed equal to the one used to measure cost at, say 10000 trials.  The cost function includes a squared error term for how close to the target the mean hits are. At few trials, the typical distance to the target will be greater than for many trials. So, as number of trials grows, cost falls, even under these conditions.  

Conversely, for a badly trained network, we expect that evaluating it with few trials used to train it is equivalent to picking out that noise instance that gives you the lowest cost-- in this case, increasing number of trials will lead to increased cost, because of overfitting.

In sum, this means that our strategy of using cost thresholds to switch between passes may be flawed.

<img src=costs_vs_ntrials.jpg>

#### 2020-05-30 5pm Recruiting Spock

Now running the old `opto_reduced_farmC32.jl` code, as in `proanti002`, on Spock. Takes about 4min30 per 10 iterations instead of 6 min on the Google VMs. A little faster but not hugely so.  Run the following four times while logged into spock.
```
sbatch --array=0-43 spockFarm.sh opto_reduced_farmC32.jl 
```
This will use all 44 cores on each CPU. Between each of those 4 runs, I waited until all runs had been assigned a core (so filenames in the Reports directory don't get confusing), by running this
```
squeue | grep spockFar | grep None
```
and waiting until it came back empty.

### 2020-05-30  Running with much higher number of maxIter

It seems that it equally often gets stuck because of low gradients, in which case it properly aborts, or because it ran out of iterations. Maybe it's best to let it have lots of iterations, trust that if it is really off on the wrong path it'll get to a local minimum and abort.

To look at this, starting `proanti009`, `proanti010` and `proanti011` with `softwallHessianFarming.jl` and 
```julia
extra_pars[:nPasses]                   =  8       # of pass blocks below

extra_pars[:pass1NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass1NIter]                = 4000     # maximum iterations in first pass
extra_pars[:pass1CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass1RnD]                  = [1.2]    # rule and delay period range

extra_pars[:pass2NTrials]              = 25        # number of trials to use in first pass
extra_pars[:pass2NIter]                = 4000      # maximum iterations in first pass
extra_pars[:pass2CostThreshold]        = 0         # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass2RnD]                  = [1.15  1.2]    # rule and delay period range

extra_pars[:pass3NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass3NIter]                = 4000      # maximum iterations in first pass
extra_pars[:pass3CostThreshold]        = 0         # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass3RnD]                  = [1.1   1.2]    # rule and delay period range

extra_pars[:pass4NTrials]              = 25        # number of trials to use in first pass
extra_pars[:pass4NIter]                = 4000      # maximum iterations in first pass
extra_pars[:pass4CostThreshold]        = 0         # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass4RnD]                  = [1.05   1.2]    # rule and delay period range

extra_pars[:pass5NTrials]              = 25        # number of trials to use in first pass
extra_pars[:pass5NIter]                = 4000      # maximum iterations in first pass
extra_pars[:pass5CostThreshold]        = -0.0001   # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass5RnD]                  = [1.0   1.2]    # rule and delay period range

extra_pars[:pass6NTrials]              = 50        # number of trials to use in first pass
extra_pars[:pass6NIter]                = 4000      # maximum iterations in first pass
extra_pars[:pass6CostThreshold]        = -0.0001   # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass6RnD]                  = [1.0   1.2]    # rule and delay period range

extra_pars[:pass7NTrials]              = 100        # number of trials to use in first pass
extra_pars[:pass7NIter]                = 4000       # maximum iterations in first pass
extra_pars[:pass7CostThreshold]        = -0.00015   # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass7RnD]                  = [1.0   1.2]    # rule and delay period range

extra_pars[:pass8NTrials]              = 4000       # number of trials to use in first pass
extra_pars[:pass8NIter]                = 400        # maximum iterations in first pass
extra_pars[:pass8CostThreshold]        = -0.00028   # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass8RnD]                  = [1.0   1.2]    # rule and delay period range
```


### 2020-05-26  Was using the wrong cost limit!  Previous "good solutions" were at -0.0001

#### 2020-05-26 : 3.21pm -- Stopping `proanti001` and `proanti003` to save money

Both of these machines have only a single process that is worthwhile. On `proanti001` it can be found in `julia_outs_proanti001_13` and on `proanti003` it is `julia_outs_proanti003_38`. Only the latter has a trace of the parameters values, so only the latter can be recovered. To recover, start the VM first.


### 2020-05-24

#### 10am

Not clear old code is faster. Starting up `proanti006`, `proanti007`, and `proanti008` on `softwallHessianFarming.jl` (with Optim package), on the grow-backwards settings from `proanti002` (see 2020-05-23 7:30am) :
```julia
# Enough trials, iters for a real run:
extra_pars[:nPasses]                   =  8       # of pass blocks below

extra_pars[:pass1NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass1NIter]                = 400      # maximum iterations in first pass
extra_pars[:pass1CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass1RnD]                  = [1.2]    # rule and delay period range

extra_pars[:pass2NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass2NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass2CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass2RnD]                  = [1.15  1.2]    # rule and delay period range

extra_pars[:pass3NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass3NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass3CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass3RnD]                  = [1.1   1.2]    # rule and delay period range

extra_pars[:pass4NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass4NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass4CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass4RnD]                  = [1.05   1.2]    # rule and delay period range

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

#### 2am

When we run only one process, it runs at 3 min/10 iters. 32 processes, 6 min/10 iters.  48 processes, 9 min/10 iters. 56 processes, 11 min/10 iters.  Increasing memory does not appear to make a difference.

### 2020-05-23

#### 11:59pm

Approach has produced only one viable candidate on `proanti002` and `proanti003`. Clearly not workable. Dropping Optim.jl, upgraded old optimization code to Julia 1.0.5, and running `setup_C32.jl` / `opto_reduced_farm_C32.jl`, old code, on `proanti002`, `proanti004`, and `proanti005`.

`proanti001` and `proanti005` have one viable candidate each (runs 13 and 38, respectively), and am letting those single ones play out.

#### 4pm

Pure gradient on `proanti004` definitely not cutting it. Closing that down, turning off VM.

#### 7:30am

Pure gradient seems extremely slow. `proanti002` and `proanti003`, with Hessian and gradual growth, also slow but kind of coming along. Switched `proanti002` to not expand rule-and-delay variability centered, but from long towards shorter:
```julia
extra_pars[:nPasses]                   =  8       # of pass blocks below

extra_pars[:pass1NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass1NIter]                = 400      # maximum iterations in first pass
extra_pars[:pass1CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass1RnD]                  = [1.2]    # rule and delay period range

extra_pars[:pass2NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass2NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass2CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass2RnD]                  = [1.15  1.2]    # rule and delay period range

extra_pars[:pass3NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass3NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass3CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass3RnD]                  = [1.1   1.2]    # rule and delay period range

extra_pars[:pass4NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass4NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass4CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass4RnD]                  = [1.05   1.2]    # rule and delay period range

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

#### 1:00am

Running a pure gradient, no Hessian attempt on `proanti004`.

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
