This directory contains 65 solutions to the 6-node network problem.

File `solutions6.mat` has:

  *  "costs" => vector of final costs for each solution
  * "fnames"=> corresponding vector of filenames
  * "paramvals"=> solutions-by-nparams matrix of parameter values
  * "argnames"=> nparams vector of strings, indicating what each parameter value represents (i.e., these are the names of each of the columns of paramvals)
  * "randseeds"=> vector of Float64s, the random seeds with which each solution was run
  * "weightMatrixMap"=> 6-by-6 cell of strings, with each entry indicating the name (in argnames) of the weight in that entry
   
