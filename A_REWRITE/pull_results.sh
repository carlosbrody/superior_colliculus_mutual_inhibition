#! /bin/csh

# This script goes into a running gcloud VM, adds, commits, and pushes the latest farm results, and then comes
# back to the local machine and executes git pull to get all the results.

if ( $#argv < 2 ) then
    echo Sorry, need two arguments, the root of the .csv file to be pulled and the name of the gcloud machine
    echo For example:  ./pull_results.sh negCosts proanti003    will pull negCosts_proanti003.csv from VM proanti003
    echo "If a third argument is present, it is taken to be the directory where all this will take place; default is ."
    exit
endif

if ( $#argv >= 3 ) then
    set farmdir = $argv[3]
else
    set farmdir = "."
endif

set rootdir = superior_colliculus_mutual_inhibition/A_REWRITE
if ( $argv[2] == "proanti024" ) then
    set zone = "us-east1-d"
else
    set zone = "us-east1-c"
endif


gcloud compute ssh --zone=$zone --command="cd $rootdir/; git pull; git add $farmdir/$argv[1]_$argv[2].csv; git commit -m "'"'"The latest from $farmdir ::$argv[1]_$argv[2],csv"'"'"  ; git push origin master" carlosbrody@$argv[2]

git pull
