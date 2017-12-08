#! /bin/csh

# This script goes into a running gcloud VM, adds, commits, and pushes the latest farm results, and then comes
# back to the local machine and executes git pull to get all the results.

if ( $#argv < 2 ) then
    echo Sorry, need two arguments, the name of the gcloud machine and the letter of the farm being created there
    echo For example:  ./pull_results.sh proanti003 D 
    echo "If a third argument is present, it is taken to be the directory where all this will take place; default is FarmFields"
    exit
endif

if ( $#argv >= 3 ) then
    set farmdir = $argv[3]
else
    set farmdir = "FarmFields"
endif

set rootdir = superior_colliculus_mutual_inhibition
if ( $argv[1] == "proanti024" ) then
    set zone = "us-east1-d"
else
    set zone = "us-east1-c"
endif


gcloud compute ssh --zone=$zone --command="cd $rootdir/; git pull; git add $farmdir/farm_$argv[2]_*; git commit -m "'"'"The latest from $farmdir ::$argv[2]"'"'"  ; git push origin master" carlosbrody@$argv[1] 

git pull
 
