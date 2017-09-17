#! /bin/csh

# This script goes into a running gcloud VM, adds, commits, and pushes the latest farm results, and then comes
# back to the local machine and executes git pull to get all the results.

if ( $#argv != 2 ) then
    echo Sorry, need two arguments, the name of the gcloud machine and the letter of the farm being created there
    echo For example:  ./pull_results.sh proanti003 D
    exit
endif

gcloud compute ssh --zone="us-east1-c" --command="cd superior_colliculus_mutual_inhibition/; git pull; git add FarmFields/farm_$argv[2]_*; git commit -m "'"'"The latest from farm $argv[2]"'"'"  ; git push origin master" carlosbrody@$argv[1] 

git pull
 
