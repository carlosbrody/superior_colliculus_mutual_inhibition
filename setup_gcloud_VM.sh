#! /bin/csh

# This script sets up the software on a brand new, but existing, gcloud compute VM

if ( ! ($#argv >= 1) ) then
    echo Sorry, need at least one argument, the name of the gcloud machine
    echo
    echo For example:  ./setup_gcloud_VM.sh proanti028  optional:my_google_username
    echo 
    echo If a second argument is present, it is taken as the google username to apply. Default is carlosbrody.
    exit
endif

if ( $#argv >= 2 ) then
    set my_google_username = $argv[2]
else
    set my_google_username = carlosbrody
endif
set my_instance        = $argv[1]
if ( $argv[1] == "proanti009" ) then
    set zone = "us-east1-b"
else
    set zone = "us-east1-c"
endif

echo "gcloud compute scp --zone $zone package_installs.sh julia_start.jl grow_ForwardDiffChunks.sh $my_google_username@$my_instance":
gcloud compute scp --zone $zone package_installs.sh julia_start.jl grow_ForwardDiffChunks.sh "$my_google_username@$my_instance":

echo "gcloud compute ssh --zone $zone --command ""./package_installs.sh"" $my_google_username@$my_instance"
gcloud compute ssh --zone $zone --command "./package_installs.sh" $my_google_username@$my_instance
