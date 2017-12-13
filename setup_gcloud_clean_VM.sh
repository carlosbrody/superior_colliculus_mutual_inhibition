#! /bin/csh

# This script erases everything on a gcloud compute VM

if ( ! ($#argv >= 1) ) then
    echo Sorry, need at least one argument, the name of the gcloud machine
    echo
    echo For example:  ./setup_gcloud_clean_VM.sh proanti028  optional:my_google_username
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
set zone               = "us-east1-c"


echo "gcloud compute ssh --zone $zone --command ""sudo \rm -rf * ; sudo \rm -rf .julia ; sudo \rm -rf ../marinopagan/* ; sudo \rm -rf ../alexpiet/* ; sudo \rm -rf ../alex.piet/"" $my_google_username@$my_instance"

gcloud compute ssh --zone $zone --command "sudo \rm -rf * ; sudo \rm -rf .julia ; sudo \rm -rf ../marinopagan/* ; sudo \rm -rf ../alexpiet/* ; sudo \rm -rf ../alex.piet/*" $my_google_username@$my_instance
