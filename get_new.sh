#! /bin/csh

# This script goes into a directory in running gcloud VM, finds any new files there NOT also in a corresponding local
# directory, and copies the new ones to the local directory.


if ( $#argv != 2 ) then
    echo Sorry, need two arguments, the name of the gcloud machine and the directory to look in. 
    echo For example:  ./get_new.sh proanti025 ../NewFarms
    exit
endif

set rootdir = superior_colliculus_mutual_inhibition
if ( $argv[1] == "proanti024" ) then
    set zone = "us-east1-d"
else
    set zone = "us-east1-c"
endif

# Get a local list of all the remote files
echo "gcloud compute ssh  --zone=$zone --command=cd $rootdir/; ls $argv[2] carlosbrody@$argv[1] | sort | uniq > temp_remote.txt"
gcloud compute ssh  --zone=$zone --command="cd $rootdir/; ls $argv[2]" carlosbrody@$argv[1] | sort | uniq > temp_remote.txt

# Make local list of local files
echo "ls $argv[2] | sort | uniq > temp_local.txt"
ls $argv[2] | sort | uniq > temp_local.txt

# Make local list of all files in remote that aren't in local. These are the files we want to copy
echo "comm -2 -3 temp_remote.txt temp_local.txt > temp_newfiles.txt"
comm -2 -3 temp_remote.txt temp_local.txt > temp_newfiles.txt

# Copy the wanted list to the desired directory in remote machine 
# echo "gcloud compute scp  --zone=$zone temp_newfiles.txt carlosbrody@$argv[1]: $rootdir/$argv[2]"
gcloud compute scp --zone=$zone temp_newfiles.txt "carlosbrody@$argv[1]":$rootdir/$argv[2]

# ssh into remote machine and tar up the files  then clean up by deleting the list after having used it
echo "gcloud compute ssh  --zone=$zone --command=cd $rootdir/; cd $argv[2] ; tar cvfz temp_tarball.tz -T temp_newfiles.txt ; \rm temp_newfiles.txt carlosbrody@$argv[1]"
gcloud compute ssh  --zone=$zone --command="cd $rootdir/; cd $argv[2] ; tar cvfz temp_tarball.tz -T temp_newfiles.txt ; \rm temp_newfiles.txt" carlosbrody@$argv[1]

# Copy the tarball to the local machine
echo "gcloud compute scp  --zone=$zone carlosbrody@$argv[1] :$rootdir/$argv[2]/temp_tarball.tz $argv[2]"
gcloud compute scp  --zone=$zone "carlosbrody@$argv[1]":$rootdir/$argv[2]/temp_tarball.tz $argv[2]

# Change to desired dir, unpack then remove the tarball, then return to regular dir
echo "pushd $argv[2] ; tar xvfz temp_tarball.tz ; \rm temp_tarball.tz ; popd"
pushd $argv[2] ; tar xvfz temp_tarball.tz ; \rm temp_tarball.tz ; popd

# Clean up locally
echo "\rm temp_remote.txt temp_local.txt temp_newfiles.txt"
\rm temp_remote.txt temp_local.txt temp_newfiles.txt

# Clean up remotely
echo "gcloud compute ssh  --zone=$zone --command cd $rootdir/$argv[2] ; \rm temp_tarball.tz carlosbrody@$argv[1]"
gcloud compute ssh  --zone=$zone --command "cd $rootdir/$argv[2] ; \rm temp_tarball.tz" carlosbrody@$argv[1]
