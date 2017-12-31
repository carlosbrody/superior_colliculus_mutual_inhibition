#! /bin/csh

# This script goes into a directory in running gcloud VM, finds any new files there NOT also in a corresponding local
# directory, and copies the new ones to the local directory.


if ( ! ($#argv >= 2) ) then
    echo Sorry, need at least three arguments, the name of the gcloud machine, the target file 
    echo pattern and the directory to look in. 
    echo
    echo For example:  ./get_new.sh proanti025 _C17_ ../NewFarms  optional: local_dir
    echo ""
    echo if a fourth argument is provided, then the third is interpreted as the cloud dir, 
    echo and the fourth as local dir.If no fourth argument is provided, then the third is
    echo interpreted as name of both local and remote dir
    exit
endif

set rootdir = superior_colliculus_mutual_inhibition
if ( $argv[1] == "proanti024" ) then
    set zone = "us-east1-d"
else
    if ( $argv[1] == "proanti009" ) then
	set zone = "us-east1-b"
    else
	set zone = "us-east1-c"
    endif
endif

set patt   = $argv[2]
set remdir = $argv[3]

if ( $#argv >= 4 ) then
    set localdir = $argv[4]
else
    set localdir = $argv[3]
endif


# Get a local list of all the remote files
echo "gcloud compute ssh  --zone=$zone --command=cd $rootdir/$remdir ; ls *$patt* carlosbrody@$argv[1] | sort | uniq > temp_remote.txt"
gcloud compute ssh  --zone=$zone --command="cd $rootdir/$remdir ; ls *$patt*" carlosbrody@$argv[1] | sort | uniq > temp_remote.txt


# Make local list of local files
echo "pushd $localdir ; ls *$patt* | sort | uniq > temp_local.txt ; popd"
pushd $localdir ; ls *$patt* | sort | uniq > temp_local.txt ; popd

# Make local list of all files in remote that aren't in local. These are the files we want to copy
echo "comm -2 -3 temp_remote.txt temp_local.txt > temp_newfiles.txt"
comm -2 -3 temp_remote.txt $localdir/temp_local.txt > temp_newfiles.txt

# Copy the wanted list to the desired directory in remote machine 
# echo "gcloud compute scp  --zone=$zone temp_newfiles.txt carlosbrody@$argv[1]: $rootdir/$remdir"
gcloud compute scp --zone=$zone temp_newfiles.txt "carlosbrody@$argv[1]":$rootdir/$remdir

# ssh into remote machine and tar up the files  then clean up by deleting the list after having used it
echo "gcloud compute ssh  --zone=$zone --command=cd $rootdir; cd $remdir ; tar cvfz temp_tarball.tz -T temp_newfiles.txt ; \rm temp_newfiles.txt carlosbrody@$argv[1]"
gcloud compute ssh  --zone=$zone --command="cd $rootdir; cd $remdir ; tar cvfz temp_tarball.tz -T temp_newfiles.txt ; \rm temp_newfiles.txt" carlosbrody@$argv[1]

# Copy the tarball to the local machine
echo "gcloud compute scp  --zone=$zone carlosbrody@$argv[1] :$rootdir/$remdir/temp_tarball.tz $localdir"
gcloud compute scp  --zone=$zone "carlosbrody@$argv[1]":$rootdir/$remdir/temp_tarball.tz $localdir

# Change to desired dir, unpack then remove the tarball, then return to regular dir
echo "pushd $localdir ; tar xvfz temp_tarball.tz ; \rm temp_tarball.tz ; popd"
pushd $localdir ; tar xvfz temp_tarball.tz ; \rm temp_tarball.tz ; popd

# Clean up locally
echo "\rm temp_remote.txt $localdir/temp_local.txt temp_newfiles.txt"
\rm temp_remote.txt $localdir/temp_local.txt temp_newfiles.txt

# Clean up remotely
echo "gcloud compute ssh  --zone=$zone --command cd $rootdir/$remdir ; \rm temp_tarball.tz carlosbrody@$argv[1]"
gcloud compute ssh  --zone=$zone --command "cd $rootdir/$remdir ; \rm temp_tarball.tz" carlosbrody@$argv[1]
