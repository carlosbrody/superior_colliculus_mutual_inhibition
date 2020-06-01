#!/bin/csh

echo gcloud compute ssh --command "\rm -rf Reportz.tz ; tar cvfz Reports.tz Reports" carlosbrody@$argv[1]
gcloud compute ssh --command "\rm -rf Reportz.tz ; tar cvfz Reports.tz Reports" carlosbrody@$argv[1]

echo "gcloud compute scp carlosbrody@${argv[1]}:Reports.tz ../../Reports.tz"
gcloud compute scp carlosbrody@${argv[1]}:Reports.tz ../../Reports.tz

echo "pushd ../.. ; tar xvfz Reports.tz ; popd"
pushd ../.. ; tar xvfz Reports.tz ; popd
