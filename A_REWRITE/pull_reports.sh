#!/bin/csh

echo ssh brody@spock.princeton.edu "\rm -rf Reportz.tz ; tar cvfz Reports.tz Reports"
ssh brody@spock.princeton.edu "\rm -rf Reportz.tz ; tar cvfz Reports.tz Reports"

echo "scp brody@spock.princeton.edu:Reports.tz ../../Reports.tz"
scp brody@spock.princeton.edu:Reports.tz ../../Reports.tz

echo "pushd ../.. ; tar xvfz Reports.tz ; popd"
pushd ../.. ; tar xvfz Reports.tz ; popd


# ----------

# echo gcloud compute ssh --command "\rm -rf Reportz.tz ; tar cvfz Reports.tz Reports" carlosbrody@$argv[1]
# gcloud compute ssh --command "\rm -rf Reportz.tz ; tar cvfz Reports.tz Reports" carlosbrody@$argv[1]
#
# echo "gcloud compute scp carlosbrody@${argv[1]}:Reports.tz ../../Reports.tz"
# gcloud compute scp carlosbrody@${argv[1]}:Reports.tz ../../Reports.tz
#
# echo "pushd ../.. ; tar xvfz Reports.tz ; popd"
# pushd ../.. ; tar xvfz Reports.tz ; popd
