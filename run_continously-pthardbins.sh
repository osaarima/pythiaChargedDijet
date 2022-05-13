#!/bin/bash 
############################

PROG=`basename $0`
if [ $# -ne 3 ]
then
    echo "Usage: $PROG pythiaconfig.cmnd firstSeed name"
    exit 0;
fi

function setenv(){ export $1=$2; }

setenv pyconfigname $1
setenv firstSeed $2
setenv oname $3
setenv trackingIneff 0.0
setenv runs 1000

setenv Disk         `pwd`
setenv Out_DIR      $Disk/output/${oname}/data
setenv LOG_DIR      $Disk/output/${oname}/logs

mkdir -p $Out_DIR
mkdir -p $LOG_DIR

echo Running $oname $runs runs with first seed no $firstSeed.

for ((i=0;i<$runs; i++))
do
    currentSeed=$((firstSeed+i))
    echo Running: runSeed$currentSeed.root
    ./pythiaChargedDijet $pyconfigname 30 40 20 ${Out_DIR}/030_040runSeed$(($currentSeed)).root $currentSeed $trackingIneff >& ${LOG_DIR}/030_040runSeed$currentSeed.log
    ./pythiaChargedDijet $pyconfigname 40 65 20 ${Out_DIR}/040_065runSeed$(($currentSeed)).root $currentSeed $trackingIneff >& ${LOG_DIR}/040_065runSeed$currentSeed.log
    ./pythiaChargedDijet $pyconfigname 65 90 20 ${Out_DIR}/065_090runSeed$(($currentSeed)).root $currentSeed $trackingIneff >& ${LOG_DIR}/065_090runSeed$currentSeed.log
    ./pythiaChargedDijet $pyconfigname 90 120 20 ${Out_DIR}/090_120runSeed$(($currentSeed)).root $currentSeed $trackingIneff >& ${LOG_DIR}/090_120runSeed$currentSeed.log
    ./pythiaChargedDijet $pyconfigname 120 150 20 ${Out_DIR}/120_150runSeed$(($currentSeed)).root $currentSeed $trackingIneff >& ${LOG_DIR}/120_150runSeed$currentSeed.log
    ./pythiaChargedDijet $pyconfigname 150 200 20 ${Out_DIR}/150_200runSeed$(($currentSeed)).root $currentSeed $trackingIneff >& ${LOG_DIR}/150_200runSeed$currentSeed.log
    ./pythiaChargedDijet $pyconfigname 200 -1 20 ${Out_DIR}/200_-1runSeed$(($currentSeed)).root $currentSeed $trackingIneff >& ${LOG_DIR}/200_-1runSeed$currentSeed.log
done

echo "All jobs finished!"
