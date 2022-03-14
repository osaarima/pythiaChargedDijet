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
    ./pythiaChargedDijet $pyconfigname 10 -1 20 ${Out_DIR}/runSeed$currentSeed.root $currentSeed $trackingIneff >& ${LOG_DIR}/runSeed$currentSeed.log
done

echo "All jobs finished!"
