#!/bin/bash 
############################

PROG=`basename $0`
if [ $# -ne 3 ]
then
    echo "Usage: $PROG pythiaconfig.cmnd trackingIneff seedNo"
    exit;
fi

function setenv(){ export $1=$2; }

setenv pyconfigname $1
setenv trackingIneff $2
setenv seedNo $3

setenv oname        deltaM

setenv Disk         `pwd`
setenv Out_DIR      $Disk/output/${oname}/data
setenv LOG_DIR      $Disk/output/${oname}/logs

mkdir -p $Out_DIR
mkdir -p $LOG_DIR

#./pythiaChargedDijet pythia.config pTHatMin pTHatMax dijetLeadingPt <output.root> [random_seed] [tracking inefficiency]
./pythiaChargedDijet $pyconfigname 11 21 20 ${Out_DIR}/11-21${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/11-21${oname}_seed$seedNo.log
./pythiaChargedDijet $pyconfigname 21 36 20 ${Out_DIR}/21-36${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/21-36${oname}_seed$seedNo.log
./pythiaChargedDijet $pyconfigname 36 57 20 ${Out_DIR}/36-57${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/36-57${oname}_seed$seedNo.log
./pythiaChargedDijet $pyconfigname 57 84 20 ${Out_DIR}/57-84${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/57-84${oname}_seed$seedNo.log
./pythiaChargedDijet $pyconfigname 84 117 20 ${Out_DIR}/84-117${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/84-117${oname}_seed$seedNo.log
./pythiaChargedDijet $pyconfigname 117 152 20 ${Out_DIR}/117-152${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/117-152${oname}_seed$seedNo.log
./pythiaChargedDijet $pyconfigname 152 191 20 ${Out_DIR}/152-191${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/152-191${oname}_seed$seedNo.log
./pythiaChargedDijet $pyconfigname 191 234 20 ${Out_DIR}/191-234${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/191-234${oname}_seed$seedNo.log
./pythiaChargedDijet $pyconfigname 234 300 20 ${Out_DIR}/234-300${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/234-300${oname}_seed$seedNo.log
./pythiaChargedDijet $pyconfigname 300 -1 20 ${Out_DIR}/300--1${oname}_seed$seedNo.root $seedNo $trackingIneff >& ${LOG_DIR}/300--1${oname}_seed$seedNo.log
