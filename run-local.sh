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

setenv oname        4percentage${trackingIneff}

setenv Disk         `pwd`
setenv Out_DIR      $Disk/output/${oname}/data
setenv LOG_DIR      $Disk/output/${oname}/logs

mkdir -p $Out_DIR
mkdir -p $LOG_DIR

#./pythiaChargedDijet pythia.config pTHatMin pTHatMax dijetLeadingPt <output.root> [random_seed] [tracking inefficiency]
./pythiaChargedDijet $pyconfigname 11 21 20 ${Out_DIR}/updated11-21TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated11-21TrakcingIneff${trackingIneff}.log
./pythiaChargedDijet $pyconfigname 21 36 20 ${Out_DIR}/updated21-36TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated21-36TrakcingIneff${trackingIneff}.log
./pythiaChargedDijet $pyconfigname 36 57 20 ${Out_DIR}/updated36-57TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated36-57TrakcingIneff${trackingIneff}.log
./pythiaChargedDijet $pyconfigname 57 84 20 ${Out_DIR}/updated57-84TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated57-84TrakcingIneff${trackingIneff}.log
./pythiaChargedDijet $pyconfigname 84 117 20 ${Out_DIR}/updated84-117TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated84-117TrakcingIneff${trackingIneff}.log
./pythiaChargedDijet $pyconfigname 117 152 20 ${Out_DIR}/updated117-152TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated117-152TrakcingIneff${trackingIneff}.log
./pythiaChargedDijet $pyconfigname 152 191 20 ${Out_DIR}/updated152-191TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated152-191TrakcingIneff${trackingIneff}.log
./pythiaChargedDijet $pyconfigname 191 234 20 ${Out_DIR}/updated191-234TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated191-234TrakcingIneff${trackingIneff}.log
./pythiaChargedDijet $pyconfigname 234 300 20 ${Out_DIR}/updated234-300TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated234-300TrakcingIneff${trackingIneff}.log
./pythiaChargedDijet $pyconfigname 300 -1 20 ${Out_DIR}/updated300--1TrackingIneff${trackingIneff}.root $seedNo $trackingIneff >& ${LOG_DIR}/updated300--1TrakcingIneff${trackingIneff}.log
