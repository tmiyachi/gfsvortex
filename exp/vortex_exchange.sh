#!/bin/bash
#
# DESCRIPTION
#  shell cript for Vortex Relocation
# 
# PARAMETER
#  clon_obs,clat_obs - Observed TC center to use first-guess TC center of input data
#  clon_new,clat_new - New TC center to relocate TC center
#
# INPUT FILE
#  vortex_move_namelist - namelist file
#  input.sigma          - input gfs sigma file
#
# OUTPUT FILE
#  output.sigma         - output gfs sigma file  
#

TOPDIR=$HOME/work/gfsvortex
BESTTRACKDIR=/mnt/drobo1/Public/cxmldata/track/besttrack

TCNAME=${TCNAME:-17PARMA}
IDATE=${IDATE:-2009093012}

JCAP=${JCAP:-190}
MN=${MN:-}
EXPENV=ncep
EXPVRT=ganal
OEXP=${EXPENV}env_${EXPVRT}vrt

TSTMP="${IDATE:0:4}-${IDATE:4:2}-${IDATE:6:2} ${IDATE:8:2}:00:00"
clon_obs=`cat ${BESTTRACKDIR}/${IDATE:0:4}/${TCNAME}.tsv | grep "${TSTMP}" | cut -f3`
clat_obs=`cat ${BESTTRACKDIR}/${IDATE:0:4}/${TCNAME}.tsv | grep "${TSTMP}" | cut -f4`
relocate='vrt'

if [ $EXPVRT = 'ganal' -o $EXPVRT = 'ganal' ]; then
    o3merge=.false.
    cwcmerge=.false.
fi

if [ $EXPENV = 'ncep' -a $JCAP = 382 ]; then
    COMINENV=${COMINENV:-/mnt/drobo1/Public/gfsinit_cfsr/${IDATE:0:6}/siganl${MN}.${IDATE}.t${JCAP}}
else
    COMINENV=${COMINENV:-$HOME/gfsinit/${EXPENV}/${IDATE}/siganl${MN}.${IDATE}.t${JCAP}}
fi
if [ $EXPVRT = 'ncep' -a $JCAP = 382 ]; then
    COMINVRT=${COMINVRT:-/mnt/drobo1/Public/gfsinit_cfsr/${IDATE:0:6}/siganl${MN}.${IDATE}.t${JCAP}}
else
    COMINVRT=${COMINVRT:-$HOME/gfsinit/${EXPVRT}/${IDATE}/siganl${MN}.${IDATE}.t${JCAP}}
fi

COMOUT=${COMOUT:-$TOPDIR/gfsinit/${OEXP}/${IDATE}/siganl${MN}.${IDATE}.t${JCAP}}
DATA=/tmp/$(whoami)/vortex_exchange
PGM=$TOPDIR/exec/vortex_exchange

if [ ! -f ${COMINENV} ]; then
    echo "cannot find ${COMINENV}"
    exit 98
fi
if [ ! -f ${COMINVRT} ]; then
    echo "cannot find ${COMINVRT}"
    exit 97
fi

OUTDIR=${COMOUT%/*}
if [ ! -d ${OUTDIR} ]; then
    mkdir -p ${OUTDIR}
fi

if [ ! -d $DATA ]; then
    mkdir -p $DATA
fi

cd $DATA
ln -f -s $COMINENV $DATA/inputenv.sigma 
ln -f -s $COMINVRT $DATA/inputvortex.sigma 
ln -f -s $COMOUT $DATA/output.sigma

cat <<EOF > vortex_move_namelist
 &tcinfo
  clon_obs=${clon_obs},clat_obs=${clat_obs} 
 /
 &param
  clon_new=${clon_new:--9999},clat_new=${clat_new:--9999},relocate="${relocate}",
  o3merge=${o3merge:-.true.},cwcmerge=${cwcmerge:-.true.}
 /
 &end
EOF

ulimit -s unlimited
$PGM
if [ $? -ne 0 ]; then
    exit 99
fi
#rm -r -f $DATA
