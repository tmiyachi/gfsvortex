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

TCNAME=17PARMA
IDATE=2009093012

JCAP=${JCAP:-190}
MN=${MN:-}
EXP=ncep
OEXP=r${EXP}

TSTMP="${IDATE:0:4}-${IDATE:4:2}-${IDATE:6:2} ${IDATE:8:2}:00:00"
clon_obs=`cat ${BESTTRACKDIR}/${IDATE:0:4}/${TCNAME}.tsv | grep "${TSTMP}" | cut -f3`
clat_obs=`cat ${BESTTRACKDIR}/${IDATE:0:4}/${TCNAME}.tsv | grep "${TSTMP}" | cut -f4`
#clon_new=`echo "scale=2; ${clon_obs}+2" | bc`
#clat_new=`echo "scale=2; ${clat_obs}+2" | bc`

EXP=${EXP:-ncep}
MN=${MN:-}
OEXP=${OEXP:-r${EXP}}

if [ $JCAP = 382 -a $EXP = 'ncep' ]; then
    COMIN=${COMIN:-$HOME/gfsinit/${EXP}/${IDATE}/siganl.${IDATE}.t${JCAP}}
else
    COMIN=${COMIN:-/mnt/drobo1/Public/gfsinit_cfsr/${IDATE:0:6}/siganl${MN}.${IDATE}.t${JCAP}}
fi

COMOUT=${COMOUT:-$TOPDIR/gfsinit/${OEXP}/${IDATE}/siganl${MN}.${IDATE}.t${JCAP}}
DATA=/tmp/$(whoami)/vortex_move
PGM=$TOPDIR/exec/vortex_move

if [ ! -f ${COMIN} ]; then
    echo "cannot find ${COMIN}"
    exit 98
fi

OUTDIR=${COMOUT%/*}
if [ ! -d ${OUTDIR} ]; then
    mkdir -p ${OUTDIR}
fi

if [ ! -d $DATA ]; then
    mkdir -p $DATA
fi

cd $DATA
ln -f -s $COMIN $DATA/input.sigma 
ln -f -s $COMOUT $DATA/output.sigma

cat <<EOF > vortex_move_namelist
 &tcinfo
  clon_obs=${clon_obs},clat_obs=${clat_obs} 
 /
 &param
  clon_new=${clon_new:--9999},clat_new=${clat_new:--9999}
 /
 &end
EOF

ulimit -s unlimited
$PGM
if [ $? -ne 0 ]; then
    exit 99
fi
rm -r -f $DATA
