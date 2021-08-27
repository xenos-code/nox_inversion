#!/bin/bash

emisdir=${1}
basedir=/work/MOD3EVAL/jeast/NO2ASSIM/CMAQ/input/2019_hemi/emis/2019_inversion/

cd ${basedir}/${emisdir}
for f in $(ls *.nc)
do
    echo $f
    nccopy -7 -d 1 -c TSTEP/1,LAY/1,ROW/187,COL/187 $f ${f}.tmp
    mv $f ${f}.uncompressed
    mv ${f}.tmp $f
done

