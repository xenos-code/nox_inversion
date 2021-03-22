#!/bin/bash

# A script to scale emissions based on delta_E inference
#have to make it a little special for first iteration bc different directories
baseemisdir=${1}
emisnew=${2}
analysisfile=${3}

cd ${baseemisdir} #base emis dir
cd scripts

# only anthropogenic
infs1=$(ls repemis_mole_all.2017.Jul.???.nc)
#infs2=$(ls repemis_mole_all_shipping.2017.Jul.???.nc)
#cd ../2018
#infs3=$(ls emis_mole_3d_finnfires_201807??_HEMI_108k.ncf)
#infs4=CAMS-108NHEMI2-SOIL_Glb_0.5x0.5_soil_nox_v1.1_2015-07-15.nc
#infs5=emis_mole_all_20170704_HEMI_108k_nobeis_withrwc_2017ga_hemi_cb6_17jh.ncf #for july 4th

cd /work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/scripts/hemi/antbe0

ncap2 -O -v -s 'scalefactor=EMISDELR+1' ${analysisfile} scalefactor.nc #create scalefactor var and extract only that var

#for emisfile in $infs1 $infs2 $infs3 $infs4 $infs5
for emisfile in $infs1
do
echo '-----'
echo $emisfile
echo '-----'


case $emisfile in

    CAMS*)
    inpath=${baseemisdir}/2015
    ncks -O -v NO ${inpath}/${emisfile} ${emisfile}.tmp #extract NOx for multiplication
    ncks -A -v scalefactor scalefactor.nc ${emisfile}.tmp  #should append vars from scalefactor.nc into new emis.tmp file
    ncap2 -O -s 'where(NO>0) NO=NO*scalefactor' ${emisfile}.tmp ${emisfile}.tmp2 #scale NOx vars in tmp file
    cp ${inpath}/${emisfile} ${baseemisdir}/${emisnew}/${emisfile} 
    ncks -A -h -v NO ${emisfile}.tmp2 ${baseemisdir}/${emisnew}/${emisfile} # overwrite scaled NOx vars to emis file
    rm ${emisfile}.tmp
    rm ${emisfile}.tmp2
    ;;

    *20170704*)
    inpath=${baseemisdir}/2017
    ;;

    *finnfires*)
    inpath=${baseemisdir}/2018
    ;;

    *repemis*)
    inpath=${baseemisdir}/scripts
    ;;

    *)
    echo Uh oh!
    ;;

esac

echo ${emisfile}
if [[ ${emisfile} != *"CAMS"* ]]; then # if it's not the CAMS file then
ncks -O -v NO,NO2,HONO ${inpath}/${emisfile} ${emisfile}.tmp #extract NOx for multiplication
ncks -A -v scalefactor scalefactor.nc ${emisfile}.tmp #should append vars from scalefactor.nc into new emis.tmp file
# Dont create emissions where they dont exist
# Fixes nans in output
cat > scalenox.nco << EOF
where(NO>0) NO=NO*scalefactor;
where(NO2>0) NO2=NO2*scalefactor;
where(HONO>0) HONO=HONO*scalefactor;
EOF
ncap2 -S scalenox.nco ${emisfile}.tmp #scale NOx vars in tmp file
cp ${inpath}/${emisfile} ${baseemisdir}/${emisnew}/${emisfile} 
ncks -A -h -v NO,NO2,HONO ${emisfile}.tmp ${baseemisdir}/${emisnew}/${emisfile} # overwrite scaled NOx vars to emis file
rm ${emisfile}.tmp
fi

done
