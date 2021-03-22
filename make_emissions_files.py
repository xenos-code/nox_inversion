# James East

import sys
sys.path.insert(0, '.')
from beta_calc_levs import emissions_sums
from datetime import date, timedelta

#OPTIONS
start_date = date(2018,7,1) # first day to get emissions
end_date = date(2018,7,31) # last day to get emissions
emisdir = '2018' # emissions dir name. not full path.
#emisdir = 'antbe0_posterior'


#SCRIPT -- shouldn't have to change anything below here --
basedir='/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/'

for d in range(int( (end_date-start_date).days )+1):
    yyyymmdd = (start_date + timedelta(d)).strftime("%Y%m%d")
    print(f'Getting date: {yyyymmdd}', flush=True)
    epathsbase = {
        'anth': basedir+'input_2018_hemi/emis/'+emisdir+'/repemis_mole_all.'+yyyymmdd+'.nc',
        'soil': basedir+'input_2018_hemi/emis/2018/CAMS-108NHEMI2-SOIL_Glb_0.5x0.5_soil_nox_v1.1_2015-07-15.nc',
        'ship': basedir+'input_2018_hemi/emis/2018/repemis_mole_all_shipping.'+yyyymmdd+'.nc',
        'fire': basedir+'input_2018_hemi/emis/2018/emis_mole_3d_finnfires_'+yyyymmdd+'_HEMI_108k.ncf',
        'lnox': basedir+'input_2018_hemi/emis/2018/emis_mole_lightning_20170711_HEMI_108k_cmaq_cb6_2017ga_hemi_cb6_17jh.ncf'
    }

    noxemis = emissions_sums(epathsbase)
    noxemisout.append(noxemis)

emisbasef = f'noxemis_{emisdir}_{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}.nc'
(xr.concat([xr.merge([earray]) for earray in noxemisbaseout], dim='TSTEP')).to_netcdf(emisbasef)




