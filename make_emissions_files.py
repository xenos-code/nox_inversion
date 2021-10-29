# James East

import sys
sys.path.append('.')
from beta_calc_levs import emissions_sums
from datetime import date, datetime, timedelta
import xarray as xr
from calendar import monthrange
from glob import glob

#OPTIONS
#start_date = date(2018,7,1) # first day to get emissions
#end_date = date(2018,7,31) # last day to get emissions
#emisdir = '2018' # emissions dir name. not full path.
#emisdir = 'antbe0_posterior'
#emisdir = 'antbe2_posterior'
#emisdir = 'antbe3_posterior'
#emisdir = 'antbe4_posterior'
#basedir='/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/'

emisdir = sys.argv[1]
datadir = sys.argv[2]
basedir = sys.argv[3]
yyyymm = sys.argv[4]

#SCRIPT -- shouldn't have to change anything below here --

start_date = datetime.strptime(yyyymm,'%Y%m')
mlength = monthrange(start_date.year, start_date.month)[1]
end_date = (start_date + timedelta(mlength-1))
mstr=start_date.strftime('%m')

noxemisout = []
#for d in range(int( (end_date-start_date).days )+1):
for d in range(mlength):
    yyyymmdd = (start_date + timedelta(d)).strftime("%Y%m%d")
    print(f'Getting date: {yyyymmdd}', flush=True)
    epathsbase = {
        'anth': basedir+'input/2019_hemi/emis/2019_inversion/'+emisdir+'/repemis_mole_all.'+yyyymmdd+'.nc',
        'soil': basedir+'input/2019_hemi/emis/2019_inversion/'+emisdir+f'/CAMS-108NHEMI2-SOIL_Glb_0.5x0.5_soil_nox_v2.1_2018-{mstr}-15.nc',
        'ship': basedir+'input/2019_hemi/emis/2019_inversion/'+emisdir+'/repemis_mole_all_shipping.'+yyyymmdd+'.nc',
        'fire': basedir+'input/2019_hemi/emis/2019_inversion/'+emisdir+'/emis_mole_3d_finnfires_'+yyyymmdd+'_HEMI_108k.ncf',
        'lnox': glob(basedir+'input/2019_hemi/emis/2019_inversion/'+emisdir+f'/emis_mole_lightning_2017{mstr}??_HEMI_108k_cmaq_cb6_2017ga_hemi_cb6_17jh*ncf')[0]
        #'anth': basedir+'input_2018_hemi/emis/'+emisdir+'/repemis_mole_all.'+yyyymmdd+'.nc',
        #'soil': basedir+'input_2018_hemi/emis/2018/CAMS-108NHEMI2-SOIL_Glb_0.5x0.5_soil_nox_v1.1_2015-07-15.nc',
        #'ship': basedir+'input_2018_hemi/emis/2018/repemis_mole_all_shipping.'+yyyymmdd+'.nc',
        #'fire': basedir+'input_2018_hemi/emis/2018/emis_mole_3d_finnfires_'+yyyymmdd+'_HEMI_108k.ncf',
        #'lnox': basedir+'input_2018_hemi/emis/2018/emis_mole_lightning_20170711_HEMI_108k_cmaq_cb6_2017ga_hemi_cb6_17jh.ncf'
    }

    noxemis = emissions_sums(epathsbase)
    noxemisout.append(noxemis)

emisbasef = f'{datadir}/noxemis_{emisdir}_{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}.nc'
(xr.concat([xr.merge([earray]) for earray in noxemisout], dim='TSTEP')).to_netcdf(emisbasef)




