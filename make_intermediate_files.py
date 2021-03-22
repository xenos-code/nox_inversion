# James East

import sys
sys.path.insert(0, '/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/scripts/hemi')
from beta_calc_levs import *
from datetime import date

basedir = '/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/'

metcro2df = glob(basedir+'input_2018_hemi/mcip/METCRO2D.108NHEMI2.44L.1807??')
metcro2df.sort()
dmet2d = xr.open_mfdataset(metcro2df, combine='nested', concat_dim='TSTEP')
dmet2d = dmet2d.load()

def save_vcds(dirname, toplev, outname):
    '''
    Get all of July or June
    '''
    files = glob(basedir+'output_2018_hemi/'+dirname+'/CCTM_CONC_*v532*GSI_201807*.nc')
    files.sort()
    d = xr.open_mfdataset(files, combine='nested', concat_dim='TSTEP')
    vcd = tovcd_partial(d.NO2,dmet2d,d,0,toplev)

    # variable attrs
    vcd.attrs = d.NO2.attrs
    vcd.attrs['units'] = 'molecules cm-2'
    vcd.attrs['var_desc'] = f'NO2 vertical column density, levels 0 to {toplev}'

    # add TFLAG
    dout = xr.merge([vcd,d.TFLAG.isel(VAR=slice(0,1))])

    # dataset attrs
    dout = dout.rename({'NO2':'NO2_VCD'})
    dout.attrs = d.attrs
    dout.attrs['SDATE'] = np.int32(2018182) #2018182=July / June or July 1, 2018
    dout.attrs['NLAYS'] = np.int32(1)
    dout.attrs['NVARS'] = np.int32(1)
    dout.attrs['VAR-LIST'] = 'NO2_VCD'
    dout.attrs['file_creation_date'] = date.today().strftime('%Y%m%d')
    dout.attrs['file_source_script'] = '/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/scripts/hemi/make_intermediate_files.py'
    dout.attrs['source_conc_files'] = files

    # save file
    dout.to_netcdf(outname)

#stdL-2
#antcutL-2
#std
#lnoxcut
#antbe0
#lnoxbe0

#save_vcds(dirname='stdL-2', toplev=10, outname='vcd_partial_10L_stdL-2_201807.nc')
#save_vcds(dirname='antcutL-2', toplev=10, outname='vcd_partial_10L_antcutL-2_201807.nc')
#save_vcds(dirname='antbe0', toplev=10, outname='vcd_partial_10L_antbe0_201807.nc')
#
#save_vcds(dirname='stdL-2', toplev=15, outname='vcd_partial_15L_stdL-2_201807.nc')
#save_vcds(dirname='antcutL-2', toplev=15, outname='vcd_partial_15L_antcutL-2_201807.nc')
#save_vcds(dirname='antbe0', toplev=15, outname='vcd_partial_15L_antbe0_201807.nc')
#
#save_vcds(dirname='stdL-2', toplev=20, outname='vcd_partial_20L_stdL-2_201807.nc')
#save_vcds(dirname='antcutL-2', toplev=20, outname='vcd_partial_20L_antcutL-2_201807.nc')
#save_vcds(dirname='antbe0', toplev=20, outname='vcd_partial_20L_antbe0_201807.nc')
#
#save_vcds(dirname='std', toplev=44, outname='vcd_partial_44L_std_201807.nc')
#save_vcds(dirname='lnoxcut', toplev=44, outname='vcd_partial_44L_lnoxcut_201807.nc')
#save_vcds(dirname='lnoxbe0', toplev=44, outname='vcd_partial_44L_lnoxbe0_201807.nc')

#save_vcds(dirname='stdL-2', toplev=44, outname='vcd_partial_44L_stdL-2_201807.nc')
#save_vcds(dirname='antcutL-2', toplev=44, outname='vcd_partial_44L_antcutL-2_201807.nc')
#save_vcds(dirname='antbe0', toplev=44, outname='vcd_partial_44L_antbe0_201807.nc')


#save_vcds(dirname='stdL-2', toplev=10, outname='vcd_partial_10L_stdL-2_201806.nc')
#save_vcds(dirname='antcutL-2', toplev=10, outname='vcd_partial_10L_antcutL-2_201806.nc')
#save_vcds(dirname='antbe0', toplev=10, outname='vcd_partial_10L_antbe0_201806.nc')

#save_vcds(dirname='stdL-2', toplev=15, outname='vcd_partial_15L_stdL-2_201806.nc')
#save_vcds(dirname='antcutL-2', toplev=15, outname='vcd_partial_15L_antcutL-2_201806.nc')
#save_vcds(dirname='antbe0', toplev=15, outname='vcd_partial_15L_antbe0_201806.nc')

#save_vcds(dirname='stdL-2', toplev=20, outname='vcd_partial_20L_stdL-2_201806.nc')
#save_vcds(dirname='antcutL-2', toplev=20, outname='vcd_partial_20L_antcutL-2_201806.nc')
#save_vcds(dirname='antbe0', toplev=20, outname='vcd_partial_20L_antbe0_201806.nc')

#save_vcds(dirname='std', toplev=44, outname='vcd_partial_44L_std_201806.nc')
#save_vcds(dirname='lnoxcut', toplev=44, outname='vcd_partial_44L_lnoxcut_201806.nc')
#save_vcds(dirname='lnoxbe0', toplev=44, outname='vcd_partial_44L_lnoxbe0_201806.nc')

#save_vcds(dirname='stdL-2', toplev=44, outname='vcd_partial_44L_stdL-2_201806.nc')
#save_vcds(dirname='antcutL-2', toplev=44, outname='vcd_partial_44L_antcutL-2_201806.nc')
#save_vcds(dirname='antbe0', toplev=44, outname='vcd_partial_44L_antbe0_201806.nc')

#save_vcds(dirname='std2', toplev=20, outname='vcd_partial_20L_std2_201807.nc')
save_vcds(dirname='antbe1', toplev=20, outname='vcd_partial_20L_antbe1_201807.nc')
