#!/usr/bin/env python
# coding: utf-8

# Calculate beta and create files for scaling NOx emissions

import sys
sys.path.insert(0, '/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/scripts/hemi')
from datetime import date
from beta_calc_levs import *

basedir = '/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/'
case = 'antbe1'
base = 'std2'

emisdir='antbe0_posterior' # name of the directory containing the emissions used
    
def do_nox_inversion(month):
    '''
    month: int of month (e.g. 6 for June)
    
    creates plots, saved to cwd 
    returns: analysis file for scaling emissions as dataset
    '''
    monthcurrent = date(2018,month,1)
    monthnext = date(2018,month+1,1)
    lenmonth = (monthnext - monthcurrent).days
    # How many days do you want?
    days=lenmonth
    edayi=0
    edayj=days*25
    
    metcro2df = sorted(glob(basedir+f'input_2018_hemi/mcip/METCRO2D.108NHEMI2.44L.180{month}??'))
    dmet2d = xr.open_mfdataset(metcro2df[-days:], combine='nested', concat_dim='TSTEP')
    
    metcro3df = sorted(glob(basedir+f'input_2018_hemi/mcip/METCRO3D.108NHEMI2.44L.180{month}??'))
    dmetcro3d = xr.open_mfdataset(metcro3df[-days:], combine='nested', concat_dim='TSTEP')
    
    gridcro2df = sorted(glob(basedir+f'input_2018_hemi/mcip/GRIDCRO2D.108NHEMI2.44L.180{month}01'))
    gridcro2d = xr.open_mfdataset(gridcro2df[-days:], combine='nested', concat_dim='TSTEP')
    
    st = lenmonth-days+1
    start_date = date(2018,month,st)
    end_date = date(2018,month,lenmonth)
    datestr=f'{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}'
    # Looks like surface differences go up to about layer 10 to 15, 20.
    
    # Open file created in monthly_beta() call
    
    beta20 = beta_monthly(start_date,
                          end_date,
                          lok=0,
                          hik=20,
                          concdir='stdL-2', # base case
                          cutdir='std2', # perturbed case
                          emisbase='2018',
                          emisperturb='antbe0_posterior',
                          **{'hourly':False,
                             'cutfracfile':'antbe0_inversion_analysis.nc',
                             'slimit':True,
                             's_lim':0.05})
    
    emisfname = f'../noxemis_{emisdir}_{datestr}.nc'
    emisfname2 = f'./noxemis_{emisdir}_{datestr}.nc'
    try:
        noxemis = open_emis(emisfname, start_date, end_date)
    except FileNotFoundError:
        noxemis = open_emis(emisfname2, start_date, end_date)
        
    
    # get VCDs and the difference
    levs=20
    gsivcd = get_columns(f'../vcd_partial_{levs}L_{case}_20180{month}.nc',start_date, end_date)
    basevcd = get_columns(f'../vcd_partial_{levs}L_{base}_20180{month}.nc', start_date, end_date)
    diff = (gsivcd-basevcd).mean(dim='TSTEP').squeeze()
    diffrel = diff/(basevcd).mean(dim='TSTEP').squeeze()
    
    dEE = beta20*diffrel
    print(f'type(dEE): {type(dEE)}')
    dEE = np.where(dEE.mask, 0, dEE) # make all NaNs zero, there should not be NaNs in emis
    print(f'type(dEE): {type(dEE)}')
                                      # file. So CMAQ doesn't fail when reading emis file
    print(f'type(noxemis.anth.squeeze()): {type(noxemis.anth.squeeze())}')
    dE = dEE*noxemis.anth.mean(dim='TSTEP').squeeze()
    print(f'type(dE): {type(dE)}')
    
    beta20a = beta20.copy()
    beta20 = xr.DataArray(beta20)
    beta20 = beta20.rename('BETA')
    beta20 = beta20.rename({'dim_0':'ROW', 'dim_1':'COL'})
    beta20.attrs['units'] = 'unitless'
    beta20.attrs['var_desc'] = 'unitless scaling factor beta'
    
    diff = diff.rename('NO2INC')
    diff.attrs['units'] = 'molecules cm-2'
    diff.attrs['var_desc'] = 'month average analysis increment'
    
    diffrel = diffrel.rename('NO2INCR') 
    diffrel.attrs['units'] = 'fraction'
    diffrel.attrs['var_desc'] = 'month average relative analysis increment'
    
    dEE = xr.DataArray(dEE)
    dEE = dEE.rename('EMISDELR')
    dEE = dEE.rename({'dim_0':'ROW', 'dim_1':'COL'})
    dEE.attrs['units'] = 'fraction'
    dEE.attrs['var_desc'] = 'relative emissions change'
    
    dE = xr.DataArray(dE)
    dE = dE.rename('EMISDEL') 
    print('')
    print('dE array:')
    print(dE)
    #dE = dE.rename({'dim_0':'ROW', 'dim_1':'COL'})
    dE.attrs['units'] = 'moles/s'
    dE.attrs['var_desc'] = 'mass emissions change'
    
    outf = xr.merge([beta20, diff, diffrel, dEE, dE])
    datestr=f'{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}'
    monthstr =f'{start_date.strftime("%Y%m")}' 
    #outf.to_netcdf(path=f'./{case}_inversion_analysis_{monthstr}.nc')
    
    # Plots
    btitle=rf'$\beta$, {beta20a.mean():0.2f} $\pm$ {beta20a.std():0.2f} '\
           rf'({beta20a.min():0.2f}, {beta20a.max():0.2f})'
    plot_hemi(beta20, title=btitle, lo=0.3, hi=1.7, cmap='Spectral_r',
              save=True, fname=f'beta_{case}_{datestr}.png')
    
    plot_hemi(diff,title=f'$\Delta \Omega$ mean filtered for overpass',
              lo=-0.5e14, hi=0.5e14, cmap='RdBu_r', units='molecules cm-2',
              save=True, fname=f'analysis_increment_{case}_{datestr}.png')
    
    plot_hemi(diffrel,title=r'$\Delta \Omega / \Omega$',
              lo=-0.5,hi=0.5, cmap='RdBu_r', units='fraction',
              save=True, fname=f'analysis_increment_rel_{case}_{datestr}.png')
    
    # Calculate and plot the emissions change
    plot_hemi(dEE,title=r'$\frac{\Delta E}{E}$',units='fraction',
              lo=-0.1,hi=0.1,cmap='RdBu_r',
              save=True, fname=f'emis_change_rel_{case}_{datestr}.png')
    
    # emis change mass
    plot_hemi(dE, title=r'$\Delta$ E', units='moles/s',
              lo=-2,hi=2, cmap='RdBu_r',
              save=True, fname=f'emis_change_{case}_{datestr}.png')

    return outf

#outf_june = do_nox_inversion(6)
outf_july = do_nox_inversion(7)
#outff = xr.concat([outf_june, outf_july],dim='TSTEP')
#outff.to_netcdf(path=f'./{case}_inversion_analysis.nc')
outf_july.to_netcdf(path=f'./{case}_inversion_analysis.nc')


