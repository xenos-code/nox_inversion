#!/usr/bin/env python
# coding: utf-8

# Calculate beta and create files for scaling NOx emissions

import sys
sys.path.insert(0, '/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/scripts/hemi')
from datetime import date
from beta_calc_levs import *

# excessive number of sys args...
outdir   = sys.argv[1] # output directory for figures and analysis file
datadir   = sys.argv[2] # dir where pre-computed files are stored
emisdir  = sys.argv[3] # name of emissions created by previous iteration, like antbe5_posterior
case     = sys.argv[4] # name (APPL) of this iteration
base     = sys.argv[5] # name of non-GSI run in this iteration
betafile = sys.argv[6] # file used for beta IF calcbeta=False
month    = int(sys.argv[7]) # number of the month
toplev   = int(sys.argv[8])
calcbeta = sys.argv[9].lower() == 'true' # true or false, case insensitive

#optional
if calcbeta:
    cutfracfile = sys.argv[10] # analysis file of previous iteration with path, like /[path]/antbe4_inversion_analysis.nc
    emisbase    = sys.argv[11] # name of emissions used as input to previous iteration, like antbe4_posterior
    concdir     = sys.argv[12] # name of non-GSI run in previous iteration, this run uses emisbase as input (e.g. std6)


#basedir = '/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/'
#case = 'antbe1'
#base = 'std2'
#emisdir='antbe0_posterior' # name of the directory containing the emissions used
    
def do_nox_inversion(month, calcbeta=False):
    '''
    month: int of month (e.g. 6 for June)
    
    creates plots, saved to cwd 
    returns: analysis file for scaling emissions as dataset
    '''
    monthcurrent = date(2018,month,1)
    monthnext = date(2018,month+1,1)
    lenmonth = (monthnext - monthcurrent).days
    
    days=lenmonth # How many days do you want?
    edayi=0
    edayj=days*25
    
    st = lenmonth-days+1
    start_date = date(2018,month,st)
    end_date = date(2018,month,lenmonth)
    datestr=f'{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}'
    
    emisfname = f'{datadir}/noxemis_{emisdir}_{datestr}.nc'
    noxemis = open_emis(emisfname, start_date, end_date)
    
    # get VCDs and the difference
    levs=toplev
    gsivcd = get_columns(f'{datadir}/vcd_partial_{levs}L_{case}_2018{month:02}.nc',start_date, end_date)
    basevcd = get_columns(f'{datadir}/vcd_partial_{levs}L_{base}_2018{month:02}.nc', start_date, end_date)
    diff = (gsivcd-basevcd).mean(dim='TSTEP').squeeze()
    diffrel = diff/(basevcd).mean(dim='TSTEP').squeeze()

    # Calculate a new beta?
    if calcbeta:
        print('Calculating a new beta!', flush=True)
        # Open file created in monthly_beta() call
        beta = beta_monthly(start_date,
                              end_date,
                              lok=0,
                              hik=levs,
                              ltng=True,
                              concdir=concdir, #'stdL-2', # base case
                              cutdir=base, #'std2', # perturbed case
                              emisbase=emisbase, #'2018',
                              emisperturb=emisdir,
                              **{'hourly':False,
                                 'cutfrac':0.15, # cutfrac for lnox=0.15
                                 'slimit':True,
                                 's_lim':0.01,
                                 'min_limit': 0.01,
                                 'max_limit': 10
                                }
                              )

    else:
        print('Using static beta. NOT calculating a new beta!', flush=True)
        d = xr.open_dataset(betafile)
        beta = np.ma.masked_where(np.isnan(d.BETA),d.BETA)
    
    dEE = beta*diffrel
    #print(f'type(dEE): {type(dEE)}')
    dEE = np.where(dEE.mask, 0, dEE) # make all NaNs zero, there should not be NaNs in emis
    #print(f'type(dEE): {type(dEE)}')
                                      # file. So CMAQ doesn't fail when reading emis file
    #print(f'type(noxemis.anth.squeeze()): {type(noxemis.anth.squeeze())}')
    dE = dEE*noxemis.anth.mean(dim='TSTEP').squeeze()
    #print(f'type(dE): {type(dE)}')
    
    betaa = beta.copy()
    beta = xr.DataArray(beta)
    beta = beta.rename('BETA')
    beta = beta.rename({'dim_0':'ROW', 'dim_1':'COL'})
    beta.attrs['units'] = 'unitless'
    beta.attrs['var_desc'] = 'unitless scaling factor beta'
    
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
    #print('')
    #print('dE array:')
    #print(dE)
    #dE = dE.rename({'dim_0':'ROW', 'dim_1':'COL'})
    dE.attrs['units'] = 'moles/s'
    dE.attrs['var_desc'] = 'mass emissions change'
    
    outf = xr.merge([beta, diff, diffrel, dEE, dE])
    datestr=f'{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}'
    monthstr =f'{start_date.strftime("%Y%m")}' 
    #outf.to_netcdf(path=f'./{case}_inversion_analysis_{monthstr}.nc')
    #print('made it here 1',flush=True)
    
    # Plots
    btitle=rf'$\beta$, {betaa.mean():0.2f} $\pm$ {betaa.std():0.2f} '\
           rf'({betaa.min():0.2f}, {betaa.max():0.2f})'
    plot_hemi(beta, title=btitle, lo=0.3, hi=1.7, cmap='Spectral_r',
              save=True, fname=f'{outdir}/beta_{case}_{datestr}.png')
    
    plot_hemi(diff,title=f'$\Delta \Omega$ mean filtered for overpass',
              lo=-0.5e14, hi=0.5e14, cmap='RdBu_r', units='molecules cm-2',
              save=True, fname=f'{outdir}/analysis_increment_{case}_{datestr}.png')
    
    plot_hemi(diffrel,title=r'$\Delta \Omega / \Omega$',
              lo=-0.5,hi=0.5, cmap='RdBu_r', units='fraction',
              save=True, fname=f'{outdir}/analysis_increment_rel_{case}_{datestr}.png')
    
    # Calculate and plot the emissions change
    plot_hemi(dEE,title=r'$\frac{\Delta E}{E}$',units='fraction',
              lo=-0.1,hi=0.1,cmap='RdBu_r',
              save=True, fname=f'{outdir}/emis_change_rel_{case}_{datestr}.png')
    
    # emis change mass
    plot_hemi(dE, title=r'$\Delta$ E', units='moles/s',
              lo=-2,hi=2, cmap='RdBu_r',
              save=True, fname=f'{outdir}/emis_change_{case}_{datestr}.png')

    #return outf

    outpath = f'{outdir}/{case}_inversion_analysis.nc'

    outf.attrs['file_creation_date'] = date.today().strftime('%Y%m%d')
    outf.attrs['file_source_script'] = '/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/scripts/hemi/nox_inversion/do_inversion_cli.py'

    outf.to_netcdf(path=outpath)
    print(f'Created file: {outpath}', flush=True)

#outf_july = do_nox_inversion(7)
do_nox_inversion(month, calcbeta=calcbeta)
#outf_july.to_netcdf(path=outpath)
