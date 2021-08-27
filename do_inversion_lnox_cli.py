#!/usr/bin/env python
# coding: utf-8

# Calculate beta and create files for scaling NOx emissions

import sys
sys.path.insert(0, '/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/scripts/hemi')
from datetime import date, datetime
from beta_calc_levs import *
from calendar import monthrange

# excessive number of sys args...
outdir   = sys.argv[1] # output directory for figures and analysis file
datadir   = sys.argv[2] # dir where pre-computed files are stored
emisdir  = sys.argv[3] # name of unperturbed emissions
base     = sys.argv[4] # name of non-GSI/std run 
cut     = sys.argv[5] # name (APPL) of perturbed/cut run
gsirun   = sys.argv[6] # name of GSI run
#betafile = sys.argv[6] # file used for beta IF calcbeta=False
month    = int(sys.argv[7]) # number of the month
toplev   = int(sys.argv[8])
#calcbeta = sys.argv[9].lower() == 'true' # true or false, case insensitive


def do_nox_inversion(month):
    '''
    month: int of month (e.g. 6 for June)
    
    creates plots, saved to cwd 
    returns: analysis file for scaling emissions as dataset
    '''
    #monthcurrent = date(2019,month,1)
    #monthnext = date(2019,month+1,1)
    #lenmonth = (monthnext - monthcurrent).days
    lenmonth = monthrange(2019,month)[1]
    
    days=lenmonth # How many days do you want?
    edayi=0
    edayj=days*25
    
    st = lenmonth-days+1
    start_date = date(2019,month,st)
    end_date = date(2019,month,lenmonth)
    datestr=f'{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}'
    
    print('Calculating a new beta!', flush=True)
    # Open file created in monthly_beta() call
    beta = beta_monthly(
        start_date,
        end_date,
        datadir,
        lok=0,
        hik=toplev,
        ltng=True,
        concdir=base, # base case
        cutdir=cut, # perturbed case
        emisbase=emisdir,
        emisperturb=None,
        **{'hourly':False,
           'cutfrac':0.15, # cutfrac for lnox=0.15
           'slimit':False,#True,
           #'s_lim':0.01,
           'min_limit': 0.01,
           'max_limit': 10
          }
    )

    ##
    # Calc LNOX emissions update
    ##

    emisfname = f'{datadir}/noxemis_{emisdir}_{datestr}.nc'
    noxemis = open_emis(emisfname, start_date, end_date)

    # get VCDs and the difference
    levs=toplev
    gsirunvcd = get_columns(f'{datadir}/vcd_partial_{levs}L_{gsirun}_2019{month:02}.nc',start_date, end_date)
    basevcd = get_columns(f'{datadir}/vcd_partial_{levs}L_{base}_2019{month:02}.nc', start_date, end_date)
    diff = (gsirunvcd-basevcd).mean(dim='TSTEP').squeeze()
    diffrel = diff/(basevcd).mean(dim='TSTEP').squeeze()
    dEE = beta*diffrel
    dEE = np.where(dEE.mask, 0, dEE) # make all NaNs zero, there should not be NaNs in emis
                                     # file. So CMAQ doesn't fail when reading emis file

    dE = dEE*noxemis['lnox'].mean(dim='TSTEP').squeeze()
    
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
    dE.attrs['units'] = 'moles/s'
    dE.attrs['var_desc'] = 'mass emissions change'
    
    outf = xr.merge([beta, diff, diffrel, dEE, dE])
    datestr = f'{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}'
    monthstr = f'{start_date.strftime("%Y%m")}' 
    
    # Plots
    btitle=rf'$\beta$, {betaa.mean():0.2f} $\pm$ {betaa.std():0.2f} '\
           rf'({betaa.min():0.2f}, {betaa.max():0.2f})'
    plot_hemi(beta, title=btitle, lo=0.3, hi=1.7, cmap='Spectral_r',
              save=True, fname=f'{outdir}/figs/beta_{gsirun}_{datestr}.png')
    
    plot_hemi(diff,title=f'$\Delta \Omega$ mean filtered for overpass',
              lo=-0.5e14, hi=0.5e14, cmap='RdBu_r', units='molecules cm-2',
              save=True, fname=f'{outdir}/figs/analysis_increment_{gsirun}_{datestr}.png')
    
    plot_hemi(diffrel,title=r'$\Delta \Omega / \Omega$',
              lo=-0.5,hi=0.5, cmap='RdBu_r', units='fraction',
              save=True, fname=f'{outdir}/figs/analysis_increment_rel_{gsirun}_{datestr}.png')
    
    # Calculate and plot the emissions change
    plot_hemi(dEE,title=r'$\frac{\Delta E}{E}$',units='fraction',
              lo=-0.1,hi=0.1,cmap='RdBu_r',
              save=True, fname=f'{outdir}/figs/emis_change_rel_{gsirun}_{datestr}.png')
    
    # emis change mass
    plot_hemi(dE, title=r'$\Delta$ E', units='moles/s',
              lo=-2,hi=2, cmap='RdBu_r',
              save=True, fname=f'{outdir}/figs/emis_change_{gsirun}_{datestr}.png')

    #return outf

    outpath = f'{outdir}/{gsirun}_{start_date.strftime("%Y%m")}_inversion_analysis.nc'

    outf.attrs['file_creation_date'] = datetime.today().strftime('%Y-%m-%d_%H:%M:%S')
    outf.attrs['file_source_script'] = '/work/MOD3EVAL/jeast/NO2ASSIM/CMAQ/scripts/hemi/nox_inversion/do_inversion_lnox_cli.py'

    outf.to_netcdf(path=outpath)
    print(f'Created file: {outpath}', flush=True)

do_nox_inversion(month)
