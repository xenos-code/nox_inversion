#!/usr/bin/env python
# coding: utf-8
import xarray as xr
import numpy as np
import sys
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from datetime import date, timedelta
import os.path

# Land borders for mapping:
# CNO to get boundaries for mapping
# code snippet copied from https://pypi.org/project/pycno/
import pycno
import pyproj
proj = pyproj.Proj(
  (
    '+proj=stere +lat_0=90 +lat_ts=45 +lon_0=-98 ' +
    '+x_0=10098000 +y_0=10098000 +R=6370000 +to_meter=108000 ' +
    '+no_defs'
  ),
  preserve_units=True
)
cno = pycno.cno(proj=proj, xlim=(0, 187), ylim=(0, 187))
dtz=0; dpop=0
haslon=False
lon=0

def beta_monthly(start_date, end_date, lok=0, hik=20, ltng=False, **kwargs):
    '''
    args:
        start_date: datetime date obj
        end_date: datetime date obj
        lok=0: lowest model level to include
        hik=20: highest model level to include
        ltng=False: whether this is a lightning inversion (usually not)
    
    kwargs:
        concdir=None:  base case (bad design, required but is kwarg with no default)
        cutdir=None:  perturbed case
        emisbase= : emissions directory used in base (concdir) runs
        emisperturb= : emissions directory used in perturbed (cutdir) runs
        hourly=True: whether to calc on hourly basis (vs daily)
        slimit=False: whether to limit delta Omega
        s_lim=0.01: what fraction to limit delta omega
        cutfrac=None: if given, should be 0.15, used for initial beta calc
        cutfracfile=None: if given, shoudl be inversion analysis file containing
                          relative emissions change array used for beta
        anthonly=False
        debug=False: if True, save out some netcdfs for debugging
        min_limit=None (float): optional lower limit for beta. Applied at daily level
        max_limit=None (float): optional upper limit for beta. Applied at daily level
    '''
    global dtz
    global dpop

    #basedir='/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/'
    basedir='/work/MOD3EVAL/jeast/NO2ASSIM/CMAQ/'
    
    tzf = '/work/ROMO/gis_rasters/tz_world/ioapi/tz_world_hours_HEMIS.IOAPI.nc'
    popf = '/work/ROMO/gis_rasters/gpw-v4/gpw_v4_une_atotpopbt_densy_HEMIS.IOAPI.nc'
    gadmpath = '/work/ROMO/gis_rasters/gadm36/ioapi/gadm36_HEMIS.IOAPI.nc'
    dtz = xr.open_dataset(tzf)
    dpop = xr.open_dataset(popf)
    dgadm = xr.open_dataset(gadmpath)
    id2name = eval(dgadm['ID_0'].description)
    name2id = {v: k for k, v in id2name.items()}
    isus = dgadm['ID_0'][0, 0] == name2id['United States']

    betas=[]
    noxemisbaseout=[]
    noxemisperturbout=[]
    global d
    for d in range(int( (end_date-start_date).days )+1):
        yyyymmdd = (start_date + timedelta(d)).strftime("%Y%m%d")
        yymmdd = (start_date + timedelta(d)).strftime("%y%m%d")
        yyyymm = (start_date + timedelta(d)).strftime("%Y%m")
        print(f'Current date: {yyyymmdd}', flush=True)

        concdir = kwargs.get('concdir',None)
        concf = basedir+f'scripts/hemi/vcd_partial_{hik}L_{concdir}_{yyyymm}.nc'
        concdtmp = xr.open_dataset(concf)
        concd = concdtmp.isel(TSTEP=slice((start_date.day-1+d)*25, (start_date.day+d)*25-1)) # pick out 'today'

        cutdir = kwargs.get('cutdir',None)
        cutf = basedir+f'scripts/hemi/vcd_partial_{hik}L_{cutdir}_{yyyymm}.nc'
        cutdtmp = xr.open_dataset(cutf)
        cutd = cutdtmp.isel(TSTEP=slice((start_date.day-1+d)*25, (start_date.day+d)*25-1)) # pick out 'today'

        #metcro2df = basedir+'input_2018_hemi/mcip/METCRO2D.108NHEMI2.44L.'+yymmdd
        metcro2df = basedir+f'input/2019_hemi/mcip/METCRO2D_{yyyymmdd}.nc4'
        

        anthonly = kwargs.get('anthonly',False)
        
        emisbase = kwargs.get('emisbase',None)
        emisbasef = f'noxemis_{emisbase}_{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}.nc'

        # if nox emis file doesn't exist, ERROR
        if os.path.isfile(emisbasef):
            noxemisbase = open_emis(emisbasef, start_date+timedelta(days=d), start_date+timedelta(days=d))
        else:
            sys.exit(f'File {emisbasef} is missing, it is required, create with make_emissions_file.py!')

        if not ltng:
            emisperturb=kwargs.get('emisperturb',None)

            emisperturbf = f'noxemis_{emisperturb}_{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}.nc'

            if os.path.isfile(emisperturbf):
                noxemisperturb = open_emis(emisperturbf, start_date+timedelta(days=d), start_date+timedelta(days=d))
            else:
                sys.exit(f'File {emisperturbf} is missing, it is required, create with make_emissions_file.py!')

        else: #lightning case
            noxemisperturb=None
    
            
        betas.append(
            beta_1day(concd,
                      cutd,
                      metcro2df,
                      noxemisbase,
                      noxemisperturb,
                      lok,
                      hik,
                      ltng=ltng,
                      anthonly=anthonly,
                      **kwargs)
                    )

    daily_betas = np.ma.masked_array(betas)
    beta = daily_betas.mean(0) # data array
    return beta


def beta_1day(concd,
              cutd,
              metcro2df,
              noxemisbase,
              noxemisperturb,
              lok,
              hik,
              ltng=False,
              anthonly=False,
              **kwargs):

    dmet2d=xr.open_dataset(metcro2df)
    
    if ltng:
        frac = lnox_frac(noxemisbase) #fraction of emissions that are lnox
        isvalid = True # use all cells
        isoverpass = overpass_filter(concd.NO2_VCD)#, dtz)
        isclear = cloud_filter(dmet2d)
        isvalid = isvalid & isoverpass & isclear # clear sky overpass times only
        cutfrac = 0.15
    else:
        # here better to use base or perturb emis? depends on case...
        frac, ismajorityant =  antnox_filter(noxemisperturb, uselnox=False) #fraction of emissionsthat are ant
        isoverpass = overpass_filter(concd.NO2_VCD)#, dtz)
        isurban = urban_filter(dpop)
        isclear = cloud_filter(dmet2d)

        # All filters here, sequentially:
        isany=np.ones_like(ismajorityant,dtype='bool')
        filterseq = [
            ('Base',  isany & isurban),
            ('Majority Ant Nox', ismajorityant),
            ('TROPOMI overpass', isoverpass),
            ('Clear Sky', isclear)
        ]
        
        isvalid = ismajorityant & isoverpass & isurban & isclear
    
    return calc_beta(base=concd,
                     perturb=cutd,
                     #dmet2d=dmet2d,
                     noxemisbase=noxemisbase,
                     noxemisperturb=noxemisperturb,
                     antfrac=frac,
                     isvalid=isvalid,
                     lok=lok,
                     hik=hik,
                     anthonly=anthonly,
                     ltng=ltng,
                     **kwargs
                    )


def calc_beta(base, perturb, noxemisbase, noxemisperturb, antfrac, isvalid, lok, hik, anthonly=False, ltng=False, **kwargs):
    '''
    beta = dE/E * VCD/dVCD
    '''
    basevcd = base.NO2_VCD.where(isvalid).load()
    cutvcd = perturb.NO2_VCD.where(isvalid).load()
    
    hourly = kwargs.get('hourly', True)

    debug = kwargs.get('debug', False)

    slimit = kwargs.get('slimit', False)
    
    s_lim = kwargs.get('s_lim', 0.01) #cutoff fraction
    
    cutfracfile = kwargs.get('cutfracfile',None)
    cutfrac = kwargs.get('cutfrac',None)

    min_limit = kwargs.get('min_limit', None)
    max_limit = kwargs.get('max_limit', None)

    if cutfrac is None:
        if cutfracfile is None:
            # calculate the cutfrac. Should be identical for all hours, use day avg.
            eb = noxemisbase['anth'].mean(dim='TSTEP').squeeze().values # now has dims ROW, COL
            ep = noxemisperturb['anth'].mean(dim='TSTEP').squeeze().values
            cutfrac = (ep - eb) / eb # ROW*COL array of emis perturbation
        else:
            cutfrac = (xr.open_dataset(cutfracfile))['EMISDELR']
    else:
        print(f'Using user specified CUTFRAC: cutfrac = {cutfrac}', flush=True)
    
    def hourly_beta():
        vcddvcd = basevcd / (basevcd - cutvcd)
        if slimit:
            return (dEE * vcddvcd.where((1./vcddvcd > s_lim)|(1./vcddvcd < -s_lim))).mean(dim='TSTEP').isel(LAY=0)
            #.where(np.abs(vcddvcd)<(1/s_lim))).mean(dim='TSTEP').isel(LAY=0)
        else:
            return (dEE * vcddvcd).mean(dim='TSTEP').isel(LAY=0)
    
    def daily_beta():
        vcddvcd = (basevcd.mean(dim='TSTEP').isel(LAY=0) /
                  #(basevcd.mean(dim='TSTEP') - cutvcd.mean(dim='TSTEP')).isel(LAY=0))
                  (cutvcd.mean(dim='TSTEP') - basevcd.mean(dim='TSTEP')).isel(LAY=0)) # JDE change 3/17
        if slimit:
            if debug:
                dEE.to_netcdf(f'./outputs_debug/dEE-2018-7-{d}.nc')
                vcddvcd.where( ((1./vcddvcd > s_lim) | (1./vcddvcd < -s_lim)) ).to_netcdf(f'./outputs_debug/vcddvcd-2018-7-{d}.nc')
                (dEE * vcddvcd).where( ((1./vcddvcd > s_lim) | (1./vcddvcd < -s_lim)) ).to_netcdf(f'./outputs_debug/dailybeta-2018-7-{d}.nc')
            return (dEE * vcddvcd).where( ((1./vcddvcd > s_lim) | (1./vcddvcd < -s_lim)) ) #where(np.abs(vcddvcd)<(1/s_lim))
        else:
            return (dEE * vcddvcd)

    if anthonly:
        dEE = cutfrac
        if hourly:
            beta = hourly_beta()
        else:
            beta = daily_beta()
    else:
        if hourly:
            dEE = cutfrac * antfrac.where(isvalid)
            if debug:
                antfrac.where(isvalid).to_netcdf(f'./outputs_debug/antfrac-2018-7-{d:02d}.nc')
            beta = hourly_beta()
        else:
            if ltng:
                #tmp = lnox_frac(noxemisbase.where(isvalid).sum(dim='TSTEP').squeeze())
                # don't consider fraction of emissions that are LNOx
                # considering them will lead to anthro influence
                dEE = cutfrac
            else:
                tmp = (antnox_filter(noxemisbase.where(isvalid).sum(dim='TSTEP').squeeze(), uselnox=False))[0]
                dEE = cutfrac * tmp
                del(tmp)
            if debug:
                if not ltng:
                    tmp.to_netcdf(f'./outputs_debug/antfrac-2018-7-{d}.nc')
                (xr.DataArray(cutfrac)).to_netcdf(f'./outputs_debug/cutfrac-2018-7-{d}.nc')
            beta = daily_beta()

    # if min/mx limits given, apply here
    if min_limit:
        beta.values[beta.values<min_limit] = min_limit
    if max_limit:
        beta.values[beta.values>max_limit] = max_limit
 
    return beta.to_masked_array()

def get_columns(path, first_date, last_date):
    '''
    path: file to open (pre computed column file)
    first_date: datetime
    last_date: datetime

    returns: xr.DataArray of daily average column, filtered for overpass time, 24-hr days
             Length of TSTEP dim is ndays
    '''
    
    d = xr.open_dataset(path)
    
    vcd = []
    for iday in range( first_date.day, last_date.day+1 ):
        no2vcd_tmp = d.NO2_VCD.isel( TSTEP=slice((iday-1)*25, (iday*25)-1) ) #24 hours
        op = overpass_filter(no2vcd_tmp)
        vcd.append( no2vcd_tmp.where(op).mean(dim='TSTEP') )
        del(no2vcd_tmp)
    # concat the list to make a new dimension (length is ndays)
    d_out = xr.concat(vcd,dim='TSTEP') # now a new TSTEP dimension, days
        
    return d_out # TSTEP, LAY, ROW, COL


def open_emis(path, first_date, last_date):
    '''
    path: str: file to open (pre computed emissions file)
    first_date: datetime
    last_date: datetime

    returns: xr.DataArray of daily average emissions, 24-hr days
             Length of TSTEP dim is ndays
    '''
    d = xr.open_dataset(path)
    efiles = []
    for iday in range( first_date.day, last_date.day+1 ):
        efiles.append( d.isel( TSTEP=slice((iday-1)*24, (iday*24)) ) ) # 24 hours
    # concat the list to make a new dimension (length is ndays)
    d_out = xr.concat(efiles,dim='TSTEP') # now 24 hr days
        
    return d_out#.mean(dim='TSTEP', keepdims=True) # TSTEP, LAY, ROW, COL


def overpass_filter_old(dconc, dtz):
    '''
    Calculate hour of day in local time
    use to estimate overpass filter
    TROPOMI equatorial overpass time is 1:30
    '''
    utct = dconc.TSTEP
    utcoff=xr.zeros_like(dtz.UTCOFFSET)
    utcoff.values=dtz.UTCOFFSET.values.view('<i8')/3.6e12 # because xarray interprets units as time units
    lsth = xr.zeros_like(dconc.NO2)
    for t in range(len(utct)):
        lsth[t,:,:,:] = (utcoff[0,0,:,:] + utct[t]) % 24
    lsth = (lsth.values[:,0:1,:,:]).round(0)
    isoverpass = (lsth >= 13) & (lsth <= 14) #TROPOMI
    return isoverpass


def get_lon():
    global haslon
    global lon
    if not haslon: # need to open gridcro and get lon
        haslon = True
        basedir='/work/ROMO/users/bhenders/HAQAST/NO2ASSIM/CMAQ/'
        gridcro2df = glob(basedir+'input_2018_hemi/mcip/GRIDCRO2D.108NHEMI2.44L.180701')
        gridcro2d = xr.open_dataset(gridcro2df[0])#, combine='nested', concat_dim='TSTEP')
        lon = gridcro2d.LON[0,0,:,:].load()


def overpass_filter(dconc):
    '''
    Calculate hour of day in local time
    use to estimate overpass filter
    TROPOMI equatorial overpass time is 1:30
    
    overpass_filter(dconc)
    dconc: IOAPI-shaped xarray dataArray
    returns: isoverpass, np bool array same shape as input except nLEV=1
    '''
    get_lon() # can now use lon variable
    solartime_offset = (lon/15.).round(0)
    utct = dconc.isel(TSTEP=slice(0,24)).TSTEP
    lsth = (xr.zeros_like(dconc)).load()
    for t in range(len(utct)):
        lsth[t,:,:,:] = (solartime_offset + utct[t]) % 24
    lsth = (lsth.values[:,0:1,:,:]).round(0)
    isoverpass = (lsth >= 13) & (lsth <= 14) #TROPOMI
    return isoverpass


def emissions_sums(epaths):
    # Calculate the 2d emissions in each cell
    emisds = {ename: xr.open_dataset(efile) for ename, efile in epaths.items()}
    noxemis = {}
    for ename, d in emisds.items():
        #print(ename + ' has NO')
        noxemis[ename]=d.NO.isel(TSTEP=slice(0,24)).copy()
        if 'NO2' in d:
            #print(ename + ' has NO2')
            noxemis[ename] += d.NO2.isel(TSTEP=slice(0,24)).values
        if 'HONO' in d:
            #print(ename + ' has HONO')
            noxemis[ename] += d.HONO.isel(TSTEP=slice(0,24)).values
        noxemis[ename] = noxemis[ename].sum(dim='LAY',keepdims=True)
    return noxemis


def antnox_filter(noxemis, uselnox):
    '''
    noxemis: xarray of noxemissions (precomputed)
    uselnox: bool, false if not including LNOX in nox emis total 
    '''
    # condition for majoritynox emissions anthropogenic
    if uselnox:
        # Calculate anthro fraction with lnox included in total
        noxtot=sum(a for ename,a in noxemis.items())
        antfrac=noxemis['anth']/noxtot
    else:
        # Calculate anthro fraction with lightning not included
        noxtotl = sum(a for ename,a in noxemis.items()) - noxemis['lnox']
        antfrac = noxemis['anth']/noxtotl
    return antfrac, (antfrac > 0.5).values


def lnox_frac(noxemis):
    '''
    calc fraction of lnox that is total
    dont fileter out where lnox fraction small, because
    ultimately those are the regions
    we are interested in
    '''
    noxtot = sum(a for ename,a in noxemis.items())
    lfrac = noxemis['lnox']/noxtot
    return lfrac


def urban_filter(dpop):
    # filter for urban
    ppkm2 = dpop.DENS[3:4] # 2015, TSTEP, LAY, ROW, COL
    return (ppkm2 > 15).values


def cloud_filter(dmet2d):
    # Filter for model cloud cover
    cfrac = dmet2d.CFRAC.isel(TSTEP=slice(0,24))
    return (cfrac < 0.3).values


def tovcd(x, dmet2d, dconc):
    # Calc VCD with surface file only
    pedges = (
        (dmet2d.PRSFC.values - dconc.VGTOP) *
        dconc.VGLVLS[None,:,None,None] + dconc.VGTOP
    )
    dp = -np.diff(pedges, axis=1) / 100
    hPa_to_du = (
        10 * 1.3807e-23 * 6.022e23 / 0.02894 * 273.15 / 9.80665 / 101325.
    )
    return (x * dp).sum(dim='LAY', keepdims=True) * hPa_to_du * 2.69e16


def tovcd_partial(x, dmet2d, dconc, lok, hik):
    '''
    x = array of 3D NO2
    dmet2d = xarray metrco2d file
    dconc = xarray CONC file
    lok = int, index for bottom level used for partial column
    hik = int, index for top level used for partial column
    '''
    # Calc VCD with surface file only
    x = x[:,lok:hik,:,:]
    pedges = (
        (dmet2d.PRSFC.values - dconc.VGTOP) *
        dconc.VGLVLS[None,:,None,None] + dconc.VGTOP
    )
    pedges = pedges[:,lok:hik+1,:,:]
    dp = -np.diff(pedges, axis=1) / 100
    hPa_to_du = (
        10 * 1.3807e-23 * 6.022e23 / 0.02894 * 273.15 / 9.80665 / 101325.
    )
    return (x * dp).sum(dim='LAY', keepdims=True) * hPa_to_du * 2.69e16


# Plotting stuff below here
def plot_hemi(a, title='hemi', cmap='viridis', cbar=True, show=True, save=False, retax=False,ca=None,**kwargs):
    
    if ca is not None:
        ax = plt.sca(ca)
    else:
        plt.close('all')
    lognorm = kwargs.get('lognorm',False)
    if lognorm:
        plt.pcolormesh(a, norm=colors.LogNorm()) # set vmin, vmax manually
    else:
        plt.pcolormesh(a)
    cno.draw()
    if cbar:
        cb = plt.colorbar()
    ax=plt.gca()
    ax.set_aspect('equal')
    if 'lo' in kwargs:
        plt.clim(vmin=kwargs.get('lo',None))
    if 'hi' in kwargs:
        plt.clim(vmax=kwargs.get('hi',None))
    if 'units' in kwargs:
        cb.set_label(kwargs.get('units',None))
    elif hasattr(a,'units'):
        cb.set_label(a.units)
    if 'fname' in kwargs:
        fname = kwargs.get('fname',None)
    else:
        fname = 'hemiplot.png'
    plt.set_cmap(cmap)
    plt.title(title)
    if save:
        fig=plt.gcf()
        #fig.set_size_inches((10,10))
        plt.savefig(fname, dpi=600)
    else:
        if retax:
            ax = plt.gca()
            return ax
        else:
            if show:
                plt.show()
            else:
                return


def plot_beta(beta, cbar=True):
    #print(type(beta))
    USbeta = np.ma.masked_where(~isus, beta).mean()
    Gblbeta = float(beta[:, :].mean())
    plt.pcolormesh(
    beta,
    #vmin=0.6, vmax=1.64, cmap='jet'
    vmin=0, vmax=2.0, cmap='RdBu_r'
    )
    cno.draw()
    plt.title(f'US beta: {USbeta:.2f}; Global beta: {Gblbeta:.2f}')
    if cbar:
        plt.colorbar();


def plot_opts(base, perturb, cutfrac, antfrac, isvalid, anthonly, hourly):
    ismine = np.array([True])
    fig, axx = plt.subplots(1, 4, figsize=(16,4))
    
    for axi, (key, extra) in enumerate(filterseq):
        ismine = ismine & extra
        plt.sca(axx.ravel()[axi])
        plot_beta(calc_beta(
                base=base,
                perturb=perturb,
                cutfrac=cutfrac,
                antfrac=antfrac,
                isvalid=ismine,
                anthonly=anthonly,
                hourly=anthonly), 
            cbar=False)
        plt.text(1, 1, key)
        axx[axi].set_facecolor('gray')
    cax = fig.add_axes([0.1, 0.05, 0.8, .025])
    plt.colorbar(cax=cax, orientation='horizontal');


#plot_opts(dno2, cutno2, 0.15, antfrac, True, anthonly=False, hourly=True)
if '__name__' == '__main__':
    start_date = date(2018,7,1)
    end_date = date(2018,7,31)
    beta = beta_monthly(start_date, end_date)
    glb = beta.mean()
    USbeta = np.ma.masked_where(~isus, beta).mean()
    dbeta = xr.DataArray(beta,dims=('ROW', 'COL'), name='beta')
    dbeta.to_dataset('july_beta_std_antcut.nc')
    plot_hemi(beta, title=f'US beta: {USbeta:.2f}; Global beta: {glb:.2f}', cmap='jet', save=True, lo=0.6, hi=1.7, fname='july_beta_std_antcut.png') 
