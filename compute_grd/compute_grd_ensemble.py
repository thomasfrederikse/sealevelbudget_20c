# ----------------------------------------
# Compute the full 5000 member ensemble of
# RSL and rad solutions for all products,
# and also compute basin-mean rsl trends
# ===========
# For each ensemble member:
# 1. Compute annual GRACE grd solution for
#    - GrIS, AIS, TWS, Glaciers
# 2. Combine with ensemble members to get
#    full fingerprint
# ----------------------------------------
from netCDF4 import Dataset
import numpy as np
import os
import pySLE
import multiprocessing as mp
import datetime as dt
import mod_gentools as gentools
import ctypes as ct
from scipy.interpolate import interp1d

def main():
    global settings
    set_settings()
    prepare_grace() # Read GRACE/GRACE FO data
    prepare_glacier_grace()
    prepare_glacier_insitu()
    prepare_icesheets()
    prepare_tws()
    set_random_numbers()
    # Run ensembles
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(compute_ensemble_member, range(settings['num_ens']))
    #out  = pool.map(compute_ensemble_member, range(nstart,nstop))
    return

def set_settings():
    global settings
    settings = {}
    if os.uname().nodename == 'MT-110180': settings['nproc'] = 4
    else: settings['nproc'] = 40

    settings['test_run_ICE6G_D'] = True
    settings['dir_data']    = os.getenv('HOME')+'/Data/'

    # Directories
    settings['dir_gwd_prep'] = settings['dir_data'] + 'Budget_20c/grd_prep/'
    settings['dir_glacier_zemp'] =  settings['dir_data'] + 'Glaciers/Zemp_2019/'
    if settings['test_run_ICE6G_D']:
        settings['dir_grd_save'] = settings['dir_data'] + 'Budget_20c/grd_ICE6G/'
    else:
        settings['dir_grd_save'] = settings['dir_data'] + 'Budget_20c/grd/'

    # Files
    settings['fn_mask'] = settings['dir_data'] +'Budget_20c/grd_prep/mask.npy'
    settings['fn_love'] = settings['dir_data'] + 'Budget_20c/grd_prep/love.npy'
    if settings['test_run_ICE6G_D']:
        settings['fn_grace'] = settings['dir_data'] + 'Budget_20c/grd_prep/ewh_GRACE_annual.nc'
    else:
        settings['fn_grace'] = settings['dir_data'] + 'Budget_20c/grd_prep/ewh_GRACE_annual_noGIA.nc'
    settings['fn_mascon_coords'] = settings['dir_data'] + 'GRACE/JPL_mascon/mascon_coords.npy'
    settings['fn_gia_ewh']       = settings['dir_data'] + 'GIA/Caron/Ensemble/ewh_ens_05.nc'
    settings['fn_marzeion_t'] = settings['dir_data']+'Glaciers/Marzeion_2015/data_for_thomas.txt'
    settings['fn_marzeion_p'] = settings['dir_data']+'Glaciers/Marzeion_2015/data_marzeion_etal_update_2015_regional.txt'
    settings['fn_parkes']     = settings['dir_data']+'Glaciers/Parkes2018/annual.csv'

    settings['fn_grd_glacier_regions'] = settings['dir_data'] + 'Budget_20c/grd_prep/grd_glacier_regions.nc'
    settings['fn_kjeldsen'] = settings['dir_data']+'IceSheets/Kjeldsen/Kjeldsen.csv'
    settings['fn_bamber']   = settings['dir_data']+'IceSheets/Bamber_2018/Bamber-etal_2018.tab'
    settings['fn_imbie']    = settings['dir_data']+'IceSheets/IMBIE/IMBIE2.txt'
    settings['fn_mouginot'] = settings['dir_data']+'IceSheets/Mouginot_2019/GrIS_total.csv'

    settings['fn_grd_dam']      = settings['dir_data']+'Budget_20c/grd_prep/grd_tws_dam.nc'
    settings['fn_grd_gwd_wada'] = settings['dir_data']+'Budget_20c/grd_prep/grd_tws_gwd_wada.nc'
    settings['fn_grd_gwd_doll'] = settings['dir_data']+'Budget_20c/grd_prep/grd_tws_gwd_doll.nc'

    settings['tstart'] = dt.datetime.now().replace(microsecond=0)

    settings['time']         = np.arange(1900,2019)
    settings['time_grace']   = np.arange(2003,2019)
    settings['time_insitu']  = np.arange(1900,2004)
    settings['num_ens'] = 100
    return

# -- Preparation routines --
def prepare_grace():
    print('   Preparing GRACE data...')
    global grace, settings
    # Mask and mascon coordinates
    mask            = np.load(settings['fn_mask'],allow_pickle=True).all()
    glb_mscn_coords = np.load(settings['fn_mascon_coords'],allow_pickle=True)

    grace = {} # Global dict with GRACE data
    file_handle = Dataset(settings['fn_grace'])
    file_handle.set_auto_mask(False)
    grace_time = file_handle.variables['years'][:]
    grace_ewh  = file_handle.variables['ewh'][np.in1d(grace_time,settings['time_grace']),:,:]
    grace_ewh_ste = file_handle.variables['ewh_ste'][np.in1d(grace_time,settings['time_grace']),:,:]
    file_handle.close()


    # Save uncertainties in mascon coordinates
    grace_mscn_map     = np.zeros(mask['land'].shape, dtype=int)
    grace_ewh_ste_mscn = np.zeros([len(settings['time_grace']), len(glb_mscn_coords)],dtype=np.float32)
    for k in range(len(glb_mscn_coords)):
        lat_acc = np.where((mask['lat'] >= glb_mscn_coords[k, 0]) & (mask['lat'] < glb_mscn_coords[k, 1]))[0]
        lon_acc = np.where((mask['lon'] >= glb_mscn_coords[k, 2]) & (mask['lon'] < glb_mscn_coords[k, 3]))[0]
        grace_mscn_map[lat_acc[0]:lat_acc[-1]+1,lon_acc[0]:lon_acc[-1]+1] = k
        grace_ewh_ste_mscn[:, k] = np.nanmean(grace_ewh_ste[:,lat_acc[0]:lat_acc[-1]+1,lon_acc[0]:lon_acc[-1] + 1],axis=(1,2))

    # Store everything in global array
    grace['ewh'] = mp_filled_float(grace_ewh)
    grace['mscn_map'] = mp_filled_int(grace_mscn_map)
    grace['ewh_ste_mscn'] = mp_filled_float(grace_ewh_ste_mscn)
    grace['mscn_coords'] = mp_filled_float(glb_mscn_coords)
    return

def prepare_tws():
    # ---------------------------------------------------
    # Prepare TWS changes from res. impoundment and GWD
    # Read the GRD fingerprints and store in shared array
    # ---------------------------------------------------
    print('   Preparing TWS data...')
    global tws
    tws = {}
    file_handle = Dataset(settings['fn_grd_dam'])
    file_handle.set_auto_mask(False)
    tws['grd_dam_rad'] = mp_filled_float(file_handle.variables['rad'][:])
    tws['grd_dam_rsl'] = mp_filled_float(file_handle.variables['rsl'][:])
    file_handle.close()

    file_handle = Dataset(settings['fn_grd_gwd_wada'])
    file_handle.set_auto_mask(False)
    tws['grd_gwd_wada_rad'] = mp_filled_float(file_handle.variables['rad'][:])
    tws['grd_gwd_wada_rsl'] = mp_filled_float(file_handle.variables['rsl'][:])
    file_handle.close()

    file_handle = Dataset(settings['fn_grd_gwd_doll'])
    file_handle.set_auto_mask(False)
    tws['grd_gwd_doll_rad'] = mp_filled_float(file_handle.variables['rad'][:])
    tws['grd_gwd_doll_rsl'] = mp_filled_float(file_handle.variables['rsl'][:])
    file_handle.close()
    return

def prepare_icesheets():
    print('   Preparing ice sheet data...')
    # Read time series:
    # - GrIS: Bamber/Mouginot/Kjeldsen
    # - AIS:  Bamber/IMBIE2/Adhikari
    global icesheets, glaciers_insitu,settings
    icesheets = {}

    imbie_raw    = np.loadtxt(settings['fn_imbie'], delimiter=';',usecols=(0,3,4))
    mouginot_raw = np.loadtxt(settings['fn_mouginot'], delimiter=';')
    kjeldsen_raw = np.loadtxt(settings['fn_kjeldsen'],delimiter=',')
    bamber_raw   = np.loadtxt(settings['fn_bamber'],delimiter='	',skiprows=21,usecols=(0,1,2,3,4,5,6))

    # 1. AIS
    icesheets['AIS_ens'] = mp_empty_float([settings['num_ens'],len(settings['time_insitu'])])
    rnd_AIS_model    = np.random.randint(0,2,settings['num_ens']) # Randomly select IMBIE/Bamber
    rnd_IS_insitu   = np.random.normal(0,1,[settings['num_ens'], len(settings['time_insitu'])]) # Random pertubation with the in situ uncertainties
    rnd_AIS_longterm = np.random.normal(0.05, 0.04, settings['num_ens'])
    acc_b_a = np.in1d(settings['time_insitu'], bamber_raw[:,0])
    acc_b_b = np.in1d(bamber_raw[:,0], settings['time_insitu'])

    # IMBIE from monhly cumsum to annual rate
    imbie_yrs =np.arange(1992,2017)
    imbie_ann = np.zeros(len(imbie_yrs))
    imbie_ste = np.zeros(len(imbie_yrs))

    for idx,yr in enumerate(imbie_yrs):
        yr_idx = np.floor(imbie_raw[:,0])==yr
        imbie_ann[idx] = imbie_raw[yr_idx,1].mean()
        imbie_ste[idx] = imbie_raw[yr_idx,2].mean()
    imbie_rate_yrs = imbie_yrs[:-1]
    imbie_rate_ann = -np.diff(imbie_ann)
    imbie_rate_ste = np.sqrt(np.abs(np.diff(imbie_ste ** 2)))
    acc_i_a = np.in1d(settings['time_insitu'], imbie_rate_yrs)
    acc_i_b = np.in1d(imbie_rate_yrs, settings['time_insitu'])
    for ens in range(settings['num_ens']):
        rate_lcl = np.ones(len(settings['time_insitu'])) * rnd_AIS_longterm[ens]
        if rnd_AIS_model[ens] == 0:
            rate_lcl[acc_b_a] = (bamber_raw[acc_b_b,3] + bamber_raw[acc_b_b,5])/-362 + rnd_IS_insitu[ens,acc_b_a] * np.sqrt(bamber_raw[acc_b_b,4]**2 + bamber_raw[acc_b_b,6]**2)/362
        elif rnd_AIS_model[ens] == 1:
            rate_lcl[acc_i_a] = -imbie_rate_ann[acc_i_b] + imbie_rate_ste[acc_i_b] * rnd_IS_insitu[ens,acc_i_a]
        icesheets['AIS_ens'][ens,:] = np.cumsum(rate_lcl)

    # 2. GrIS
    icesheets['GrIS_ens'] = mp_empty_float([settings['num_ens'],len(settings['time_insitu'])])
    acc_b_a = np.in1d(settings['time_insitu'], bamber_raw[:,0])
    acc_b_b = np.in1d(bamber_raw[:,0], settings['time_insitu'])
    acc_m_a = np.in1d(settings['time_insitu'], mouginot_raw[:,0])
    acc_m_b = np.in1d(mouginot_raw[:,0], settings['time_insitu'])
    acc_m_k = np.in1d(kjeldsen_raw[:,0], settings['time_insitu'])
    rnd_GrIS_model = np.random.randint(0,3,settings['num_ens']) # Randomly select Mouginot/Kjeldsen/Bamber

    for ens in range(settings['num_ens']):
        rate_lcl = (-kjeldsen_raw[acc_m_k,1] + rnd_IS_insitu[ens,:] * kjeldsen_raw[acc_m_k,2])/362
        if rnd_GrIS_model[ens] == 1:
            rate_lcl[acc_b_a] = (-bamber_raw[acc_b_b,1] + rnd_IS_insitu[ens,acc_b_a] * bamber_raw[acc_b_b,2])/362
        elif rnd_GrIS_model[ens] == 2:
            rate_lcl[acc_m_a] = (-mouginot_raw[acc_m_b,1] + rnd_IS_insitu[ens,acc_m_a] * mouginot_raw[acc_m_b,2])/362
        icesheets['GrIS_ens'][ens,:] = np.cumsum(rate_lcl)

    # Add peripheral glaciers to GrIS mass balance
    icesheets['GrIS_ens'] += glacier_insitu['Greenland_periphery']
    return

def prepare_glacier_grace():
    print('   Preparing Zemp et al. glacier data for GRACE...')
    global glacier_grace, settings
    # ----------------------------------------------------------------------------------------------
    # Glacier mass balance from Zemp2019: use to separate glaciers and TWS from GRACE/GRACEFO data
    # ----------------------------------------------------------------------------------------------
    glacier_grace={}
    mask = np.load(settings['fn_mask'],allow_pickle=True).all()
    zemp_time = np.arange(2002,2017)
    zemp_rate_mean  = np.zeros([len(mask['glacier_num_insitu']),15],dtype=np.float32) # 2002 - 2016
    zemp_rate_sterr = np.zeros([len(mask['glacier_num_insitu']),15],dtype=np.float32)
    zemp_flist = os.listdir(settings['dir_glacier_zemp'])
    for reg,num in enumerate(mask['glacier_num_insitu']):
        fname = settings['dir_glacier_zemp']+[i for i in zemp_flist if '_'+str(num)+'_' in i][0]
        raw_data = np.loadtxt(fname,skiprows=28,delimiter=',',usecols=(0,10,18))
        years      = raw_data[:,0]
        acc_idx = np.in1d(years,zemp_time)
        zemp_rate_mean[reg,:]  = raw_data[acc_idx,1] # Rate in gigatons
        zemp_rate_sterr[reg,:] = raw_data[acc_idx,2] # Uncertainty in gigatons
    glacier_grace['zemp_time'] = mp_filled_int(zemp_time)
    glacier_grace['zemp_rate_mean'] = mp_filled_float(zemp_rate_mean)
    glacier_grace['zemp_rate_sterr'] = mp_filled_float(zemp_rate_sterr)
    return

def prepare_glacier_insitu():
    print('   Preparing glacier data for in-situ estimates....')
    # ----------------------------------------------------------------------------------------------
    # Glaciers: Read the data from Zemp, Marzeion, Parkes for quick processing
    # ----------------------------------------------------------------------------------------------
    global glacier_insitu, settings
    glacier_insitu={}
    mask = np.load(settings['fn_mask'],allow_pickle=True).all()

    marzeion_mean  = np.zeros([18, len(settings['time_insitu'])],dtype=np.float32)
    marzeion_sterr = np.zeros([18, len(settings['time_insitu'])],dtype=np.float32)

    # Process Marzeion
    data_t = np.loadtxt(settings['fn_marzeion_t'], delimiter='  ')
    data_p = np.loadtxt(settings['fn_marzeion_p'])
    glac_mn = -(data_p[:, 1:]).T
    glac_se = (data_t[:, 20:-1]).T
    for reg in range(18):
        unc_lcl  = glac_se[reg, :]
        mean_lcl = glac_mn[reg,:]
        unc_rate = np.sqrt(np.abs(np.diff(unc_lcl**2)))
        marzeion_mean[reg,:] = np.hstack((mean_lcl[0], mean_lcl[0], mean_lcl))[:-10]
        marzeion_sterr[reg,:]  = np.hstack((unc_rate[0], unc_rate[0], unc_rate[0], unc_rate))[:-10]

    # Process Zemp
    zemp_mean  = marzeion_mean.copy()
    zemp_sterr = marzeion_sterr.copy()
    zemp_flist = os.listdir(settings['dir_glacier_zemp'])
    for reg in range(18):
        fname = settings['dir_glacier_zemp']+[i for i in zemp_flist if '_'+str(reg+1)+'_' in i][0]
        raw_data = np.loadtxt(fname,skiprows=28,delimiter=',',usecols=(0,10,18))
        years      = raw_data[:,0]
        acc_a = np.in1d(settings['time_insitu'],years)
        acc_b = np.in1d(years,settings['time_insitu'])
        zemp_mean[reg,acc_a]  = -raw_data[acc_b,1]/362 # Rate in gigatons
        zemp_sterr[reg,acc_a] = raw_data[acc_b,2]/362 # Rate in gigatons

    # Process Parkes
    data_parkes = np.loadtxt(settings['fn_parkes'], delimiter=',', skiprows=1)
    rate_parkes_global = np.zeros(len(settings['time_insitu']),dtype=np.float32)
    rate_parkes_global[2:] = -data_parkes[:-12,3]
    rate_parkes_global[:2] = -data_parkes[:10, 3].mean()

    # Generate ensemble
    glacier_insitu['ensemble'] = mp_empty_float([settings['num_ens'],17,len(settings['time_insitu'])])
    glacier_insitu['Greenland_periphery'] = mp_empty_float([settings['num_ens'],len(settings['time_insitu'])])

    glac_model   = np.random.randint(0,2,settings['num_ens']) # Randomly select Zemp/Marzeion
    ptb_parkes   = np.random.uniform(12.3/42.7, 1, settings['num_ens']) # Uniform sampling of Parkes uncertainty
    ptb_insitu   = np.random.normal(0, 1, [settings['num_ens'], len(settings['time_insitu'])]) # Random pertubation with the in situ uncertainties
    for ens in range(settings['num_ens']):
        if glac_model[ens] == 0:
            rate_lcl = zemp_mean + ptb_insitu[ens, :][np.newaxis, :] * zemp_sterr
        elif glac_model[ens] == 1:
            rate_lcl = marzeion_mean + ptb_insitu[ens, :][np.newaxis, :] * marzeion_sterr
        else: print('GURU PANIC')
        region_split = rate_lcl.mean(axis=1) / rate_lcl.mean(axis=1).sum()
        rate_lcl += region_split[:,np.newaxis] * (ptb_parkes[ens] * rate_parkes_global)[np.newaxis,:]
        glacier_insitu['ensemble'][ens,:,:] = np.cumsum(rate_lcl[mask['glacier_num']-1,:],axis=1)
        glacier_insitu['Greenland_periphery'][ens,:] = np.cumsum(rate_lcl[4,:])
    # Normalized glacier fingerprints
    file_handle = Dataset(settings['fn_grd_glacier_regions'])
    file_handle.set_auto_mask(False)
    glacier_insitu['fp_rsl'] = mp_filled_float(file_handle.variables['rsl'][:])
    glacier_insitu['fp_rad'] = mp_filled_float(file_handle.variables['rad'][:])
    file_handle.close()
    return

def set_random_numbers():
    print('   Generating random numbers for MC...')
    # Pre-define all random numbers and store in global arrays
    global random_numbers, grace, glacier_grace, settings
    random_numbers = {}
    random_numbers['tws_dam_scale']   = mp_filled_float(np.random.normal(1,0.2,[settings['num_ens']]))
    random_numbers['tws_gwd_scale']   = mp_filled_float(np.random.normal(1,0.2,[settings['num_ens']]))
    random_numbers['tws_natural_num'] = mp_filled_int(np.random.randint(0, 100, settings['num_ens']))
    random_numbers['tws_wada_doll']   = mp_filled_int(np.random.randint(0, 2, settings['num_ens']))
    random_numbers['glacier_zemp']    = mp_filled_float(np.random.normal(0,1,[settings['num_ens'],glacier_grace['zemp_rate_mean'].shape[0],glacier_grace['zemp_rate_mean'].shape[1]]))
    return

# -- Computation routines --
def compute_ensemble_member(ens):
    # ---------------------------------------------------------------------------------------
    # Generate a single ensemble member of the load and resulting GRD patterns
    # - Start with a GIA ensemble member and compute the resulting GRACE mass changes
    # - Separate these mass changes into individual contributors and solve SLE for each
    # - Merge with SLE and load from all in-situ processes to get full 1900-2018 time series
    # - Save load and GRD
    # ---------------------------------------------------------------------------------------
    global grace, random_numbers, settings
    print('   Ensemble '+str(ens+1)+'/'+str(settings['num_ens'])+' Elapsed time: ',(dt.datetime.now().replace(microsecond=0) - settings['tstart']))
    mask = np.load(settings['fn_mask'],allow_pickle=True).all()

    # Generate GRACE load realization using perturbations and GIA
    grace_load = comp_grace_grd(mask, ens)

    # Separate GRACE data into individual processes
    grace_grd = grd_contrib(grace_load, mask, ens)

    # Compute fingerprints for individual processes
    grd_ens = {}
    grd_ens['glac']  = comp_glac_grd(grace_grd, mask, ens)
    grd_ens['GrIS']  = comp_GrIS_grd(grace_grd, mask, ens)
    grd_ens['AIS']   = comp_AIS_grd(grace_grd, mask, ens)
    grd_ens['tws']   = comp_tws_grd(grace_grd, mask, ens)

    # Compute total
    grd_ens_total = {}
    grd_ens_total['barystatic'] = np.zeros(len(settings['time']))
    grd_ens_total['qglb'] = np.zeros(len(settings['time']),dtype=np.float32)
    grd_ens_total['rsl']  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    grd_ens_total['rad']  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    for process in grd_ens:
        for qty in grd_ens_total:
            grd_ens_total[qty] += grd_ens[process][qty]
    # Save ensemble
    save_ens_indiv(mask, grd_ens, ens)
    save_ens_total(mask, grd_ens_total, ens)
    return

def comp_grace_grd(mask, ens):
    global grace, settings
    np.random.seed()
    grace_ptb = np.random.normal(0,1,(grace['ewh_ste_mscn'].shape[0],grace['ewh_ste_mscn'].shape[1]))
    grace_load =  mask['land'][np.newaxis,...]*(grace['ewh']+(grace_ptb*grace['ewh_ste_mscn'])[:,grace['mscn_map']])
    if settings['test_run_ICE6G_D']==False:
        # Add GIA ewh ensemble member
        gia      = read_GIA_ens(ens, settings)
        gia_mscn = np.zeros(mask['land'].shape,dtype=np.float32)
        for k in range(len(grace['mscn_coords'])):
            lat_acc = np.where((mask['lat'] >= grace['mscn_coords'][k, 0]) & (mask['lat'] < grace['mscn_coords'][k, 1]))[0]
            lon_acc = np.where((mask['lon'] >= grace['mscn_coords'][k, 2]) & (mask['lon'] < grace['mscn_coords'][k, 3]))[0]
            weight = np.cos(np.deg2rad(mask['lat'][lat_acc])) / np.mean(np.cos(np.deg2rad(mask['lat'][lat_acc])))  # Weight by cos lat
            gia_mscn[lat_acc[0]:lat_acc[-1]+1,lon_acc[0]:lon_acc[-1]+1] = np.nanmean(weight[:,np.newaxis] * gia['ewh'][lat_acc[0]:lat_acc[-1]+1,lon_acc[0]:lon_acc[-1]+1])
        gia_mscn = gia_mscn * mask['land']
        grace_load-=(settings['time_grace'] - settings['time_grace'].mean())[:,np.newaxis,np.newaxis] * gia_mscn
    return(grace_load)

def save_ens_indiv(mask, grd_ens, ens):
    global settings
    for process in grd_ens:
        fn = settings['dir_grd_save'] + 'grd_'+process+'_'+str(ens)+'.nc'
        file_handle = Dataset(fn, 'w')
        file_handle.createDimension('x', len(mask['lon']))
        file_handle.createDimension('y', len(mask['lat']))
        file_handle.createDimension('t', len(settings['time']))
        file_handle.createVariable('x', 'f4', ('x',),zlib=True)[:] = mask['lon']
        file_handle.createVariable('y', 'f4', ('y',),zlib=True)[:] = mask['lat']
        file_handle.createVariable('t', 'i2', ('t',),zlib=True)[:] = settings['time']
        file_handle.createVariable('rsl', 'i2', ('t', 'y', 'x',),zlib=True,complevel=4)[:] = 20*grd_ens[process]['rsl']
        file_handle.createVariable('rad', 'i2', ('t', 'y', 'x',),zlib=True,complevel=4)[:] = 20*grd_ens[process]['rad']
        file_handle.createVariable('barystatic', 'i2', ('t',), zlib=True, complevel=4)[:] = 20*grd_ens[process]['barystatic']
        file_handle.createVariable('qglb', 'i2', ('t',), zlib=True, complevel=4)[:] = 20*grd_ens[process]['qglb']
        file_handle.variables['rsl'].setncattr('scale_factor', 0.05)
        file_handle.variables['rad'].setncattr('scale_factor', 0.05)
        file_handle.variables['barystatic'].setncattr('scale_factor', 0.05)
        file_handle.variables['qglb'].setncattr('scale_factor', 0.05)
        if process=='tws':
            file_handle.createVariable('barystatic_dam', 'i2', ('t',), zlib=True, complevel=4)[:] = 20*grd_ens[process]['barystatic_dam']
            file_handle.createVariable('barystatic_gwd', 'i2', ('t',), zlib=True, complevel=4)[:] = 20*grd_ens[process]['barystatic_gwd']
            file_handle.createVariable('barystatic_natural', 'i2', ('t',), zlib=True, complevel=4)[:] = 20*grd_ens[process]['barystatic_natural']
            file_handle.createVariable('qglb_dam', 'i2', ('t',), zlib=True, complevel=4)[:] = 20*grd_ens[process]['qglb_dam']
            file_handle.createVariable('qglb_gwd', 'i2', ('t',), zlib=True, complevel=4)[:] = 20*grd_ens[process]['qglb_gwd']
            file_handle.createVariable('qglb_natural', 'i2', ('t',), zlib=True, complevel=4)[:] = 20*grd_ens[process]['qglb_natural']
            file_handle.variables['barystatic_dam'].setncattr('scale_factor', 0.05)
            file_handle.variables['barystatic_gwd'].setncattr('scale_factor', 0.05)
            file_handle.variables['barystatic_natural'].setncattr('scale_factor', 0.05)
            file_handle.variables['qglb_dam'].setncattr('scale_factor', 0.05)
            file_handle.variables['qglb_gwd'].setncattr('scale_factor', 0.05)
            file_handle.variables['qglb_natural'].setncattr('scale_factor', 0.05)
        file_handle.close()
    return

def save_ens_total(mask, grd_ens_total, ens):
    global settings
    fn = settings['dir_grd_save'] + 'grd_'+str(ens)+'.nc'
    file_handle = Dataset(fn, 'w')
    file_handle.createDimension('x', len(mask['lon']))
    file_handle.createDimension('y', len(mask['lat']))
    file_handle.createDimension('t', len(settings['time']))
    file_handle.createVariable('x', 'f4', ('x',),zlib=True)[:] = mask['lon']
    file_handle.createVariable('y', 'f4', ('y',),zlib=True)[:] = mask['lat']
    file_handle.createVariable('t', 'i2', ('t',),zlib=True)[:] = settings['time']
    file_handle.createVariable('rsl', 'i2', ('t', 'y', 'x',),zlib=True,complevel=4)[:] = 20*grd_ens_total['rsl']
    file_handle.createVariable('rad', 'i2', ('t', 'y', 'x',),zlib=True,complevel=4)[:] = 20*grd_ens_total['rad']
    file_handle.createVariable('barystatic', 'i2', ('t',), zlib=True, complevel=4,least_significant_digit=4)[:] = 20*grd_ens_total['barystatic']
    file_handle.createVariable('qglb', 'i2', ('t',), zlib=True, complevel=4,least_significant_digit=4)[:] = 20*grd_ens_total['qglb']
    file_handle.variables['rsl'].setncattr('scale_factor', 0.05)
    file_handle.variables['rad'].setncattr('scale_factor', 0.05)
    file_handle.variables['barystatic'].setncattr('scale_factor', 0.05)
    file_handle.variables['qglb'].setncattr('scale_factor', 0.05)
    file_handle.close()
    return

# Individual contributors
def grd_contrib(grace_load,mask,ens):
    # ------------------------------------------------------------
    # Compute the fingerprints of the individual contributors over
    # the GRACE period
    # ------------------------------------------------------------
    global grace, glacier_grace
    grace_load_GrIS = grace_load * 1.0 * mask['GrIS'][np.newaxis,:,:]
    grace_load_AIS  = grace_load * 1.0 * mask['AIS'][np.newaxis,:,:]
    grace_load_glac =  grace_load * 1.0 * mask['glacier_mask_grace'][np.newaxis,:,:] # GRACE contribution

    # Compute Zemp et al perturbed estimate
    zemp_real = glacier_grace['zemp_rate_mean'] + glacier_grace['zemp_rate_sterr'] * random_numbers['glacier_zemp'][ens,:,:]
    zemp_real = interp1d(glacier_grace['zemp_time'],np.cumsum(zemp_real,axis=1),axis=1,kind='linear', fill_value='extrapolate')(settings['time_grace'])
    zemp_real = (zemp_real - np.mean(zemp_real,axis=1)[:,np.newaxis]) # To mm
    for reg in range(len(mask['glacier_num_insitu'])):
        reg_idx = (mask['glacier_num'] == mask['glacier_num_insitu'][reg])
        grace_load_glac = grace_load_glac + zemp_real[reg,:][:,np.newaxis,np.newaxis] * mask['glacier_scale'][reg_idx,...].squeeze()[np.newaxis,:,:]
    grace_load_tws = grace_load - grace_load_GrIS - grace_load_AIS - grace_load_glac

    # Solve SLE
    love = np.load(settings['fn_love'],allow_pickle=True).all()
    grace_grd_land = pySLE.solver(lat=mask['lat'], lon=mask['lon'], load=grace_load, time=settings['time_grace'], slm=1.0-mask['land'], love=love, lmax=359)
    grace_grd_land.solve()
    grace_grd_GrIS = pySLE.solver(lat=mask['lat'], lon=mask['lon'], load=grace_load_GrIS, time=settings['time_grace'], slm=1.0-mask['land'], love=love, lmax=359)
    grace_grd_GrIS.solve()
    grace_grd_AIS = pySLE.solver(lat=mask['lat'], lon=mask['lon'], load=grace_load_AIS, time=settings['time_grace'], slm=1.0-mask['land'], love=love, lmax=359)
    grace_grd_AIS.solve()
    grace_grd_glac = pySLE.solver(lat=mask['lat'], lon=mask['lon'], load=grace_load_glac, time=settings['time_grace'], slm=1.0-mask['land'], love=love, lmax=359)
    grace_grd_glac.solve()
    grace_grd_tws = pySLE.solver(lat=mask['lat'], lon=mask['lon'], load=grace_load_tws, time=settings['time_grace'], slm=1.0-mask['land'], love=love, lmax=359)
    grace_grd_tws.solve()

    grace_grd = {}
    grace_grd['land'] = grace_grd_land.result
    grace_grd['GrIS'] = grace_grd_GrIS.result
    grace_grd['AIS'] = grace_grd_AIS.result
    grace_grd['glac'] = grace_grd_glac.result
    grace_grd['tws'] = grace_grd_tws.result
    return(grace_grd)

def comp_glac_grd(grace_grd, mask, ens):
    # Compute glacier contribution
    global settings, random_numbers, glacier_insitu
    grd_ens = {}
    area = gentools.grid_area(mask['lat'],mask['lon'])
    rsl  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    rad  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)

    rsl[np.in1d(settings['time'],settings['time_grace']),:,:] = grace_grd['glac']['rsl'] * 1000
    rad[np.in1d(settings['time'],settings['time_grace']),:,:] = grace_grd['glac']['rad'] * 1000

    glac_insitu_rsl = np.zeros([len(settings['time_insitu']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    glac_insitu_rad = np.zeros([len(settings['time_insitu']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    #
    for reg in range(17):
        glac_insitu_rsl+= glacier_insitu['fp_rsl'][reg,...]*glacier_insitu['ensemble'][ens,reg,:][:,np.newaxis,np.newaxis]
        glac_insitu_rad+= glacier_insitu['fp_rad'][reg,...]*glacier_insitu['ensemble'][ens,reg,:][:,np.newaxis,np.newaxis]

    ovl_time = (settings['time'] == 2003)
    insitu_time =   settings['time'] < 2003

    rsl[insitu_time,:,:] = glac_insitu_rsl[:-1,:,:] - glac_insitu_rsl[-1,:,:][np.newaxis,:,:] + rsl[ovl_time,:,:][np.newaxis,:,:]
    rad[insitu_time,:,:] = glac_insitu_rad[:-1,:,:] - glac_insitu_rad[-1,:,:][np.newaxis,:,:] + rad[ovl_time,:,:][np.newaxis,:,:]
    grd_ens['rsl'] = rsl - rsl[-19:,:,:].mean(axis=0)[np.newaxis,...]
    grd_ens['rad'] = rad - rad[-19:,:,:].mean(axis=0)[np.newaxis,...]
    area_ocn = ((area * (1 - mask['land']))[np.newaxis, ...])
    area_ocn_tot = area_ocn.sum()
    grd_ens['barystatic'] = (area_ocn*rsl).sum(axis=(1,2)) / area_ocn_tot
    area_ocn = ((area * np.isfinite(mask['basin']))[np.newaxis, ...])
    area_ocn_tot = area_ocn.sum()
    grd_ens['qglb'] = (area_ocn*rsl).sum(axis=(1,2)) / area_ocn_tot
    return(grd_ens)

def comp_AIS_grd(grace_grd, mask, ens):
    grd_ens = {}
    global settings, icesheets
    area = gentools.grid_area(mask['lat'],mask['lon'])
    rsl  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    rad  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    rsl[np.in1d(settings['time'],settings['time_grace']),:,:] = grace_grd['AIS']['rsl'] * 1000
    rad[np.in1d(settings['time'],settings['time_grace']),:,:] = grace_grd['AIS']['rad'] * 1000

    # Compute uniform grd fingerprints
    rsl_uniform = gentools.field_trend(settings['time_grace'],grace_grd['AIS']['rsl'])
    scale_factor = ((rsl_uniform * (1 - mask['land']) * area).sum() / ((1 - mask['land']) * area).sum())
    rsl_uniform = rsl_uniform / scale_factor
    rad_uniform = gentools.field_trend(settings['time_grace'],grace_grd['AIS']['rad']) / scale_factor

    AIS_insitu_rsl =  icesheets['AIS_ens'][ens,:][:,np.newaxis,np.newaxis] * rsl_uniform
    AIS_insitu_rad =  icesheets['AIS_ens'][ens,:][:,np.newaxis,np.newaxis] * rad_uniform
    ovl_time = (settings['time'] == 2003)
    insitu_time =   settings['time'] < 2003
    rsl[insitu_time,:,:] = AIS_insitu_rsl[:-1,:,:] - AIS_insitu_rsl[-1,:,:][np.newaxis,:,:] + rsl[ovl_time,:,:][np.newaxis,:,:]
    rad[insitu_time,:,:] = AIS_insitu_rad[:-1,:,:] - AIS_insitu_rad[-1,:,:][np.newaxis,:,:] + rad[ovl_time,:,:][np.newaxis,:,:]
    grd_ens['rsl'] = rsl - rsl[-19:,:,:].mean(axis=0)[np.newaxis,...]
    grd_ens['rad'] = rad - rad[-19:,:,:].mean(axis=0)[np.newaxis,...]

    area_ocn = ((area * (1 - mask['land']))[np.newaxis, ...])
    area_ocn_tot = area_ocn.sum()
    grd_ens['barystatic'] = (area_ocn*rsl).sum(axis=(1,2)) / area_ocn_tot
    area_ocn = ((area * np.isfinite(mask['basin']))[np.newaxis, ...])
    area_ocn_tot = area_ocn.sum()
    grd_ens['qglb'] = (area_ocn*rsl).sum(axis=(1,2)) / area_ocn_tot
    return(grd_ens)

def comp_GrIS_grd(grace_grd, mask, ens):
    global settings, glaciers, icesheets
    grd_ens = {}
    area = gentools.grid_area(mask['lat'],mask['lon'])
    rsl  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    rad  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    rsl[np.in1d(settings['time'],settings['time_grace']),:,:] = grace_grd['GrIS']['rsl'] * 1000
    rad[np.in1d(settings['time'],settings['time_grace']),:,:] = grace_grd['GrIS']['rad'] * 1000

    # Compute uniform grd fingerprints
    rsl_uniform = gentools.field_trend(settings['time_grace'],grace_grd['GrIS']['rsl'])
    scale_factor = ((rsl_uniform * (1 - mask['land']) * area).sum() / ((1 - mask['land']) * area).sum())
    rsl_uniform = rsl_uniform / scale_factor
    rad_uniform = gentools.field_trend(settings['time_grace'],grace_grd['GrIS']['rad']) / scale_factor

    GrIS_insitu_rsl =  icesheets['GrIS_ens'][ens,:][:,np.newaxis,np.newaxis] * rsl_uniform
    GrIS_insitu_rad =  icesheets['GrIS_ens'][ens,:][:,np.newaxis,np.newaxis] * rad_uniform
    ovl_time = (settings['time'] == 2003)
    insitu_time =   settings['time'] < 2003
    rsl[insitu_time,:,:] = GrIS_insitu_rsl[:-1,:,:] - GrIS_insitu_rsl[-1,:,:][np.newaxis,:,:] + rsl[ovl_time,:,:][np.newaxis,:,:]
    rad[insitu_time,:,:] = GrIS_insitu_rad[:-1,:,:] - GrIS_insitu_rad[-1,:,:][np.newaxis,:,:] + rad[ovl_time,:,:][np.newaxis,:,:]
    grd_ens['rsl'] = rsl - rsl[-19:,:,:].mean(axis=0)[np.newaxis,...]
    grd_ens['rad'] = rad - rad[-19:,:,:].mean(axis=0)[np.newaxis,...]
    area_ocn = ((area * (1 - mask['land']))[np.newaxis, ...])
    area_ocn_tot = area_ocn.sum()
    grd_ens['barystatic']   = (area_ocn*rsl).sum(axis=(1,2)) / area_ocn_tot
    area_ocn = ((area * np.isfinite(mask['basin']))[np.newaxis, ...])
    area_ocn_tot = area_ocn.sum()
    grd_ens['qglb'] = (area_ocn*rsl).sum(axis=(1,2)) / area_ocn_tot
    return(grd_ens)

def comp_tws_grd(grace_grd,mask, ens):
    global settings, tws, random_numbers
    grd_ens = {}
    area = gentools.grid_area(mask['lat'],mask['lon'])
    rsl  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)
    rad  = np.zeros([len(settings['time']),len(mask['lat']),len(mask['lon'])],dtype=np.float32)

    rsl[np.in1d(settings['time'],settings['time_grace']),:,:] = grace_grd['tws']['rsl'] * 1000
    rad[np.in1d(settings['time'],settings['time_grace']),:,:] = grace_grd['tws']['rad'] * 1000

    # Compute natural, grd, and dam over Thompson mask
    tws_natural_rsl,tws_natural_rad = read_tws_natural(ens)
    tws_dams_rsl = random_numbers['tws_dam_scale'][ens] * tws['grd_dam_rsl']
    tws_dams_rad = random_numbers['tws_dam_scale'][ens] * tws['grd_dam_rad']
    if random_numbers['tws_wada_doll'][ens] == 0:
        tws_gwd_rsl  = random_numbers['tws_gwd_scale'][ens] * tws['grd_gwd_wada_rsl']
        tws_gwd_rad  = random_numbers['tws_gwd_scale'][ens] * tws['grd_gwd_wada_rad']
    else:
        tws_gwd_rsl  = random_numbers['tws_gwd_scale'][ens] * tws['grd_gwd_doll_rsl']
        tws_gwd_rad  = random_numbers['tws_gwd_scale'][ens] * tws['grd_gwd_doll_rad']

    # Combine three sources
    tws_insitu_rsl = tws_dams_rsl + tws_gwd_rsl + tws_natural_rsl
    tws_insitu_rad = tws_dams_rad + tws_gwd_rad + tws_natural_rad
    ovl_time = (settings['time'] == 2003)
    insitu_time =   settings['time'] < 2003
    rsl[insitu_time,:,:] = tws_insitu_rsl[:-1,:,:] - tws_insitu_rsl[-1,:,:][np.newaxis,:,:] + rsl[ovl_time,:,:][np.newaxis,:,:]
    rad[insitu_time,:,:] = tws_insitu_rad[:-1,:,:] - tws_insitu_rad[-1,:,:][np.newaxis,:,:] + rad[ovl_time,:,:][np.newaxis,:,:]
    grd_ens['rsl'] = rsl - rsl[-19:,:,:].mean(axis=0)[np.newaxis,...]
    grd_ens['rad'] = rad - rad[-19:,:,:].mean(axis=0)[np.newaxis,...]

    # Save barystatic terms, as well as individual terms
    area_ocn = ((area * (1 - mask['land']))[np.newaxis, ...])
    area_ocn_tot = area_ocn.sum()
    grd_ens['barystatic']         = (area_ocn*rsl).sum(axis=(1,2)) / area_ocn_tot
    grd_ens['barystatic_dam']     = np.hstack([(area_ocn*tws_dams_rsl).sum(axis=(1,2))    / area_ocn_tot, np.zeros(15,dtype=np.float32)*np.nan])
    grd_ens['barystatic_gwd']     = np.hstack([(area_ocn*tws_gwd_rsl).sum(axis=(1,2))/ area_ocn_tot, np.zeros(15,dtype=np.float32)*np.nan])
    grd_ens['barystatic_natural'] = np.hstack([(area_ocn*tws_natural_rsl).sum(axis=(1,2))   / area_ocn_tot, np.zeros(15,dtype=np.float32)*np.nan])

    # Save quasi-global estimate (i.e. Thompsons mask)
    area_ocn = ((area * np.isfinite(mask['basin']))[np.newaxis, ...])
    area_ocn_tot = area_ocn.sum()
    grd_ens['qglb']         = (area_ocn*rsl).sum(axis=(1,2)) / area_ocn_tot
    grd_ens['qglb_dam']     = np.hstack([(area_ocn*tws_dams_rsl).sum(axis=(1,2))    / area_ocn_tot, np.zeros(15,dtype=np.float32)*np.nan])
    grd_ens['qglb_gwd']     = np.hstack([(area_ocn*tws_gwd_rsl).sum(axis=(1,2))/ area_ocn_tot, np.zeros(15,dtype=np.float32)*np.nan])
    grd_ens['qglb_natural'] = np.hstack([(area_ocn*tws_natural_rsl).sum(axis=(1,2))   / area_ocn_tot, np.zeros(15,dtype=np.float32)*np.nan])

    # Same baseline
    grd_ens['barystatic_dam']     = grd_ens['barystatic_dam']     - grd_ens['barystatic_dam'][103]     + grd_ens['barystatic'][103]
    grd_ens['barystatic_gwd']     = grd_ens['barystatic_gwd']     - grd_ens['barystatic_gwd'][103]     + grd_ens['barystatic'][103]
    grd_ens['barystatic_natural'] = grd_ens['barystatic_natural'] - grd_ens['barystatic_natural'][103] + grd_ens['barystatic'][103]

    grd_ens['qglb_dam']     = grd_ens['qglb_dam']     - grd_ens['qglb_dam'][103]     + grd_ens['qglb'][103]
    grd_ens['qglb_gwd']     = grd_ens['qglb_gwd']     - grd_ens['qglb_gwd'][103]     + grd_ens['qglb'][103]
    grd_ens['qglb_natural'] = grd_ens['qglb_natural'] - grd_ens['qglb_natural'][103] + grd_ens['qglb'][103]
    return(grd_ens)

def read_tws_natural(ens):
    # Read ensemble member of natural TWS
    global settings
    file_handle = Dataset(settings['dir_gwd_prep'] + 'grd_tws_natural_'+str(random_numbers['tws_natural_num'][ens])+'.nc')
    file_handle.set_auto_mask(False)
    tws_natural_rsl = file_handle.variables['rsl'][:]
    tws_natural_rad = file_handle.variables['rad'][:]
    file_handle.close()
    return(tws_natural_rsl,tws_natural_rad)

def read_GIA_ens(ens,settings):
    gia = {}
    # rsl
    file_handle = Dataset(settings['fn_gia_ewh'], 'r')
    file_handle.set_auto_mask(False)
    gia['ewh'] = file_handle.variables['ewh'][ens,:,:]
    file_handle.close()
    return(gia)

# ----- Multiprocessing routines
def mp_empty_float(shape):
    shared_array_base = mp.RawArray(ct.c_float, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_empty_int(shape):
    shared_array_base = mp.RawArray(ct.c_int, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_empty_bool(shape):
    shared_array_base = mp.RawArray(ct.c_bool, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_float(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_float, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_int(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_int, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_bool(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_bool, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

if __name__ == '__main__':
    main()
