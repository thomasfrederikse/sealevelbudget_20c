# -------------------------------------------
# Read region list and generate ensembles
# of time series for each region
#
# These time series are subsequently merged
# into virtual stations in
# compute_virtual_station,py
#
# Full model specification
# RSL_basin = RSL_local + resvlm + dGIA +dGRD
# -------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
import multiprocessing as mp
import ctypes as ct
import mod_gentools as gentools
import glob

def main():
    set_settings()
    merge_stations()
    merge_vlm_observations()
    prepare_ensemble_storage()
    compute_region_ensembles()
    save_data()
    return

def set_settings():
    print('Defining settings...')
    global settings, mask
    settings = {}
    settings['select_latest_statlist'] = False
    settings['test_run_ICE6G_D'] = True
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 4
    else:
        settings['nproc'] = 40
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_mask'] = settings['dir_data'] +'Budget_20c/grd_prep/mask.npy'
    if settings['test_run_ICE6G_D']:
        settings['dir_grd'] = settings['dir_budget'] + 'grd_ICE6G/'
        settings['fn_gia_rsl'] = settings['dir_data']+'GIA/ICE6G_D/ICE6G_D_05.nc'
        settings['fn_alttg_data'] = settings['dir_budget'] + 'vlm/alttg_for_virstat_ice6g.npy'
        settings['fn_gps_data'] = settings['dir_budget'] + 'vlm/gps_data_ice6g.npy'
        settings['fn_region_ensembles']   = settings['dir_budget']+'region_data/region_ensembles_ice6g.npy'
    else:
        settings['dir_grd'] = settings['dir_budget'] + 'grd/'
        settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
        settings['fn_alttg_data'] = settings['dir_budget'] + 'vlm/alttg_for_virstat.npy'
        settings['fn_gps_data'] = settings['dir_budget'] + 'vlm/gps_data.npy'
        settings['fn_region_ensembles']   = settings['dir_budget']+'region_data/region_ensembles.npy'
    settings['fn_station_data'] = settings['dir_budget']+'tg/station_data.npy'
    settings['years'] = np.arange(1900, 2019)
    settings['num_ens'] = 100
    if settings['select_latest_statlist']:
        flist = glob.glob(settings['dir_budget']+'region_data/region_list*')
        cdate = np.zeros(len(flist))
        for idx, file in enumerate(flist): cdate[idx] = os.path.getmtime(file)
        settings['fn_region_list'] = flist[np.argmax(cdate)]
    else:
        settings['fn_region_list'] = settings['dir_budget'] + 'region_data/region_list_beta_march_9.npy'
    mask = np.load(settings['fn_mask'],allow_pickle=True).all()
    return

def save_data():
    global settings, region_info, region_ensemble_storage
    region_ensembles = np.zeros(6,dtype = object)
    for basin in range(6):
        region_ensembles[basin] = {}
        region_ensembles[basin]['coords']        = region_info[basin]['coords']
        region_ensembles[basin]['station_names'] = region_info[basin]['station_names']
        region_ensembles[basin]['height']        = region_info[basin]['height']
        region_ensembles[basin]['resvlm_ens']        = region_info[basin]['resvlm_ens']
        region_ensembles[basin]['tg_no_corrections'] = region_ensemble_storage[basin]['tg_no_corrections']
        region_ensembles[basin]['tg_gia']            = region_ensemble_storage[basin]['tg_gia']
        region_ensembles[basin]['tg_gia_grd']        = region_ensemble_storage[basin]['tg_gia_grd']
        region_ensembles[basin]['tg_resvlm']         = region_ensemble_storage[basin]['tg_resvlm']
        region_ensembles[basin]['tg_full']           = region_ensemble_storage[basin]['tg_full']
    np.save(settings['fn_region_ensembles'],region_ensembles)
    return

def merge_stations():
    print('Merging individual stations into regions...')
    # ------------------------------------------------------------------------
    # Merge individual stations within each region into region estimate
    # For each region estimate determine lag-1-autocorrelation for pertubation
    # ------------------------------------------------------------------------
    global settings, region_info, gia, mask
    region_list  = np.load(settings['fn_region_list'],allow_pickle=True)
    station_data = np.load(settings['fn_station_data'],allow_pickle=True).all()
    region_info = []

    for basin in range(len(region_list)):
        region_info.append({})
        region_info[basin]['station_names'] = np.zeros(len(region_list[basin]['list']),dtype=object)
        region_info[basin]['sample_coords'] = mp_empty_int([len(region_list[basin]['list']), 2])
        region_info[basin]['coords']        = mp_empty_float([len(region_list[basin]['list']),2])
        region_info[basin]['height']        = mp_empty_float([len(region_list[basin]['list']),len(settings['years'])])
        region_info[basin]['AR_sigma']      = mp_empty_float([len(region_list[basin]['list'])])
        region_info[basin]['AR_rho']        = mp_empty_float([len(region_list[basin]['list'])])
        for region in range(len(region_list[basin]['list'])):
            region_info[basin]['station_names'][region] = []
            # Merge stations sea-level data
            rsl_in_region = np.zeros([len(settings['years']),len(region_list[basin]['list'][region]['id'])])
            for station in range(len(region_list[basin]['list'][region]['id'])):
                idx = station_data['id'] == region_list[basin]['list'][region]['id'][station]
                rsl_in_region[:,station] = station_data['height_corr'][idx,:]
                region_info[basin]['station_names'][region].append(station_data['name'][idx][0])
            region_info[basin]['height'][region,:] = merge_stations_to_region(rsl_in_region)
            if np.isfinite(region_info[basin]['height'][region,:]).sum()<20: print(str(basin)+' '+str(region))
            # Set location
            region_info[basin]['coords'][region,:] = station_data['coords'][idx,:].flatten()
            # Determine sigma and rho for AR-pertubation
            reg_lcl   = region_info[basin]['height'][region].copy()
            trend_lcl = gentools.lsqtrend(settings['years'], reg_lcl)
            reg_lcl -= trend_lcl * (settings['years'] - settings['years'].mean())
            region_info[basin]['AR_rho'][region] = np.ma.corrcoef(np.ma.masked_invalid(reg_lcl[1:]), np.ma.masked_invalid(reg_lcl[:-1]))[0, 1]
            region_info[basin]['AR_sigma'][region] = np.nanstd(reg_lcl)
        region_info[basin]['sample_coords'][:,0] = np.argmin(np.abs(region_info[basin]['coords'][:,0][np.newaxis, :] - mask['lat'][:, np.newaxis]), axis=0)
        region_info[basin]['sample_coords'][:,1] = np.argmin(np.abs(region_info[basin]['coords'][:,1][np.newaxis, :] - mask['lon'][:, np.newaxis]), axis=0)
    return

def merge_vlm_observations():
    print('Merging individual VLM observations into region ensembles...')
    # -----------------------------------------------------------
    # For each region, merge individual residual VLM observations
    # into region-mean ensemble
    # -----------------------------------------------------------
    global settings, region_info
    gps_data     = np.load(settings['fn_gps_data'],allow_pickle=True)
    alttg_data   = np.load(settings['fn_alttg_data'],allow_pickle=True).all()
    region_list  = np.load(settings['fn_region_list'],allow_pickle=True)
    gps_dict = {}
    for station in range(len(gps_data)):
        if gps_data[station]['has_gps']:
            for idx, code in enumerate(gps_data[station]['gps_name']):
                gps_dict[code] = np.array([station,idx])

    alttg_dict = {}
    for region in range(len(alttg_data['id'])):
        for station in range(len(alttg_data['id'][region])):
            alttg_dict[alttg_data['id'][region][station]] = region
    # Determine ensemble of resvlm for each region
    for basin in range(len(region_list)):
        print('   '+region_list[basin]['name'])
        region_info[basin]['resvlm_ens'] = mp_empty_float([len(region_list[basin]['list']),settings['num_ens']])
        for region in range(len(region_list[basin]['list'])):
            if len(region_list[basin]['list'][region]['vlm_id'])==0: # No resvlm
                region_info[basin]['resvlm_ens'][region,:] = np.random.normal(0,1,settings['num_ens'])
            else: # VLM available
                vlm_mean_indiv  = np.zeros([settings['num_ens'],len(region_list[basin]['list'][region]['vlm_id'])])
                vlm_sterr_indiv = np.zeros([settings['num_ens'],len(region_list[basin]['list'][region]['vlm_id'])])
                for code in range(len(region_list[basin]['list'][region]['vlm_id'])):
                    if region_list[basin]['list'][region]['vlm_id'][code] == 'ALTTG':
                        alttg_idx = alttg_dict[region_list[basin]['list'][region]['id'][0]]
                        vlm_mean_indiv[:,code]  = alttg_data['resvlm_trend_ens'][:,alttg_idx]
                        vlm_sterr_indiv[:,code] = alttg_data['resvlm_sterr_AR1'][alttg_idx]
                    else:
                        gps_idx = gps_dict[region_list[basin]['list'][region]['vlm_id'][code]]
                        vlm_mean_indiv[:,code]  = gps_data[gps_idx[0]]['resvlm_mean_ens'][gps_idx[1],:]
                        vlm_sterr_indiv[:,code] = gps_data[gps_idx[0]]['resvlm_sterr_ens'][gps_idx[1],:]
                # Compute weights:
                weights = (1/vlm_sterr_indiv**2 / (1/vlm_sterr_indiv**2).sum(axis=1)[:,np.newaxis])
                resvlm_mean  = (weights*vlm_mean_indiv).sum(axis=1)
                resvlm_sterr = np.sqrt(((weights*vlm_sterr_indiv)**2).sum(axis=1))
                region_info[basin]['resvlm_ens'][region,:] = resvlm_mean + np.random.normal(0,1,settings['num_ens'])*resvlm_sterr
    return

def prepare_ensemble_storage():
    print('Allocating region ensembles...')
    global settings, region_info, region_ensemble_storage
    region_ensemble_storage = []
    # Prepare storage of time series ensemble
    for basin in range(len(region_info)):
        region_ensemble_storage.append({})
        region_ensemble_storage[basin]['tg_no_corrections'] = mp_empty_float([settings['num_ens'],len(region_info[basin]['coords']), len(settings['years'])])
        region_ensemble_storage[basin]['tg_gia']            = mp_empty_float([settings['num_ens'],len(region_info[basin]['coords']), len(settings['years'])])
        region_ensemble_storage[basin]['tg_gia_grd']        = mp_empty_float([settings['num_ens'],len(region_info[basin]['coords']), len(settings['years'])])
        region_ensemble_storage[basin]['tg_resvlm']         = mp_empty_float([settings['num_ens'],len(region_info[basin]['coords']), len(settings['years'])])
        region_ensemble_storage[basin]['tg_full']           = mp_empty_float([settings['num_ens'],len(region_info[basin]['coords']), len(settings['years'])])
    return

def compute_region_ensembles():
    print('Computing region ensembles...')
    global gia, mask, settings
    mask['area'] = gentools.grid_area(mask['lat'],mask['lon'])
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(compute_region_ens_indiv, range(settings['num_ens']))
    return

def compute_region_ens_indiv(ens):
    print('   Ensemble '+str(ens))
    global mask, settings, region_info, region_ensemble_storage
    grd_rsl_ens = read_GRD_rsl_ens(ens, settings) # Read GRD ensemble member
    gia = read_gia_ens(ens)
    for basin in range(6):
        mask_lcl = (mask['basin']==basin)
        grd_basin = np.sum(grd_rsl_ens*(mask_lcl*mask['area'])[np.newaxis,:,:],axis=(1,2))/np.sum(mask_lcl*mask['area'])
        for region in range(len(region_info[basin]['coords'])):
            # AR1 pertubation
            np.random.seed()
            tseries_ptb = rnd_num = np.random.normal(0,region_info[basin]['AR_sigma'][region],size=len(settings['years']))
            for yr in range(len(settings['years'])):
                tseries_ptb[yr] = region_info[basin]['AR_rho'][region] * tseries_ptb[yr - 1] + rnd_num[yr]
            # Compute corrections
            dGRD_region   = grd_basin - grd_rsl_ens[:,region_info[basin]['sample_coords'][region, 0],region_info[basin]['sample_coords'][region, 1]]
            dGIA_region   = (settings['years']-settings['years'].mean()) * (gia['basin_mean'][basin] - gia['rsl'][region_info[basin]['sample_coords'][region, 0],region_info[basin]['sample_coords'][region, 1]])
            resvlm_region = (settings['years']-settings['years'].mean()) * region_info[basin]['resvlm_ens'][region,ens]
            # Compute ensemble members
            region_ensemble_storage[basin]['tg_no_corrections'][ens,region,:] = region_info[basin]['height'][region,:] + tseries_ptb
            region_ensemble_storage[basin]['tg_gia'][ens,region,:]            = region_info[basin]['height'][region,:] + dGIA_region + tseries_ptb
            region_ensemble_storage[basin]['tg_gia_grd'][ens,region,:]        = region_info[basin]['height'][region,:] + dGIA_region + dGRD_region + tseries_ptb
            region_ensemble_storage[basin]['tg_resvlm'][ens,region,:]         = region_info[basin]['height'][region,:] + resvlm_region + tseries_ptb
            region_ensemble_storage[basin]['tg_full'][ens,region,:]           = region_info[basin]['height'][region,:] + dGIA_region + dGRD_region + resvlm_region + tseries_ptb
    return

def read_gia_ens(ens):
    global mask, settings
    gia = {}
    file_handle = Dataset(settings['fn_gia_rsl'], 'r')
    file_handle.set_auto_mask(False)
    gia['lat'] =  file_handle.variables['y'][:]
    gia['lon'] =  file_handle.variables['x'][:]
    if settings['test_run_ICE6G_D']:
        gia['rsl'] =  file_handle.variables['RSL'][:]
    else:
        gia['rsl'] =  file_handle.variables['rsl'][ens,:,:]
    file_handle.close()
    # Basin-mean estimate
    area = gentools.grid_area(gia['lat'],gia['lon'])
    gia['basin_mean'] = np.zeros(6)
    for basin in range(6):
        # Basin-mean GIA
        mask_lcl = (mask['basin']==basin)
        gia['basin_mean'][basin] = np.sum(gia['rsl']*mask_lcl*area)/np.sum(mask_lcl*area)
    return(gia)

def read_GRD_rsl_ens(ens,settings):
    fn = settings['dir_grd']+'grd_'+str(ens)+'.nc'
    GRD_rsl_ens = Dataset(fn,'r').variables['rsl'][:]._get_data()
    return(GRD_rsl_ens)

def merge_stations_to_region(rsl_in_region):
    while rsl_in_region.shape[1]>1:
        # Find max number of overlaps
        n_ovl = np.zeros([rsl_in_region.shape[1],rsl_in_region.shape[1]],dtype=int)
        for i in range(rsl_in_region.shape[1]):
            for j in range(rsl_in_region.shape[1]):
                if i>j: n_ovl[i,j] = np.isfinite(rsl_in_region[:,i]*rsl_in_region[:,j]).sum()
        merge_idx = np.unravel_index(np.argmax(n_ovl),n_ovl.shape)
        merge_array_lcl = rsl_in_region[:,merge_idx]
        merge_array_lcl = gentools.merge_common_mean(merge_array_lcl)
        rsl_in_region = np.hstack([rsl_in_region,merge_array_lcl[:,np.newaxis]])
        rsl_in_region = np.delete(rsl_in_region,merge_idx,axis=1)
    rsl_in_region = rsl_in_region.flatten()
    return(rsl_in_region)

# Parallel processing routines
def mp_empty_float(shape):
    shared_array_base = mp.RawArray(ct.c_float, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_empty_int(shape):
    shared_array_base = mp.RawArray(ct.c_int, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_float(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_float, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_bool(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_bool, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

if __name__ == '__main__':
    main()
