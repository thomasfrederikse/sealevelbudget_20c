# ------------------------------------------------------------------
# Compute residual VLM from alt-tg
# This program has 2 modes:
# First mode:
# Compute alt-tg for provisional region list for during TG QC
# Compute alt-tg for final region list
# ------------------------------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
import mod_gentools as gentools
import shutil
import multiprocessing as mp
import ctypes as ct
import glob

def main():
    set_settings()
    read_altimetry()
    set_alttg_list()
    compute_correlating_points()
    compute_ts()
    compute_trend()
    save_data()
    return

def set_settings():
    print('Define settings...')
    global settings
    settings = {}
    settings['region_selection'] = False # True: read from region_selection list. # False: read from final list
    settings['select_latest_statlist'] = False
    settings['test_run_ICE6G_D'] = True
    settings['years'] = np.arange(1900,2019)
    settings['min_alt_years'] = 15  # Minimum number of years of overlap between alt/tg required to compute VLM trend
    settings['min_corr'] = 0.5      # Minimum correlation between altimetry and tide gauge
    settings['max_dist'] = 300000   # Maximum distance (m) between altimetry grid cell and tide gauge location
    settings['num_ens']  = 100     # Number of ensembles

    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_scratch'] = os.getenv('HOME') + '/Scratch/'
    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 4
        settings['fn_hector'] = os.getenv('HOME') + '/Scripts/Hector/Python/MacOS/est_trend/estimatetrend'
    else:
        settings['nproc'] = 40
        settings['fn_hector'] = os.getenv('HOME') + '/Code/Hector/estimatetrend'
    settings['dir_budget']  = settings['dir_data'] + 'Budget_20c/'

    if settings['test_run_ICE6G_D']:
        settings['dir_grd'] = settings['dir_budget'] + 'grd_ICE6G/'
        settings['fn_gia_rad'] = settings['dir_data']+'GIA/ICE6G_D/ICE6G_D_05.nc'
        settings['fn_gia_rsl'] = settings['dir_data']+'GIA/ICE6G_D/ICE6G_D_05.nc'
        settings['probability'] = np.ones(settings['num_ens'])/settings['num_ens']
    else:
        settings['dir_grd'] = settings['dir_budget'] + 'grd/'
        settings['fn_gia_rad'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rad_ens_05.nc'
        settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
        settings['probability'] = Dataset(settings['fn_gia_rad'], 'r').variables['probability'][:settings['num_ens']]._get_data()
        settings['probability'] = settings['probability'] / settings['probability'].sum()

    settings['fn_altimetry'] = settings['dir_budget']+'vlm/Altimetry_annual.nc'
    settings['fn_station_data'] = settings['dir_budget']+'tg/station_data.npy'
    settings['fn_regions_for_selection'] = settings['dir_budget']+'tg/regions_for_selection.npy'
    if settings['region_selection']:
        print('   MODE: REGION SELECTION')
        settings['fn_alttg_data'] = settings['dir_budget'] + 'vlm/alttg_for_region_selection.npy'
    else:
        print('   MODE: VIRTUAL STATION')
        if settings['test_run_ICE6G_D']:
            settings['fn_alttg_data'] = settings['dir_budget'] + 'vlm/alttg_for_virstat_ice6g.npy'
        else:
            settings['fn_alttg_data'] = settings['dir_budget'] + 'vlm/alttg_for_virstat.npy'
        # REGION LIST
        if settings['select_latest_statlist']:
            flist = glob.glob(settings['dir_budget']+'region_data/region_list*')
            cdate = np.zeros(len(flist))
            for idx, file in enumerate(flist): cdate[idx] = os.path.getmtime(file)
            settings['fn_region_list'] = flist[np.argmax(cdate)]
            print(flist[np.argmax(cdate)])
        else:
            settings['fn_region_list'] = settings['dir_budget']+'region_data/region_list_beta_march_9.npy'
    return

def set_alttg_list():
    print('Filling alt-tg-information...')
    global alttg_list, settings
    alttg_list = {}
    time_alt_idx = np.in1d(settings['years'],altimetry['time'])
    if settings['region_selection']: # Read from region_selection
        regions_for_selection = np.load(settings['fn_regions_for_selection'], allow_pickle=True).all()
        alttg_list['id']     = regions_for_selection['id'].copy()
        alttg_list['coords'] = mp_filled_float(regions_for_selection['coords'])
        alttg_list['height'] = mp_filled_float(regions_for_selection['height_corr'][:,time_alt_idx])
    else: # Compute merged stations from definitive list
        station_data = np.load(settings['fn_station_data'], allow_pickle=True).all()
        region_list  = np.load(settings['fn_region_list'], allow_pickle=True)
        alttg_id = []
        alttg_coords = []
        alttg_height = []
        for basin in range(len(region_list)):
            for region in range(len(region_list[basin]['list'])):
                if 'ALTTG' in region_list[basin]['list'][region]['vlm_id']:
                    rsl_in_region = np.zeros([len(settings['years']), len(region_list[basin]['list'][region]['id'])])
                    for station in range(len(region_list[basin]['list'][region]['id'])):
                        idx = station_data['id'] == region_list[basin]['list'][region]['id'][station]
                        rsl_in_region[:,station] = station_data['height_corr'][idx,:]
                        height_lcl = merge_stations_to_region(rsl_in_region)
                    # Store
                    alttg_id.append(region_list[basin]['list'][region]['id'])
                    alttg_coords.append(station_data['coords'][idx])
                    alttg_height.append(height_lcl)
        alttg_list['id'] = np.array(alttg_id)
        alttg_list['coords'] = mp_filled_float(np.array(alttg_coords).squeeze())
        alttg_list['height'] = mp_filled_float(np.array(alttg_height)[:,time_alt_idx])
        return

def compute_correlating_points():
    global altimetry, alttg_list, settings
    # ------------------------------------------------------
    # Determine points in altimetry that correlate with tide
    # gauge sea level and store for ensemble computation
    # ------------------------------------------------------
    print('Computing correlation points...')
    alttg_list['tg_coords']     = np.zeros([len(alttg_list['id']),2],dtype=int)
    alttg_list['weight']        = np.zeros(len(alttg_list['id']),dtype=object)
    alttg_list['has_corr']      = np.zeros(len(alttg_list['id']),dtype=bool)
    alttg_list['time_acc']      = np.zeros(len(alttg_list['id']),dtype=object)
    alttg_list['tg_tseries']    = np.zeros(len(alttg_list['id']),dtype=object)
    alttg_list['alt_coords']    = np.zeros(len(alttg_list['id']),dtype=object)
    for region in range(len(alttg_list['id'])):
        time_acc  = np.isfinite(alttg_list['height'][region,:])
        if time_acc.sum()>settings['min_alt_years']:
            alttg_list['tg_coords'][region,0] = np.argmin(np.abs(altimetry['lat'] - alttg_list['coords'][region,0]))
            alttg_list['tg_coords'][region,1] = np.argmin(np.abs(altimetry['lon'] - alttg_list['coords'][region,1]))
            # Detrend TG and altimetry data set for correlation
            amat = np.ones([time_acc.sum(),2])
            amat[:,1] = altimetry['time'][time_acc] - altimetry['time'][time_acc].mean()
            tg_detrend = alttg_list['height'][region,:][time_acc] - np.matmul(amat,np.linalg.lstsq(amat, alttg_list['height'][region,:][time_acc],rcond=None)[0])
            # Accepted points
            distance  = gentools.point_grid_distance(alttg_list['coords'][region,0],alttg_list['coords'][region,1],altimetry['lat'],altimetry['lon'])
            distance[~altimetry['slm']] = 1e9
            alt_acc = np.array(np.where(distance < settings['max_dist'])).T
            corr_array = np.zeros(len(alt_acc))
            for alt in range(len(alt_acc)):
                alt_detrend = altimetry['ssh'][:,alt_acc[alt,0],alt_acc[alt,1]][time_acc] - np.matmul(amat, np.linalg.lstsq(amat, altimetry['ssh'][:,alt_acc[alt,0],alt_acc[alt,1]][time_acc], rcond=None)[0])
                corr_array[alt] = np.corrcoef(alt_detrend,tg_detrend)[0,1]
            corr_array[np.isnan(corr_array)]=-1
            if (corr_array>settings['min_corr']).sum()>0:
                # Compute weight
                corr_array_flt = corr_array[corr_array>settings['min_corr']]
                weight = corr_array_flt/corr_array_flt.sum()
                # Store data
                alttg_list['has_corr'][region] = True
                alttg_list['weight'][region] = weight
                alttg_list['time_acc'][region]   = time_acc
                alttg_list['alt_coords'][region] = alt_acc[corr_array>settings['min_corr'],:]
    return

def compute_ts():
    global altimetry, alttg_list, alttg_data, settings
    print('Sampling GIA and GRD at altimetry points...')
    # --------------------------------------------------------
    # Compute time series of VLM and residual VLM
    # vlm_res = Alt - GSL_gia - GSL_pd - TG + RSL_gia + RSL_pd
    # --------------------------------------------------------
    alttg_data = {}
    alttg_data['vlm_ts']    = mp_filled_float(np.zeros([len(alttg_list['id']),len(altimetry['time'])])*np.nan)
    alttg_data['resvlm_ts'] = mp_filled_float(np.zeros([settings['num_ens'],len(alttg_list['id']),len(altimetry['time'])])*np.nan)
    # Full VLM time series
    for region in range(len(alttg_list['has_corr'])):
        if alttg_list['has_corr'][region]:
            # vlm_ts_lcl = altimetry - tg: Weighted average of all grid points
            vlm_ts_lcl = (alttg_list['weight'][region] * (altimetry['ssh'][:,alttg_list['alt_coords'][region][:,0], alttg_list['alt_coords'][region][:,1]] - alttg_list['height'][region,:][:, np.newaxis])).sum(axis=1)
            alttg_data['vlm_ts'][region] = vlm_ts_lcl - np.nanmean(vlm_ts_lcl)
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(resvlm_ts_ens, range(settings['num_ens']))
    return

def resvlm_ts_ens(ens):
    print(ens)
    global altimetry, alttg_list, alttg_data, settings
    time_alt_idx = np.in1d(settings['years'],altimetry['time'])
    GRD = read_GRD_ens(ens, time_alt_idx, settings)
    GIA = read_GIA_ens(ens,settings)
    # Dynamic altimetry: altimetry - GSL_GRD - GSL_GIA
    # Dynamic tide gauge = tide gauge - RSL_GIA - RSL_GRD
    altimetry_dynamic = altimetry['ssh'] - GIA['gsl'][np.newaxis, :, :] * (altimetry['time'] - altimetry['time'].mean())[:, np.newaxis, np.newaxis] - GRD['gsl']
    for region in range(len(alttg_list['has_corr'])):
        if alttg_list['has_corr'][region]:
            tidegauge_dynamic = alttg_list['height'][region] - GIA['rsl'][alttg_list['tg_coords'][region,0],alttg_list['tg_coords'][region,1]]*(altimetry['time']-altimetry['time'].mean())-GRD['rsl'][:,alttg_list['tg_coords'][region, 0],alttg_list['tg_coords'][region, 1]]
            residual_vlm_lcl = (alttg_list['weight'][region]*(altimetry_dynamic[:,alttg_list['alt_coords'][region][:,0], alttg_list['alt_coords'][region][:, 1]]-tidegauge_dynamic[:, np.newaxis])).sum(axis=1)
            alttg_data['resvlm_ts'][ens,region,:] = residual_vlm_lcl
    return

def compute_trend():
    global alttg_list, alttg_data, settings
    print('Computing ALTTG trends...')
    alttg_data['vlm_trend']         = mp_filled_float(np.zeros([len(alttg_list['id']),2])*np.nan)
    alttg_data['resvlm_trend_mean'] = mp_filled_float(np.zeros([len(alttg_list['id']),2])*np.nan)
    alttg_data['resvlm_trend_ens']  = mp_filled_float(np.zeros([settings['num_ens'],len(alttg_list['id'])])*np.nan)

    alttg_data['resvlm_sterr_AR1'] = mp_filled_float(np.zeros(len(alttg_list['id']))*np.nan)
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(compute_trend_indiv, range(len(alttg_list['id'])))
    return

def compute_trend_indiv(region):
    global altimetry, alttg_list, alttg_data, settings
    if alttg_list['has_corr'][region]:
        print('   Region ' + str(region))
        # Trend in VLM
        alttg_data['vlm_trend'][region,:] = np.array(trend_ar1(region, altimetry['time'][alttg_list['time_acc'][region]], alttg_data['vlm_ts'][region][alttg_list['time_acc'][region]]))
        # Trend in residual VLM
        # AR1 trend uncertainty
        alttg_data['resvlm_trend_mean'][region,:] = np.array(trend_ar1(region, altimetry['time'][alttg_list['time_acc'][region]],  (settings['probability'][:,np.newaxis] * alttg_data['resvlm_ts'][:,region, :]).sum(axis=0)[alttg_list['time_acc'][region]]))
        alttg_data['resvlm_sterr_AR1'][region] = alttg_data['resvlm_trend_mean'][region,1].copy()
        amat = np.ones([alttg_list['time_acc'][region].sum(), 2])
        amat[:, 1] = altimetry['time'][alttg_list['time_acc'][region]] - altimetry['time'][alttg_list['time_acc'][region]].mean()
        resvlm_ens_mean = np.zeros(settings['num_ens'])
        for ens in range(settings['num_ens']):
            resvlm_ens_mean[ens] = np.linalg.lstsq(amat, alttg_data['resvlm_ts'][ens,region, alttg_list['time_acc'][region]], rcond=None)[0][1]
        alttg_data['resvlm_trend_ens'][:,region] = resvlm_ens_mean
        resvlm_mn = (settings['probability'] * resvlm_ens_mean).sum()
        resvlm_se = np.sqrt((settings['probability'] * (resvlm_ens_mean-resvlm_mn)**2).sum())
        alttg_data['resvlm_trend_mean'][region,1] = np.sqrt(alttg_data['resvlm_sterr_AR1'][region]**2+resvlm_se**2)
    return

def save_data():
    print('Saving data...')
    global alttg_list, alttg_data, settings
    acc_idx = np.isfinite(alttg_data['vlm_trend'][:,0])
    alttg = {}
    alttg['id'] = alttg_list['id'][acc_idx]
    alttg['coords'] = alttg_list['coords'][acc_idx,:]
    alttg['vlm_trend'] = alttg_data['vlm_trend'][acc_idx,:]
    alttg['resvlm_trend_mean'] = alttg_data['resvlm_trend_mean'][acc_idx,:]
    alttg['resvlm_trend_ens']  = alttg_data['resvlm_trend_ens'][:,acc_idx]
    alttg['resvlm_sterr_AR1']  = alttg_data['resvlm_sterr_AR1'][acc_idx]
    alttg['code'] = np.zeros(acc_idx.sum(),dtype=object)
    alttg['code'][:] = 'ALTTG'
    np.save(settings['fn_alttg_data'],alttg)
    return

### HELPER FUNCTIONS
# def read_GIA(settings):
#     print('Reading GIA...')
#     global GIA
#     GIA = {}
#     file_handle = Dataset(settings['fn_gia_rad'], 'r')
#     file_handle.set_auto_mask(False)
#     GIA['lat'] =  mp_filled_float(file_handle.variables['y'][:])
#     GIA['lon'] =  mp_filled_float(file_handle.variables['x'][:])
#     GIA['probability'] = mp_filled_float(file_handle.variables['probability'][:settings['num_ens']])
#     GIA['probability'] = GIA['probability']/GIA['probability'].sum()
#     GIA['rad'] = mp_filled_float(file_handle.variables['rad'][:settings['num_ens'],:,:])
#     file_handle.close()
#
#     file_handle = Dataset(settings['fn_gia_rsl'], 'r')
#     file_handle.set_auto_mask(False)
#     GIA['rsl'] = mp_filled_float(file_handle.variables['rsl'][:settings['num_ens'],:,:])
#     file_handle.close()
#     GIA['gsl'] = mp_filled_float(GIA['rad'] + GIA['rsl'])
#     return

def read_GIA_ens(ens,settings):
    GIA = {}
    # radial deformation
    file_handle = Dataset(settings['fn_gia_rad'], 'r')
    file_handle.set_auto_mask(False)
    if settings['test_run_ICE6G_D']:
        GIA['rad'] = file_handle.variables['rad'][:]
    else:
        GIA['rad'] = file_handle.variables['rad'][ens,:,:]
    # rsl
    file_handle = Dataset(settings['fn_gia_rsl'], 'r')
    file_handle.set_auto_mask(False)
    if settings['test_run_ICE6G_D']:
        GIA['rsl'] = file_handle.variables['RSL'][:]
    else:
        GIA['rsl'] = file_handle.variables['rsl'][ens,:,:]
    file_handle.close()
    #gsl
    GIA['gsl'] = mp_filled_float(GIA['rad'] + GIA['rsl'])
    return(GIA)

def read_GRD_ens(ens,time_alt_idx,settings):
    file_handle = Dataset(settings['dir_grd']+'grd_'+str(ens)+'.nc', 'r')
    file_handle.set_auto_mask(False)
    PD = {}
    PD['rad'] = file_handle.variables['rad'][time_alt_idx,:,:]
    PD['rsl'] = file_handle.variables['rsl'][time_alt_idx,:,:]
    PD['gsl'] = PD['rad'] + PD['rsl']
    file_handle.close()
    return(PD)

def read_altimetry():
    print('Reading altimetry...')
    global altimetry, settings
    altimetry = {}
    file_handle = Dataset(settings['fn_altimetry'], 'r')
    file_handle.set_auto_mask(False)
    altimetry['lat']  =  mp_filled_float(file_handle.variables['y'][:])
    altimetry['lon']  =  mp_filled_float(file_handle.variables['x'][:])
    altimetry['time'] =  mp_filled_float(file_handle.variables['t'][:])
    altimetry['ssh'] = mp_filled_float(file_handle.variables['z'][:])
    file_handle.close()
    altimetry['slm'] = mp_filled_bool(np.isfinite(altimetry['ssh'][-1,:,:]))
    return

def trend_ar1(region,time,tseries):
    global settings
    # Determine trend and associated uncertainty (AR1) using Hector
    # 1. Save settings
    dir_lcl = settings['dir_scratch']+str(region)+'/'
    if not os.path.isdir(dir_lcl):
        os.mkdir(dir_lcl)
    config_list = []
    config_list.append('DataFile input.mom\n')
    config_list.append('DataDirectory ./\n')
    config_list.append('OutputFile trend.out\n')
    config_list.append('interpolate no\n')
    config_list.append('firstdifference no\n')
    config_list.append('PhysicalUnit m\n')
    config_list.append('DegreePolynomial 1\n')
    config_list.append('seasonalsignal no\n')
    config_list.append('halfseasonalsignal no\n')
    config_list.append('estimateoffsets no\n')
    config_list.append('NoiseModels ARMA White\n')
    config_list.append('AR_p 1\n')
    config_list.append('MA_q 0\n')
    config_list.append('RandomiseFirstGuess yes\n')
    open(dir_lcl+'estimatetrend.ctl','w+').writelines(config_list)
    out = shutil.copy2(settings['fn_hector'],dir_lcl) # Copy Hector executable to scratch
    tseries = tseries - tseries.mean()
    mjd = 365.25 * (time-time[0])
    sample_period = np.min(np.diff(mjd))
    headerline = 'sampling period '+str(sample_period)
    np.savetxt(dir_lcl+'input.mom', np.transpose([mjd, tseries]), fmt=['%.4f', '%.4f'], header=headerline)
    os.chdir(dir_lcl)
    os.system(dir_lcl + 'estimatetrend >' + dir_lcl + 'output_orig.txt')
    output_data = [line.rstrip('\n') for line in open(dir_lcl + 'output_orig.txt')]
    trend_str = next((s for s in output_data if 'trend: ' in s), None)
    trend_mean = float(trend_str.split()[1])
    trend_sterr = float(trend_str.split()[3])
    os.chdir(os.getenv('HOME')+'/Scripts/Python/')
    return(trend_mean,trend_sterr)

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