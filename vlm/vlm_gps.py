# ----------------------------------------------------
# Select all GPS stations that are
#   - 30 km from a tide-gauge region
#   - Have 4 years of data over 1985-2019.0
# Compute vlm and residual vlm for these stations
#   - Using original time series data
#   - Using time series data corrected for GIA and GRD
# Save data
#   - Trends for each GPS station and ensemble member
# ----------------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
import mod_midas_py_ens as midas_py_ens
import mod_midas_py_single as midas_py_single
from scipy.interpolate import interp1d
import multiprocessing as mp
import ctypes as ct

def main():
    set_settings()
    find_gps_stations()
    sample_GIA_GRD()
    compute_gps_trend()
    save_data()
    return

def set_settings():
    global settings
    settings = {}
    settings['test_run_ICE6G_D'] = True
    settings['dir_data']    = os.getenv('HOME') + '/Data/'
    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 4
    else:
        settings['nproc'] = 40
    settings['dir_budget']  = settings['dir_data'] + 'Budget_20c/'

    if settings['test_run_ICE6G_D']:
        settings['dir_grd'] = settings['dir_budget'] + 'grd_ICE6G/'
        settings['fn_gia_rad'] = settings['dir_data']+'GIA/ICE6G_D/ICE6G_D_05.nc'
        settings['fn_gps_data'] = settings['dir_budget'] + 'vlm/gps_data_ice6g.npy'
    else:
        settings['dir_grd'] = settings['dir_budget']+'grd/'
        settings['fn_gia_rad'] = settings['dir_data']+'GIA/ICE6G_D/Ensemble/rad_ens_05.nc'
        settings['fn_gps_data'] = settings['dir_budget'] + 'vlm/gps_data.npy'

    settings['fn_gps_tseries'] = settings['dir_budget'] + 'vlm/gps_tseries.npy'
    settings['fn_regions_for_selection'] = settings['dir_budget']+'tg/regions_for_selection.npy'
    settings['max_dist'] = 30000  # Maximum distance between GPS and TG
    settings['min_years'] = 4     # Minimum number of years for GPS trend
    settings['num_ens'] = 100    # Number of ensembles
    settings['time_grd'] = np.arange(1985,2019)
    return

def find_gps_stations():
    global gps_regions, settings
    # ---------------------------------------------------
    # Find GPS stations that:
    # 1. Are within 30 km of each region
    # 2. Are at least 4 years long between 1900-2018.9999
    # ---------------------------------------------------
    # Load data
    gps_tseries           = np.load(settings['fn_gps_tseries'],allow_pickle=True)
    regions_for_selection = np.load(settings['fn_regions_for_selection'], allow_pickle=True).all()

    # Make array with GPS corrdinates
    gps_coords = np.zeros([len(gps_tseries),2])
    for i in range(len(gps_tseries)):
        gps_coords[i,0] = gps_tseries[i]['lat']
        gps_coords[i,1] = gps_tseries[i]['lon']

    gps_regions = np.zeros(len(regions_for_selection['id']),dtype=object)
    for region in range(len(regions_for_selection['id'])):
        print(regions_for_selection['id'][region])
        # Distance matrix
        distance = 2*6371000*np.arcsin(np.sqrt(np.sin(np.deg2rad(0.5*(gps_coords[:,0].astype(float)-regions_for_selection['coords'][region,0])))**2+np.cos(np.deg2rad(regions_for_selection['coords'][region,0]))*np.cos(np.deg2rad(gps_coords[:,0].astype(float)))*np.sin(np.deg2rad(0.5*(gps_coords[:,1].astype(float)-regions_for_selection['coords'][region,1])))**2))
        distance_acc = (distance < settings['max_dist'])
        gps_regions[region] = []
        if distance_acc.sum()>0: # There are nearby GPS stations: check length
            acc_idx = np.where(distance_acc)[0]
            for i in range(len(acc_idx)):
                if (gps_tseries[acc_idx[i]]['time']<2019).sum()>(365*settings['min_years']):
                    gps_regions[region].append(gps_tseries[acc_idx[i]])
    return

def sample_GIA_GRD():
    global gps_regions, gia_grd, settings
    regions_for_selection = np.load(settings['fn_regions_for_selection'], allow_pickle=True).all()

    gia_grd = {}
    gia_grd['sample_points'] = mp_empty_int([len(gps_regions),2])
    gia_grd['gia'] = mp_empty_float([settings['num_ens'],len(gps_regions)])
    gia_grd['grd'] = mp_empty_float([settings['num_ens'],len(gps_regions),len(settings['time_grd'])])

    # load GIA
    GIA = read_GIA_rad(settings)
    # Compute sampling indices
    for region in range(len(gps_regions)):
        gia_grd['sample_points'][region,0] = np.argmin(np.abs(GIA['lat']-regions_for_selection['coords'][region][0]))
        gia_grd['sample_points'][region,1] = np.argmin(np.abs(GIA['lon']-regions_for_selection['coords'][region][1]))

    # Sample GIA
    for region in range(len(gps_regions)):
        gia_grd['gia'][:,region] =  GIA['rad'][:,gia_grd['sample_points'][region,0],gia_grd['sample_points'][region,1]]

    gia_grd['probability'] = mp_filled_float(GIA['probability'])
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(sample_grd_indiv, range(settings['num_ens']))
    return()

def sample_grd_indiv(ens):
    global gia_grd, settings
    print(ens)
    GRD_rad_ens = read_GRD_rad_ens(ens, settings)
    for region in range(gia_grd['sample_points'].shape[0]):
        gia_grd['grd'][ens, region, :] = GRD_rad_ens[:, gia_grd['sample_points'][region, 0], gia_grd['sample_points'][region, 1]]
    return

def compute_gps_trend():
    # ----------------------------------------------------------
    # Compute raw and residual GPS trends at tide gauge stations
    # - ensembles:
    # - VLM trends
    # - Residual VLM trends
    # ----------------------------------------------------------
    global gia_grd, gps_regions, gps_trends, settings
    # Data preparation for MIDAS
    gps_trends = {}
    gps_trends['vlm_trend_mean']          = []
    gps_trends['resvlm_trend_mean']          = []
    gps_trends['resvlm_trend_mean_ens']   = []
    gps_trends['resvlm_trend_sterr_ens']  = []

    for region in range(len(gps_regions)):
        if len(gps_regions[region]) > 0:
            gps_trends['vlm_trend_mean'].append(mp_empty_float([len(gps_regions[region]), 2]))
            gps_trends['resvlm_trend_mean'].append(mp_empty_float([len(gps_regions[region]), 2]))

            gps_trends['resvlm_trend_mean_ens'].append(mp_empty_float([len(gps_regions[region]), settings['num_ens']]))
            gps_trends['resvlm_trend_sterr_ens'].append(mp_empty_float([len(gps_regions[region]), settings['num_ens']]))
        else:
            gps_trends['vlm_trend_mean'].append([])
            gps_trends['resvlm_trend_mean'].append([])
            gps_trends['resvlm_trend_mean_ens'].append([])
            gps_trends['resvlm_trend_sterr_ens'].append([])
        # gps_trends['vlm'].append(mp_empty_float([len(gps_statinfo['code'][region]), 2]))
        # gps_trends['resvlm_trend_mean'].append(mp_empty_float([len(gps_statinfo['code'][region]), 2]))
        # gps_trends['resvlm_trend_ens'].append(mp_empty_float([len(gps_statinfo['code'][region]), settings['num_ens']]))
        # gps_trends['resvlm_trend_ens_sterr'].append(mp_empty_float([len(gps_statinfo['code'][region]), settings['num_ens']]))
    # Loop over all stations to compute the trends
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(compute_gps_trend_indiv, range(len(gps_regions)))
    return()

def compute_gps_trend_indiv(region):
    # ----------------------------------------------------------
    # Compute raw and residual GPS trends at tide gauge stations
    # Individual station function
    # ----------------------------------------------------------
    global gia_grd, gps_regions, gps_trends, settings
    print(region)
    if len(gps_regions[region])>0:
        for stat in range(len(gps_regions[region])):
            # Read station time, height
            time_acc = (gps_regions[region][stat]['time']<2019.0)
            time_len   = len(gps_regions[region][stat]['time'][time_acc])
            time_pad   = np.pad(gps_regions[region][stat]['time'][time_acc], (0, 9999 - time_len), 'constant', constant_values=0)
            height_pad = 1000*np.pad(gps_regions[region][stat]['height'][time_acc], (0, 9999 - time_len), 'constant', constant_values=0)
            # Find jumps
            if gps_regions[region][stat]['has_jumps']:
                njumps = len(gps_regions[region][stat]['jumps'].astype(float))
                jump_date = np.pad(gps_regions[region][stat]['jumps'].astype(float), (0, 100 - njumps), 'constant', constant_values=0)
            else:
                njumps = 0
                jump_date = np.zeros(100)
            gps_trends['vlm_trend_mean'][region][stat, :] = np.array(midas_py_single.midas_py(time_len, time_pad, height_pad, njumps, jump_date))

            # Create ensemble of perturbed solutions
            GIA_lcl = (gps_regions[region][stat]['time'][time_acc] - gps_regions[region][stat]['time'][time_acc].mean()) * gia_grd['gia'][:, region][:, np.newaxis]
            GRD_lcl = interp1d(settings['time_grd'] + 0.5, gia_grd['grd'][:, region, :], axis=1, kind='linear', fill_value='extrapolate')(gps_regions[region][stat]['time'][time_acc])
            resid_vlm_ens = np.zeros([5000, 9999])
            resid_vlm_ens[:settings['num_ens'], :time_len] = 1000*gps_regions[region][stat]['height'][time_acc] - GIA_lcl - GRD_lcl

            # Residual trends trends
            resid_vlm_trend, resid_vlm_sterr = midas_py_ens.midas_py(time_len, time_pad, resid_vlm_ens.T, settings['num_ens'], njumps, jump_date)
            resid_vlm_trend = resid_vlm_trend[:settings['num_ens']]
            resid_vlm_sterr = resid_vlm_sterr[:settings['num_ens']]

            # Generate ensemble
            gps_trends['resvlm_trend_mean_ens'][region][stat, :]  = resid_vlm_trend
            gps_trends['resvlm_trend_sterr_ens'][region][stat, :] = resid_vlm_sterr

            # Compute mean, sterr
            resid_vlm_trend_mean  = (gia_grd['probability']*resid_vlm_trend).sum()
            resid_vlm_trend_sterr = np.sqrt((gia_grd['probability']*(resid_vlm_trend-resid_vlm_trend_mean)**2).sum())
            resid_vlm_trend_sterr = np.sqrt(resid_vlm_trend_sterr**2 + resid_vlm_sterr.mean()**2)
            gps_trends['resvlm_trend_mean'][region][stat,:] = np.array([resid_vlm_trend_mean,resid_vlm_trend_sterr])
    return

# def read_jumplist(settings):
#     fn = settings['dir_gps'] + 'steps.txt'
#     jumplist = np.loadtxt(fn, usecols=(0,1,2),dtype=object)
#     jumplist[:,2] = jumplist[:,2].astype(int)
#     jumplist = jumplist[jumplist[:,2]==1]
#     # Convert date
#     decyr_lookup = np.loadtxt(settings['dir_gps'] + 'decyr.txt', delimiter=' ', dtype=object, usecols=(0,1))
#     decyr_dict = dict(zip(decyr_lookup[:,0],decyr_lookup[:,1].astype(float)))
#     for i in range(len(jumplist)):
#         jumplist[i,1] = decyr_dict[jumplist[i,1]]
#     return(jumplist)
#
# def read_gps_time(statcode,settings):
#     fn = settings['dir_gps']+statcode+'.tenv'
#     gps_time=np.loadtxt(fn,usecols=2,skiprows=1)
#     return(gps_time)
#
# def read_gps_tseries(statcode,settings):
#     fn = settings['dir_gps']+statcode+'.tenv'
#     gps_tseries=np.loadtxt(fn,usecols=(2,12),skiprows=1)
#     gps_tseries[:,1] = gps_tseries[:,1]*1000
#     return(gps_tseries)

def read_GIA_rad(settings):
    GIA = {}
    file_handle = Dataset(settings['fn_gia_rad'], 'r')
    file_handle.set_auto_mask(False)
    GIA['lat'] =  file_handle.variables['y'][:]
    GIA['lon'] =  file_handle.variables['x'][:]
    if settings['test_run_ICE6G_D']:
        GIA['probability'] = np.ones(settings['num_ens'])/settings['num_ens']
        GIA['rad'] = file_handle.variables['rad'][:,:][np.newaxis,:,:]
    else:
        GIA['probability'] = file_handle.variables['probability'][:settings['num_ens']]
        GIA['probability'] = GIA['probability']/GIA['probability'].sum()
        GIA['rad'] = file_handle.variables['rad'][:settings['num_ens'],:,:]
    file_handle.close()
    return(GIA)

def read_GRD_rad_ens(ens,settings):
    file_handle = Dataset(settings['dir_grd']+'grd_'+str(ens)+'.nc', 'r')
    file_handle.set_auto_mask(False)
    GRD_rad_ens = file_handle.variables['rad'][85:,:,:]
    file_handle.close()
    return(GRD_rad_ens)

def save_data():
    global gps_regions, gps_trends, settings
    regions_for_selection = np.load(settings['fn_regions_for_selection'], allow_pickle=True).all()
    gps_data = np.zeros(len(gps_trends['vlm_trend_mean']),dtype=object)
    for region in range(len(gps_data)):
        gps_data[region] = {}
        gps_data[region]['id']     = regions_for_selection['id'][region]
        gps_data[region]['coords'] = regions_for_selection['coords'][region]
        if len(gps_trends['vlm_trend_mean'][region]) == 0:
            gps_data[region]['has_gps'] = False
        else:
            gps_data[region]['has_gps']  = True
            gps_data[region]['gps_name'] = []
            gps_data[region]['ref_frame'] = []
            for stat in gps_regions[region]:
                gps_data[region]['gps_name'].append(stat['name'])
                gps_data[region]['ref_frame'].append(stat['ref_frame'])
            gps_data[region]['vlm_trend_mean']    = np.array(gps_trends['vlm_trend_mean'][region])
            gps_data[region]['resvlm_trend_mean'] = np.array(gps_trends['resvlm_trend_mean'][region])
            gps_data[region]['resvlm_mean_ens']   = gps_trends['resvlm_trend_mean_ens'][region]
            gps_data[region]['resvlm_sterr_ens']  = gps_trends['resvlm_trend_sterr_ens'][region]
    # # tg_tst = []
    # # for region in range(len(gps_data)):
    # #     if gps_data[region]['has_gps']:
    # #         print(gps_data[region]['id'],gps_trends['resvlm_trend_mean'][region][:,0].mean())
    # #         tg_tst.append(gps_trends['resvlm_trend_mean'][region][:, 0].mean())
    #
    # gps_data = {}
    # gps_data['id']     = regions_for_selection['id']
    # gps_data['coords'] = regions_for_selection['coords']
    #
    # gps_data['code'] = gps_statinfo['code']
    #
    # gps_data['vlm_trend'] = np.array(gps_trends['vlm'])
    # gps_data['resvlm_trend_mean'] = np.array(gps_trends['resvlm_trend_mean'])
    # gps_data['resvlm_trend_ens']  = np.array(gps_trends['resvlm_trend_ens'])
    # gps_data['resvlm_trend_ens_sterr']  = np.array(gps_trends['resvlm_trend_ens_sterr'])
    np.save(settings['fn_gps_data'],gps_data)
    return


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
