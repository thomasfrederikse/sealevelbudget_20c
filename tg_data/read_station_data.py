# ------------------------------------------------------
# Create list with all station data for quick processing
# ------------------------------------------------------
import numpy as np
import os
from netCDF4 import Dataset
import mod_gentools as gentools
import scipy.io as scio
import statsmodels.api as sm
lowess = sm.nonparametric.lowess

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_psmsl'] = settings['dir_data'] + 'TideGauges/rlr_annual/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_hogarth']   = settings['dir_data']+'TideGauges/Hogarth/A_ext_Ann_msl_2018.mat'
    settings['fn_fl_met']   = settings['dir_data']+'TideGauges/met_monthly/filelist.txt'
    settings['fn_fl_rlr']   = settings['dir_data']+'TideGauges/rlr_annual/filelist.txt'

    settings['fn_station_data'] = settings['dir_budget']+'tg/station_data.npy'
    settings['fn_nodal']   = settings['dir_data']+'Nodal/Nodal.nc'
    settings['fn_ERA']     = settings['dir_budget']+'tg/ERA.nc'
    settings['fn_ERA_slm'] = settings['dir_budget']+'tg/ERA_slm.nc'

    settings['years'] = np.arange(1900,2019)
    settings['min_years_indiv'] = 10

    station_data = {}
    station_data['id']     = []
    station_data['name']   = []
    station_data['coords'] = []
    station_data['height'] = []

    station_data = read_annual_RLR_data(station_data,settings)
    station_data = read_non_RLR_data(station_data, settings)
    station_data = remove_meteo_forcing_nodal_cycle(station_data, settings)
    save_data(station_data, settings)
    return

def read_annual_RLR_data(station_data,settings):
    print('Reading PSMSL and Pether Hogarths data:')
    # -------------------------------------------
    # Read annual RLR data from PSMSL and Hogarth
    # If Hogarth differs from PSMSL use Hogarth
    # Only stations > 20 yr unflagged data
    # Store results in list
    # -------------------------------------------

    # List of blocked coastline codes
    blocked_coast_id_baltic = list(range(50,111))
    blocked_coast_id_russian_arctic = [30]
    blocked_coast_id_blacksea = list(range(295,306))
    blocked_coast_id_antarctica = [999]
    blocked_coast_id_greece = [290]
    blocked_coast_id_japan = list(range(641,649))
    blocked_coast_id = blocked_coast_id_baltic + blocked_coast_id_russian_arctic + blocked_coast_id_blacksea + blocked_coast_id_antarctica + blocked_coast_id_greece + blocked_coast_id_japan

    statlist_met = read_fl_psmsl(settings['fn_fl_met'], settings)
    statlist_rlr = read_fl_psmsl(settings['fn_fl_rlr'], settings)
    obs_hogarth  = read_hogarth(settings)

    for idx,id in enumerate(obs_hogarth['id']):
        # Station data
        statname = statlist_met['name'][statlist_met['id'] == id][0]
        coords   = np.array([statlist_met['lat'][statlist_met['id'] == id][0],statlist_met['lon'][statlist_met['id'] == id][0]])
        in_gulf_lawrence =  ((coords[0] > 46.22) & (coords[0] < 50.30) & (coords[1] > 287.35) & (coords[1] < 294.0))
        if ((statlist_met['coast'][statlist_met['id'] == id][0] not in blocked_coast_id) and (not in_gulf_lawrence)) or (id in [130,132,133,134]):
            # Station is not coast-blocked
            rlr_idx = (statlist_rlr['id'] == id)
            if rlr_idx.sum() == 0: # Station is not in RLR list, use Hogarth
                height = obs_hogarth['height'][idx,:]
            else: # PSMSL also has the data
                rlr_tseries = read_rlr_tseries(id, settings)
                # Use data from Peter Hogarth, and merge if PSMSL has been recently updated over the last 10 years
                if np.isfinite(rlr_tseries['height'][-10:]).sum() > np.isfinite(obs_hogarth['height'][idx,-10:]).sum():
                    obs_merge    = obs_hogarth['height'][idx,:]
                    both_end_idx = np.isfinite(rlr_tseries['height'][-10:]) & np.isfinite(obs_hogarth['height'][idx,-10:])
                    obs_merge[-10:] = rlr_tseries['height'][-10:] - rlr_tseries['height'][-10:][both_end_idx].mean() + obs_hogarth['height'][idx,-10:][both_end_idx].mean()
                    height = obs_merge - np.nanmean(obs_merge)
                else:
                    height = obs_hogarth['height'][idx, :]
                height[rlr_tseries['flag']] = np.nan
            height_detrend = gentools.detrend(settings['years'],height)
            height[np.abs(height_detrend)>3*np.nanstd(height_detrend)] = np.nan
            height -= np.nanmean(height)
            station_data['id'].append(id)
            station_data['name'].append(statname)
            station_data['coords'].append(coords)
            station_data['height'].append(height)
    return(station_data)

def read_non_RLR_data(station_data,settings):
    # Add additional stations to data set beyond PSMSL RLR annual data
    station_data = process_falklands(station_data,settings)
    station_data = process_dakar(station_data,settings)
    return(station_data)

def remove_meteo_forcing_nodal_cycle(station_data,settings):
    print('Determine and remove meteo and nodal cycle:')
    # -----------------------------------------------------------------------
    # Remove nodal cycle and meteo forcing from record
    # 1. Remove self-consistent equilibrium RSL nodal cycle from each record
    # 2. Remove local wind and pressure effects using least squares estimator
    # -----------------------------------------------------------------------

    station_data['id']     = np.array(station_data['id'])
    station_data['name']   = np.array(station_data['name'])
    station_data['coords'] = np.array(station_data['coords'])
    station_data['height'] = np.array(station_data['height'])

    # Nodal cycle
    station_data['height_nodal'] = np.zeros(station_data['height'].shape)
    station_data['height_meteo'] = np.zeros(station_data['height'].shape)
    station_data['height_corr']  = np.zeros(station_data['height'].shape)

    print('   Nodal cycle...')
    nodal = read_nodal(settings)
    for region in range(len(station_data['id'])):
        lat_idx = np.argmin(np.abs(nodal['lat']-station_data['coords'][region,0]))
        lon_idx = np.argmin(np.abs(nodal['lon']-station_data['coords'][region,1]))
        station_data['height_nodal'][region,:] = -nodal['amplitude'][lat_idx,lon_idx]*np.cos(2*np.pi/18.612958*(settings['years']-1922.7))
        station_data['height'][region,:] = station_data['height'][region,:] - station_data['height_nodal'][region,:]
    # Remove meteo forcing
    print('   Meteo...')
    meteo = read_meteo_data(settings)
    for region in range(len(station_data['id'])):
        # Determine points within 250 km of stations
        distance  = gentools.point_grid_distance(station_data['coords'][region,0],station_data['coords'][region,1],meteo['lat'],meteo['lon'])
        distance[~meteo['slm']] = 1e9
        meteo_acc = np.where(distance < 500000)
        time_acc  = np.isfinite(station_data['height'][region,:])

        # Remove linear trend from sea level
        amat = np.ones([time_acc.sum(), 2])
        amat[:, 1] = settings['years'][time_acc] - settings['years'][time_acc].mean()
        sl_detrend = station_data['height'][region,time_acc] - np.matmul(amat, np.linalg.lstsq(amat, station_data['height'][region,time_acc],rcond=None)[0])

        corr_mslp = np.zeros(len(meteo_acc[0]))
        corr_uws = np.zeros(len(meteo_acc[0]))
        corr_vws = np.zeros(len(meteo_acc[0]))

        # Find time series with highest correlation
        for mp in range(len(meteo_acc[0])):
            corr_mslp[mp] = np.corrcoef(sl_detrend,meteo['mslp'][time_acc,meteo_acc[0][mp],meteo_acc[1][mp]])[0, 1]
            corr_uws[mp] = np.corrcoef(sl_detrend,meteo['uws'][time_acc,meteo_acc[0][mp],meteo_acc[1][mp]])[0, 1]
            corr_vws[mp] = np.corrcoef(sl_detrend,meteo['vws'][time_acc,meteo_acc[0][mp],meteo_acc[1][mp]])[0, 1]
        slp_lcl = meteo['mslp'][time_acc, meteo_acc[0][np.argmax(np.abs(corr_mslp))], meteo_acc[1][np.argmax(np.abs(corr_mslp))]]
        uws_lcl = meteo['uws'][time_acc, meteo_acc[0][np.argmax(np.abs(corr_uws))], meteo_acc[1][np.argmax(np.abs(corr_uws))]]
        vws_lcl = meteo['vws'][time_acc, meteo_acc[0][np.argmax(np.abs(corr_vws))], meteo_acc[1][np.argmax(np.abs(corr_vws))]]

        # Remove meteo effects from final height estimate
        amat = np.ones([time_acc.sum(), 5])
        amat[:, 1] = settings['years'][time_acc] - settings['years'][time_acc].mean()
        amat[:, 2] = slp_lcl
        amat[:, 3] = uws_lcl
        amat[:, 4] = vws_lcl
        sol_tg = np.linalg.lstsq(amat,station_data['height'][region,time_acc],rcond=None)[0]

        # Compute covariance matrix
        residual_sq = np.sum((station_data['height'][region,time_acc]-np.matmul(amat,sol_tg))** 2)
        cov_mat = np.matrix(np.matmul(amat.T,amat)).I*(1/(time_acc.sum()-9))*residual_sq
        CI = np.array(1.65*np.sqrt(cov_mat.diagonal()))
        sol_tg[(np.abs(sol_tg) - CI).flatten() < 0] = 0
        sol_tg[1] = 0
        solution = np.matmul(amat, sol_tg)
        solution = solution - np.mean(solution)
        station_data['height_meteo'][region,:] = np.zeros(len(settings['years']))*np.nan
        station_data['height_corr'][region,:] = np.zeros(len(settings['years']))*np.nan
        station_data['height_meteo'][region,time_acc] = solution
        station_data['height_corr'][region,time_acc] = station_data['height'][region,time_acc] - solution
    return(station_data)

def save_data(station_data,settings):
    np.save(settings['fn_station_data'],station_data)
    return

#### HELPER FUNCTIONS #######
def process_falklands(station_data,settings):
    data_raw = np.loadtxt(settings['dir_data']+'Paleo/Falklands.txt',usecols=(1,5,6),skiprows=1)
    acc_year = data_raw[:,0]>1800
    data_raw = data_raw[acc_year,:]
    # Sort data
    data_raw = data_raw[np.argsort(data_raw[:,0])]
    # Unique data
    uniq_yrs = np.unique(data_raw[:,0]).astype(int)
    uniq_hgt = np.zeros(len(uniq_yrs))
    for i in range(len(uniq_yrs)):
        idx = (data_raw[:,0]==uniq_yrs[i])
        weight = (1./data_raw[idx,2]**2) /((1./data_raw[idx,2]**2).sum())
        uniq_hgt[i] = 1000*(weight * data_raw[idx,1]).sum()
    # Smooth record
    z = lowess(uniq_hgt, uniq_yrs, frac=0.45, it=3)
    hgt_falklands = np.interp(settings['years'],z[:,0],z[:,1],right=np.nan)
    station_data['id'].append(3000)
    station_data['name'].append('Falklands Salt Marsh')
    station_data['coords'].append(np.array([-51.6903, 302.1351]))
    station_data['height'].append(hgt_falklands)
    return(station_data)

def process_dakar(station_data,settings):
    dakar_raw = np.loadtxt(settings['dir_data']+'TideGauges/Dakar/Dakar_monthly_completedwithPSMSL_Thomas.dat')
    dakar_t = dakar_raw[:,0] + dakar_raw[:,1]/12 - 1/24
    dakar_h = dakar_raw[:,2]
    amat = np.ones([len(dakar_t),6])
    amat[:,1] = dakar_t - dakar_t.mean()
    amat[:,2] = np.sin(2*np.pi*dakar_t)
    amat[:,3] = np.cos(2*np.pi*dakar_t)
    amat[:,4] = np.sin(4*np.pi*dakar_t)
    amat[:,5] = np.cos(4*np.pi*dakar_t)
    sol = np.linalg.lstsq(amat,dakar_h,rcond=None)[0]
    sol[1] = 0
    dakar_noseas = dakar_h - np.matmul(amat,sol)
    dakar_annual = np.zeros(len(settings['years']))*np.nan
    for idx,yr in enumerate(settings['years']):
        acc_idx = (np.floor(dakar_t).astype(int)==yr)
        if acc_idx.sum()>6:
            dakar_annual[idx] = dakar_noseas[acc_idx].mean()
    station_data['id'].append(3001)
    station_data['name'].append('Dakar rescued')
    station_data['coords'].append(np.array([14.683333, 342.583333]))
    station_data['height'].append(dakar_annual)
    return(station_data)

def process_mpl_metric(station_data,settings):
    fn = settings['dir_data'] + 'TideGauges/met_monthly/data/177.metdata'
    tgdata = np.loadtxt(fn,delimiter=';')
    acc_idx = tgdata[:,1]>-1000
    tgdata = tgdata[acc_idx,:]
    amat = np.ones([len(tgdata[:,0]),6])
    amat[:,1] = tgdata[:,0] - tgdata[:,0].mean()
    amat[:,2] = np.sin(2*np.pi*tgdata[:,0])
    amat[:,3] = np.cos(2*np.pi*tgdata[:,0])
    amat[:,4] = np.sin(4*np.pi*tgdata[:,0])
    amat[:,5] = np.cos(4*np.pi*tgdata[:,0])
    sol = np.linalg.lstsq(amat,tgdata[:,1],rcond=None)[0]
    sol[1] = 0
    mpl_noseas = tgdata[:,1] - np.matmul(amat,sol)
    mpl_annual = np.zeros(len(settings['years']))*np.nan
    for idx,yr in enumerate(settings['years']):
        acc_idx = (np.floor(tgdata[:,0]).astype(int)==yr)
        if acc_idx.sum()>6:
            mpl_annual[idx] = mpl_noseas[acc_idx].mean()
    station_data['id'].append(3002)
    station_data['name'].append('Mar del Plata metric')
    station_data['coords'].append(np.array([-38.033333, 302.48333]))
    station_data['height'].append(mpl_annual)
    return(station_data)

def read_nodal(settings):
    nodal = {}
    file_handle = Dataset(settings['fn_nodal'], 'r')
    file_handle.set_auto_mask(False)
    nodal['lat'] =  file_handle.variables['y'][:]
    nodal['lon'] =  file_handle.variables['x'][:]
    nodal['amplitude'] = file_handle.variables['rsl_eq'][:]
    file_handle.close()
    return(nodal)

def read_meteo_data(settings):
    meteo = {}
    file_handle = Dataset(settings['fn_ERA'], 'r')
    file_handle.set_auto_mask(False)
    meteo['lat'] =  file_handle.variables['y'][:]
    meteo['lon'] =  file_handle.variables['x'][:]
    meteo['years'] =  file_handle.variables['t'][:]
    meteo['mslp'] =  file_handle.variables['mslp'][:]
    meteo['uws'] =  file_handle.variables['uws'][:]
    meteo['vws'] =  file_handle.variables['vws'][:]
    file_handle.close()
    lsm_raw = np.flipud(Dataset(settings['fn_ERA_slm'], 'r').variables['lsm'][:]._get_data().squeeze())
    meteo['slm'] = np.zeros(lsm_raw.shape,dtype=bool)
    # meteo['slm'][lsm_raw>0.6] = False
    meteo['slm'][lsm_raw<=0.8] = True
    return(meteo)

def read_fl_psmsl(fn,settings):
    station_list = {}
    filelist = np.loadtxt(fn,dtype=object,delimiter=';')
    station_list['id'] = filelist[:,0].astype(int)
    station_list['coast'] = filelist[:,4].astype(int)
    station_list['name'] = filelist[:,3]
    for stat in range(len(filelist)): station_list['name'][stat] = station_list['name'][stat].strip()
    station_list['lat']   = filelist[:,1].astype(float)
    station_list['lon']   = filelist[:,2].astype(float)
    station_list['lon'][station_list['lon']<0] = station_list['lon'][station_list['lon']<0] + 360
    station_list['flag'] = np.zeros(len(filelist),dtype=bool)
    for stat in range(len(filelist)):
        if filelist[stat,6]==' Y': station_list['flag'][stat] = True
    return(station_list)

def read_hogarth(settings):
    # Read data set from Peter Hogarth and select 1900-2018
    obs_hogarth = {}
    data_raw = scio.loadmat(settings['fn_hogarth'])
    h_hogarth = data_raw['A_ext_Ann_msl']
    t_hogarth = np.arange(1750,2019)

    # Interpolate station 157 using quadratic trend
    t_acc = np.isfinite(h_hogarth[:,156])
    amat = np.ones([t_acc.sum(),3])
    amat[:,1] = t_hogarth[t_acc] - t_hogarth[t_acc].mean()
    amat[:,2] = (t_hogarth[t_acc] - t_hogarth[t_acc].mean())**2
    sol = np.linalg.lstsq(amat,h_hogarth[t_acc,156],rcond=None)[0]

    amat = np.ones([len(t_acc),3])
    amat[:,1] = t_hogarth - t_hogarth[t_acc].mean()
    amat[:,2] = (t_hogarth - t_hogarth[t_acc].mean())**2
    h_fit = np.matmul(amat,sol)

    t_int = np.isnan(h_hogarth[:,156])
    h_hogarth[t_int, 156] = h_fit[t_int]

    h_hogarth = h_hogarth[np.in1d(t_hogarth,settings['years']),:].T
    psmsl_id_hogarth = np.arange(1,len(h_hogarth)+1)
    n_points = np.sum(np.isfinite(h_hogarth),axis=1)
    acc_len = n_points>settings['min_years_indiv']
    obs_hogarth['id']     = psmsl_id_hogarth[acc_len]
    obs_hogarth['height'] = h_hogarth[acc_len,:]
    return(obs_hogarth)

def read_rlr_tseries(id,settings):
    rlr_tseries={}
    rlr_tseries['height'] = np.zeros(len(settings['years']))*np.nan
    rlr_tseries['flag'] = np.zeros(len(settings['years']),dtype=bool)
    rlr_raw = np.loadtxt(settings['dir_psmsl'] + 'data/' + str(id) + '.rlrdata', dtype=object, delimiter=';')
    rlr_t = rlr_raw[:, 0].astype(int)
    rlr_h = rlr_raw[:, 1].astype(int)
    rlr_f = rlr_raw[:, 3].astype(int)
    acc_idx = (rlr_t > 1899) & (rlr_t < 2019) & (rlr_h != -99999)
    rlr_t = rlr_t[acc_idx]
    rlr_h = rlr_h[acc_idx]
    rlr_f = rlr_f[acc_idx]
    rlr_f = (rlr_f == 1) | (rlr_f == 11)
    store_idx = np.in1d(settings['years'],rlr_t)
    rlr_tseries['height'][store_idx] = rlr_h
    rlr_tseries['flag'][store_idx] = rlr_f
    rlr_tseries['height_flagged'] = rlr_tseries['height']
    rlr_tseries['height_flagged'][rlr_tseries['flag']] = np.nan
    return(rlr_tseries)

