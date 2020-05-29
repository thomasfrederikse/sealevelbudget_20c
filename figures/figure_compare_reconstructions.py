### ---------------------------------------------
# Compare my reconstruction to other recent ones
### ---------------------------------------------
import os
import numpy as np
from netCDF4 import Dataset
import mod_hector as hector
import mod_gentools as gentools
def main():
    set_settings()
    read_gmsl_rec()
    process_gmsl()
    process_obs()
    save_data()
    return

def set_settings():
    global settings
    settings = {}
    settings['dir_data']   = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['dir_gmt'] = os.getenv('HOME') + '/Scripts/GMT/Papers/Budget_20c/Comp_rec/'
    settings['fn_obs_ensembles'] = settings['dir_budget'] + 'results/obs_basin_global_ens.npy'
    settings['fn_grd_ensembles'] = settings['dir_budget'] + 'results/grd_basin_global_ens.npy'
    settings['fn_steric_ensembles'] = settings['dir_budget'] + 'results/steric_basin_global_ens.npy'

    settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['fn_D2017'] = settings['dir_data'] + 'Reconstructions/D2017.csv'
    settings['fn_D2019'] = settings['dir_data'] + 'Reconstructions/D2019.txt'
    settings['fn_H2015'] = settings['dir_data'] + 'Reconstructions/hay2015.csv'
    settings['fn_CW2011'] = settings['dir_data'] + 'Reconstructions/CW2011.txt'
    settings['fn_J2014']  = settings['dir_data'] + 'Reconstructions/gslGPChange2014.txt'

    settings['years']     = np.arange(1900,2019)
    settings['num_ens']   = 5000
    settings['filt_width'] = 30
    settings['est_serial_corr'] = True
    settings['probability'] = Dataset(settings['fn_gia_rsl'],'r').variables['probability'][:settings['num_ens']]._get_data()
    settings['probability'] = settings['probability']/settings['probability'].sum()
    return

def read_gmsl_rec():
    # Read time series of all reconstructions
    # Regrid to annual data and store in common
    # format
    global settings, gmsl
    gmsl = {}
    gmsl['d2019']  = read_d2019_annual(settings)
    gmsl['d2017'] = {}
    gmsl['d2017']['tseries'] = np.loadtxt(settings['fn_D2017'],delimiter=';')
    gmsl['h2015'] = {}
    gmsl['h2015']['tseries'] = np.loadtxt(settings['fn_H2015'],delimiter=',')
    gmsl['cw2011'] = read_cw2011_annual(settings)
    gmsl['j2014'] = read_j2014_annual(settings)
    # Remove baseline 2000 - 2010
    for prod in gmsl:
        baseline_idx = (gmsl[prod]['tseries'][:,0]>1999)&(gmsl[prod]['tseries'][:,0]<2011)
        gmsl[prod]['tseries'][:, 1] -= gmsl[prod]['tseries'][baseline_idx,1].mean()
    return()

def read_j2014_annual(settings):
    j2014 = {}
    obs_CW2011_raw = np.loadtxt(settings['fn_J2014'])
    t_floor = np.floor(obs_CW2011_raw[:, 0])
    t_uniq = np.unique(t_floor)
    h_mean = np.zeros(len(t_uniq))
    h_sterr =np.zeros(len(t_uniq))
    for idx,yr in enumerate(t_uniq):
        acc_idx = (t_floor==yr)
        h_mean[idx] = obs_CW2011_raw[acc_idx,3].mean()
        h_sterr[idx] = obs_CW2011_raw[acc_idx,4].mean()
    j2014['tseries'] = np.zeros([len(t_uniq),3])
    j2014['tseries'][:,0] = t_uniq
    j2014['tseries'][:,1] = h_mean
    j2014['tseries'][:,2] = h_sterr
    return(j2014)


def read_cw2011_annual(settings):
    cw2011 = {}
    obs_CW2011_raw = np.loadtxt(settings['fn_CW2011'])
    t_floor = np.floor(obs_CW2011_raw[:, 0])
    t_uniq = np.unique(t_floor)
    h_mean = np.zeros(len(t_uniq))
    h_sterr =np.zeros(len(t_uniq))
    for idx,yr in enumerate(t_uniq):
        acc_idx = (t_floor==yr)
        h_mean[idx] = obs_CW2011_raw[acc_idx,1].mean()
        h_sterr[idx] = obs_CW2011_raw[acc_idx,2].mean()
    cw2011['tseries'] = np.zeros([len(t_uniq),3])
    cw2011['tseries'][:,0] = t_uniq
    cw2011['tseries'][:,1] = h_mean
    cw2011['tseries'][:,2] = h_sterr
    return(cw2011)

def read_d2019_annual(settings):
    d2019 = {}
    obs_D2019_raw = np.loadtxt(settings['fn_D2019'],skiprows=1)
    t_floor = np.floor(obs_D2019_raw[:, 0])
    t_uniq = np.unique(t_floor)
    h_mean = np.zeros(len(t_uniq))
    h_sterr =np.zeros(len(t_uniq))
    for idx,yr in enumerate(t_uniq):
        acc_idx = (t_floor==yr)
        h_mean[idx] = obs_D2019_raw[acc_idx,1].mean()
        h_sterr[idx] = obs_D2019_raw[acc_idx,2].mean()
    d2019['tseries'] = np.zeros([len(t_uniq),3])
    d2019['tseries'][:,0] = t_uniq
    d2019['tseries'][:,1] = h_mean
    d2019['tseries'][:,2] = h_sterr
    return(d2019)

def process_gmsl():
    global gmsl, settings
    # Compute linear and sliding trends in each reconstruction
    # Linear trend
    for prod in gmsl:
        idx_1902_2010 = (gmsl[prod]['tseries'][:,0] > 1901) & (gmsl[prod]['tseries'][:,0] < 2011)
        idx_1957_2010 = (gmsl[prod]['tseries'][:,0] > 1957) & (gmsl[prod]['tseries'][:,0] < 2011)
        idx_1993_2010 = (gmsl[prod]['tseries'][:,0] > 1993) & (gmsl[prod]['tseries'][:,0] < 2011)
        gmsl[prod]['trend_1902_2010'] = gentools.lsqtrend(gmsl[prod]['tseries'][idx_1902_2010, 0], gmsl[prod]['tseries'][idx_1902_2010, 1])
        gmsl[prod]['trend_1957_2010'] = gentools.lsqtrend(gmsl[prod]['tseries'][idx_1957_2010, 0], gmsl[prod]['tseries'][idx_1957_2010, 1])
        gmsl[prod]['trend_1993_2010'] = gentools.lsqtrend(gmsl[prod]['tseries'][idx_1993_2010, 0], gmsl[prod]['tseries'][idx_1993_2010, 1])
    # Sliding trend
    fhalf = np.int(settings['filt_width']/2)
    for prod in gmsl:
        trend_smooth = np.zeros(len(settings['years']))*np.nan
        for idx, yr in enumerate(gmsl[prod]['tseries'][:,0]):
            idx_start = np.max([0,idx-fhalf])
            idx_stop  = np.min([len(gmsl[prod]['tseries'][:,0]),idx+fhalf])
            if (idx_stop-idx_start) > settings['filt_width']*0.75:
                yr_acc=gmsl[prod]['tseries'][idx_start:idx_stop,0]
                amat = np.ones([len(yr_acc), 2])
                amat[:, 1] = yr_acc
                if yr>=settings['years'][0]:
                    acc = np.in1d(settings['years'],yr)
                    trend_smooth[acc] = np.linalg.lstsq(amat, gmsl[prod]['tseries'][idx_start:idx_stop,1], rcond=None)[0][1]
        gmsl[prod]['gliding_trend'] = trend_smooth
    return()


def process_obs():
    # Process my own reconstruction
    global settings, obs, budget
    obs_ensembles = np.load(settings['fn_obs_ensembles'], allow_pickle=True).all()
    obs = compute_stats_indiv(obs_ensembles['global']['tg_full'][:settings['num_ens']], remove_baseline=True)

    grd_ensembles = np.load(settings['fn_grd_ensembles'], allow_pickle=True).all()
    steric_ensembles = np.load(settings['fn_steric_ensembles'], allow_pickle=True).all()

    budget_ens = grd_ensembles['global']['total'][:settings['num_ens'],:] + steric_ensembles['global'][:settings['num_ens'],:]
    budget = compute_stats_indiv(budget_ens, remove_baseline=True)
    return

def compute_stats_indiv(ensemble, remove_baseline=False):
    # ---------------------------
    # Compute ensemble statistics
    # for time series ensembles:
    # - Time series
    # - Sliding trends
    # - Linear trends over
    #   * 1902-2010
    #   * 1957-2010
    #   * 1993-2010
    # ---------------------------
    global settings
    stats = {}
    if remove_baseline:
        baseline_idx = (settings['years']>1999)&(settings['years']<2011)
        ensemble -= ensemble[:, baseline_idx].mean(axis=1)[:, np.newaxis]  # Remove baseline
    # 1. Time series
    stats['tseries'] = np.zeros([len(settings['years']), 3])
    stats['tseries'][:, 1] = np.sum(settings['probability'][:, np.newaxis] * ensemble, axis=0)
    for t in range(len(settings['years'])):
        sort_idx = np.argsort(ensemble[:, t])
        sort_cdf = np.cumsum(settings['probability'][sort_idx])
        stats['tseries'][t, 0] = ensemble[sort_idx, t][np.argmin(np.abs(sort_cdf - 0.05))]
        stats['tseries'][t, 2] = ensemble[sort_idx, t][np.argmin(np.abs(sort_cdf - 0.95))]
    # 2. Sliding trends
    stats['sliding_trend'] = estimate_sliding_trend_from_ensemble(ensemble)
    # 3. Linear trends
    stats['trend_1902_2010'] = estimate_linear_trend_ensemble(1902, 2010, ensemble)
    stats['trend_1957_2010'] = estimate_linear_trend_ensemble(1957, 2010, ensemble)
    stats['trend_1993_2010'] = estimate_linear_trend_ensemble(1993, 2010, ensemble)
    return (stats)

# --------------------
# Statistics functions
# --------------------
def estimate_sliding_trend_from_ensemble(ensemble):
    # Estimate sliding trend from ensemble
    global settings
    time_acc_idx = np.isfinite(ensemble[0, :])
    ensemble = ensemble[:, time_acc_idx]
    time_lcl = settings['years'][time_acc_idx]
    fhalf = np.int(settings['filt_width'] / 2)
    trend_ens = np.zeros(ensemble.shape) * np.nan
    for idx, yr in enumerate(time_lcl):
        idx_start = np.max([0, idx - fhalf])
        idx_stop = np.min([len(time_lcl), idx + fhalf])
        yr_acc = time_lcl[idx_start:idx_stop]
        amat = np.ones([len(yr_acc), 2])
        amat[:, 1] = yr_acc
        for ens in range(len(ensemble)):
            trend_ens[ens, idx] = np.linalg.lstsq(amat, ensemble[ens, idx_start:idx_stop], rcond=None)[0][1]
    sliding_trend_lcl = np.zeros([len(time_lcl), 3]) * np.nan
    sliding_trend_lcl[:, 1] = np.sum(settings['probability'][:, np.newaxis] * trend_ens, axis=0)
    for t in range(len(time_lcl) - 1):
        sort_idx = np.argsort(trend_ens[:, t])
        sort_cdf = np.cumsum(settings['probability'][sort_idx])
        sliding_trend_lcl[t, 0] = trend_ens[sort_idx, t][np.argmin(np.abs(sort_cdf - 0.05))]
        sliding_trend_lcl[t, 2] = trend_ens[sort_idx, t][np.argmin(np.abs(sort_cdf - 0.95))]
    sliding_trend_lcl[:np.rint(settings['filt_width'] / 4).astype(int), :] = np.nan
    sliding_trend_lcl[-np.rint(settings['filt_width'] / 4).astype(int):, :] = np.nan
    sliding_trend_stats = np.zeros([len(settings['years']), 3]) * np.nan
    sliding_trend_stats[time_acc_idx, :] = sliding_trend_lcl
    return (sliding_trend_stats)

def estimate_linear_trend_ensemble(ystart, ystop, ensemble):
    # Estimate mean and standard error of ensemble trends over specified years
    global settings
    linear_trend_stats = np.zeros(3)
    trend_ens = np.zeros(settings['num_ens'])
    acc_idx = np.in1d(settings['years'], np.arange(ystart, ystop + 1))
    amat = np.ones([acc_idx.sum(), 2])
    amat[:, 1] = settings['years'][acc_idx] - settings['years'][acc_idx].mean()
    amat_T = amat.T
    amat_sq = np.linalg.inv(np.dot(amat_T, amat))
    for ens in range(settings['num_ens']): trend_ens[ens] = np.dot(amat_sq, np.dot(amat_T, ensemble[ens, acc_idx]))[1]
    # Estimate trend error due to autocorrelation
    if settings['est_serial_corr']:
        tseries_mean = np.sum(settings['probability'][:, np.newaxis] * ensemble, axis=0)
        trend_ens += np.random.normal(0, 1, settings['num_ens']) * hector.est_trend(settings['years'], tseries_mean, NoiseModel='ARMA', AR=1, est_acc=False)['trend'][1]
    linear_trend_stats[0] = (settings['probability'] * trend_ens).sum()
    sort_idx = np.argsort(trend_ens)
    sort_cdf = np.cumsum(settings['probability'][sort_idx])
    linear_trend_stats[1] = trend_ens[sort_idx][np.argmin(np.abs(sort_cdf - 0.05))] - linear_trend_stats[0]
    linear_trend_stats[2] = trend_ens[sort_idx][np.argmin(np.abs(sort_cdf - 0.95))] - linear_trend_stats[0]
    return (linear_trend_stats)

def save_data():
    global gmsl, obs, budget, settings
    # Save time series
    for idx, prod in enumerate(gmsl):
        acc_idx = np.isfinite(gmsl[prod]['tseries'][:, 0])
        time_ci = np.append(gmsl[prod]['tseries'][acc_idx,0], np.flipud(gmsl[prod]['tseries'][acc_idx,0]))
        save_array_ci  = np.array([time_ci, np.append(gmsl[prod]['tseries'][acc_idx,1]+1.65*gmsl[prod]['tseries'][acc_idx,2], np.flipud(gmsl[prod]['tseries'][acc_idx,1]-1.65*gmsl[prod]['tseries'][acc_idx,2]))]).T
        save_array_mean = np.array([gmsl[prod]['tseries'][acc_idx,0], gmsl[prod]['tseries'][acc_idx,1]]).T
        np.savetxt(settings['dir_gmt'] + prod + '_tseries_m.txt', save_array_mean, fmt='%4.3f;%4.3f')
        np.savetxt(settings['dir_gmt'] + prod + '_tseries_ci.txt', save_array_ci, fmt='%4.3f;%4.3f')
    gentools.gmt_save_tseries_ci(settings['dir_gmt'],'obs_tseries',settings['years'],obs['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'],'budget_tseries',settings['years'],budget['tseries'])

    # Save sliding trends
    for idx, prod in enumerate(gmsl):
        acc_idx = np.isfinite(gmsl[prod]['gliding_trend'])
        save_array_mean = np.array([settings['years'][acc_idx], gmsl[prod]['gliding_trend'][acc_idx]]).T
        np.savetxt(settings['dir_gmt'] + prod + '_sliding_m.txt', save_array_mean, fmt='%4.3f;%4.3f')
    gentools.gmt_save_tseries_ci(settings['dir_gmt'],'obs_sliding',settings['years'],obs['sliding_trend'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'],'budget_sliding',settings['years'],budget['sliding_trend'])

    # Save linear trends
    for idx, prod in enumerate(gmsl):
        gentools.gmt_save_trend(settings['dir_gmt'] + prod + '_trend_1902.txt', idx+2, np.ones(3)*gmsl[prod]['trend_1902_2010'])
        gentools.gmt_save_trend(settings['dir_gmt'] + prod + '_trend_1957.txt', idx+9, np.ones(3)*gmsl[prod]['trend_1957_2010'])
        gentools.gmt_save_trend(settings['dir_gmt'] + prod + '_trend_1993.txt', idx+16, np.ones(3)*gmsl[prod]['trend_1993_2010'])
    gentools.gmt_save_trend(settings['dir_gmt'] + 'obs_trend_1902.txt', 1, obs['trend_1902_2010'])
    gentools.gmt_save_trend(settings['dir_gmt'] + 'obs_trend_1957.txt', 8, obs['trend_1957_2010'])
    gentools.gmt_save_trend(settings['dir_gmt'] + 'obs_trend_1993.txt', 15, obs['trend_1993_2010'])
    gentools.gmt_save_trend(settings['dir_gmt'] + 'budget_trend_1902.txt', 0, budget['trend_1902_2010'])
    gentools.gmt_save_trend(settings['dir_gmt'] + 'budget_trend_1957.txt', 7, budget['trend_1957_2010'])
    gentools.gmt_save_trend(settings['dir_gmt'] + 'budget_trend_1993.txt', 14, budget['trend_1993_2010'])
    return

