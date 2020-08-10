# ------------------------------------------
# Compute basin-scale time series and trends
# and save for plotting with GMT
# ------------------------------------------
import os
import numpy as np
from netCDF4 import Dataset
import mod_hector as hector
import matplotlib.pyplot as plt

def main():
    set_settings()
    compute_stats()
    #save_data()
    return

def set_settings():
    global settings
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['test_run_ICE6G_D'] = False
    settings['num_ens'] = 5000
    settings['fn_steric_ensembles'] = settings['dir_budget'] + 'results/steric_basin_global_ens.npy'
    if settings['test_run_ICE6G_D']:
        settings['fn_gia_ensembles'] = settings['dir_budget'] + 'results/gia_basin_global_ens_ice6g.npy'
        settings['fn_grd_ensembles'] = settings['dir_budget'] + 'results/grd_basin_global_ens_ice6g.npy'
        settings['fn_obs_ensembles'] = settings['dir_budget'] + 'results/obs_basin_global_ens_ice6g.npy'
        settings['fn_alt_ensembles'] = settings['dir_budget'] + 'results/alt_basin_global_ens_ice6g.npy'
        settings['fn_save_stats'] = settings['dir_budget'] + 'results/stats_basin_global_ice6g.npy'
        settings['probability'] = np.ones(settings['num_ens'])/settings['num_ens']
    else:
        settings['fn_gia_ensembles'] = settings['dir_budget'] + 'results/gia_basin_global_ens.npy'
        settings['fn_grd_ensembles'] = settings['dir_budget'] + 'results/grd_basin_global_ens.npy'
        settings['fn_obs_ensembles'] = settings['dir_budget'] + 'results/obs_basin_global_ens.npy'
        settings['fn_alt_ensembles'] = settings['dir_budget'] + 'results/alt_basin_global_ens.npy'
        settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
        settings['fn_save_stats'] = settings['dir_budget'] + 'results/stats_basin_global.npy'
        settings['probability'] = Dataset(settings['fn_gia_rsl'], 'r').variables['probability'][:settings['num_ens']]._get_data()
        settings['probability'] = settings['probability'] / settings['probability'].sum()
    settings['years'] = np.arange(1900, 2019)
    settings['filt_width'] = 30  # Moving window when estimating running trends
    settings['est_serial_corr'] = True  # Estimate serial correlation in linear trends
    return

def save_data():
    print('Saving...')
    global stats, settings
    np.save(settings['fn_save_stats'], stats)
    return

def compute_stats():
    print('Computing statistics...')
    global stats, settings
    stats = {}
    stats['basin'] = np.zeros(6, dtype=object)
    stats['global'] = {}
    # Data
    obs_ensembles = np.load(settings['fn_obs_ensembles'], allow_pickle=True).all()
    steric_ensembles = np.load(settings['fn_steric_ensembles'], allow_pickle=True).all()
    gia_ensembles = np.load(settings['fn_gia_ensembles'], allow_pickle=True).all()
    grd_ensembles = np.load(settings['fn_grd_ensembles'], allow_pickle=True).all()
    alt_ensembles = np.load(settings['fn_alt_ensembles'], allow_pickle=True).all()

    time_gia = (settings['years'] - settings['years'][-17:].mean())[np.newaxis,:]
    for basin in range(len(stats['basin'])):
        print('   Basin ' + str(basin))
        stats['basin'][basin] = {}
        stats['basin'][basin]['obs'] = compute_stats_indiv(obs_ensembles['basin'][basin]['tg_full'][:settings['num_ens']],remove_baseline=True)
        stats['basin'][basin]['steric'] = compute_stats_indiv(steric_ensembles['basin'][basin][:settings['num_ens'], :],remove_baseline=True)
        for term in grd_ensembles['basin'][basin]:
            stats['basin'][basin]['grd_' + term] = compute_stats_indiv(grd_ensembles['basin'][basin][term][:settings['num_ens'], :],remove_baseline=True)
        stats['basin'][basin]['grd_ice'] = compute_stats_indiv((grd_ensembles['basin'][basin]['glac']+grd_ensembles['basin'][basin]['GrIS']+grd_ensembles['basin'][basin]['AIS'])[:settings['num_ens']], remove_baseline=True)

        stats['basin'][basin]['gia'] = compute_stats_indiv(gia_ensembles['basin'][basin][:settings['num_ens']][:,np.newaxis]*time_gia,remove_baseline=True)
        stats['basin'][basin]['altimetry'] = compute_stats_indiv(alt_ensembles['basin'][basin][:settings['num_ens'], :],remove_baseline=True)
        stats['basin'][basin]['budget'] = compute_stats_indiv(steric_ensembles['basin'][basin][:settings['num_ens'], :]+grd_ensembles['basin'][basin]['total'][:settings['num_ens'], :]+gia_ensembles['basin'][basin][:settings['num_ens']][:,np.newaxis]*time_gia,remove_baseline=True)
        stats['basin'][basin]['diff'] = compute_stats_indiv(obs_ensembles['basin'][basin]['tg_full'][:settings['num_ens']]-steric_ensembles['basin'][basin][:settings['num_ens'], :]-grd_ensembles['basin'][basin]['total'][:settings['num_ens'], :]-gia_ensembles['basin'][basin][:settings['num_ens']][:,np.newaxis]*time_gia,remove_baseline=True)
        stats['basin'][basin]['obs_steric'] = compute_stats_indiv(obs_ensembles['basin'][basin]['tg_full'][:settings['num_ens']]-steric_ensembles['basin'][basin][:settings['num_ens'], :]-gia_ensembles['basin'][basin][:settings['num_ens']][:,np.newaxis]*time_gia,remove_baseline=True)

    print('   Global')
    stats['global']['obs'] = compute_stats_indiv(obs_ensembles['global']['tg_full'][:settings['num_ens']],remove_baseline=True)
    stats['global']['steric'] = compute_stats_indiv(steric_ensembles['global'][:settings['num_ens'], :])
    for term in grd_ensembles['global']:
        stats['global']['grd_' + term] = compute_stats_indiv(grd_ensembles['global'][term][:settings['num_ens'], :])
    stats['global']['grd_ice'] = compute_stats_indiv((grd_ensembles['global']['glac'] + grd_ensembles['global']['GrIS'] + grd_ensembles['global']['AIS'])[:settings['num_ens']], remove_baseline=True)

    stats['global']['altimetry'] = compute_stats_indiv(alt_ensembles['global'][:settings['num_ens'], :])
    stats['global']['budget'] = compute_stats_indiv(steric_ensembles['global'][:settings['num_ens'], :]+grd_ensembles['global']['total'][:settings['num_ens'], :])
    stats['global']['diff'] = compute_stats_indiv(obs_ensembles['global']['tg_full'][:settings['num_ens']]-steric_ensembles['global'][:settings['num_ens'], :]-grd_ensembles['global']['total'][:settings['num_ens'], :])
    stats['global']['obs_steric'] = compute_stats_indiv(obs_ensembles['global']['tg_full'][:settings['num_ens']]-steric_ensembles['global'][:settings['num_ens'], :])

    return()


def compute_stats_indiv(ensemble,remove_baseline=False):
    # ---------------------------
    # Compute ensemble statistics
    # for time series ensembles:
    # - Time series
    # - Sliding trends
    # - Linear trends over
    #   * 1900-2018
    #   * 1955-2018
    #   * 1993-2018
    # ---------------------------
    global settings
    stats = {}

    if remove_baseline: ensemble -= ensemble[:, -17:].mean(axis=1)[:, np.newaxis]  # Remove baseline

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
    stats['trend_1900_2018'] = estimate_linear_trend_ensemble(1900, 2018, ensemble)
    stats['trend_1900_1990'] = estimate_linear_trend_ensemble(1900, 1990, ensemble)
    stats['trend_1957_2018'] = estimate_linear_trend_ensemble(1957, 2018, ensemble)
    stats['trend_1971_2018'] = estimate_linear_trend_ensemble(1971, 2018, ensemble)
    stats['trend_1993_2018'] = estimate_linear_trend_ensemble(1993, 2018, ensemble)
    stats['trend_2006_2018'] = estimate_linear_trend_ensemble(2006, 2018, ensemble)
    return(stats)

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
        amat_T = amat.T
        amat_sq = np.linalg.inv(np.dot(amat_T, amat))
        for ens in range(len(ensemble)):
            trend_ens[ens, idx] = np.dot(amat_sq, np.dot(amat_T, ensemble[ens, idx_start:idx_stop]))[1]
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
        trend_ens += np.random.normal(0, 1, settings['num_ens']) * hector.est_trend(settings['years'], tseries_mean, NoiseModel='GGM', AR=1, est_acc=False)['trend'][1]
    linear_trend_stats[0] = (settings['probability'] * trend_ens).sum()
    sort_idx = np.argsort(trend_ens)
    sort_cdf = np.cumsum(settings['probability'][sort_idx])
    linear_trend_stats[1] = trend_ens[sort_idx][np.argmin(np.abs(sort_cdf - 0.05))] - linear_trend_stats[0]
    linear_trend_stats[2] = trend_ens[sort_idx][np.argmin(np.abs(sort_cdf - 0.95))] - linear_trend_stats[0]
    return (linear_trend_stats)

def plot_basin_movav():
    global stats, settings
    plt.figure(figsize=(12, 9))
    for reg in range(6):
        plt.subplot(3,2,reg+1)
        plt.plot(settings['years'],stats['global']['obs']['sliding_trend'][:,1],linewidth=2,color='black')
        plt.plot(settings['years'],stats['basin'][reg]['obs']['sliding_trend'][:,1],linewidth=2,color='C0')
    return

if __name__ == '__main__':
    main()
