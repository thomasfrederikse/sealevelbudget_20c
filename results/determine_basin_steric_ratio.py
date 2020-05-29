# Compute the ratio of steric trends between basins
import os
import numpy as np
import mod_hector as hector
from netCDF4 import Dataset
def main():
    global settings
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['dir_gmt'] = os.getenv('HOME') + '/Scripts/GMT/Papers/Budget_20c/Mass_steric_ratio/'
    settings['fn_steric_ensembles'] = settings['dir_budget'] + 'results/steric_basin_global_ens.npy'
    settings['fn_grd_ensembles'] = settings['dir_budget'] + 'results/grd_basin_global_ens.npy'
    settings['fn_obs_ensembles'] = settings['dir_budget'] + 'results/obs_basin_global_ens.npy'

    settings['est_serial_corr'] = True
    settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['years']     = np.arange(1900,2019)
    settings['num_ens']   = 5000
    settings['filt_width'] = 30
    settings['probability'] = Dataset(settings['fn_gia_rsl'], 'r').variables['probability'][:settings['num_ens']]._get_data()
    settings['probability'] = settings['probability'] / settings['probability'].sum()


    steric_ensembles = np.load(settings['fn_steric_ensembles'], allow_pickle=True).all()
    grd_ensembles = np.load(settings['fn_grd_ensembles'], allow_pickle=True).all()

    sa_vs_ep = estimate_linear_trend_ensemble(1957, 2018, steric_ensembles['basin'][2]-steric_ensembles['basin'][3])
    sa_vs_ep_ratio = estimate_linear_trend_ratio_ensemble(1957, 2018, steric_ensembles['basin'][2], steric_ensembles['basin'][3])
    steric_mass_ratio = estimate_linear_trend_ratio_ensemble(1900, 2018, grd_ensembles['global']['total'],steric_ensembles['global'])

    return

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

def estimate_linear_trend_ratio_ensemble(ystart, ystop, ensemble1,ensemble2):
    # Estimate mean and standard error of ensemble trends over specified years
    global settings

    trend_ens1 = np.zeros(settings['num_ens'])
    trend_ens2 = np.zeros(settings['num_ens'])

    acc_idx = np.in1d(settings['years'], np.arange(ystart, ystop + 1))

    amat = np.ones([acc_idx.sum(), 2])
    amat[:, 1] = settings['years'][acc_idx] - settings['years'][acc_idx].mean()
    amat_T = amat.T
    amat_sq = np.linalg.inv(np.dot(amat_T, amat))

    for ens in range(settings['num_ens']):
        trend_ens1[ens] = np.dot(amat_sq, np.dot(amat_T, ensemble1[ens, acc_idx]))[1]
        trend_ens2[ens] = np.dot(amat_sq, np.dot(amat_T, ensemble2[ens, acc_idx]))[1]

    ratio = trend_ens1/trend_ens2
    ratio_stats = np.zeros(3)
    ratio_stats[0] = (settings['probability'] * ratio).sum()
    sort_idx = np.argsort(ratio)
    sort_cdf = np.cumsum(settings['probability'][sort_idx])
    ratio_stats[1] = ratio[sort_idx][np.argmin(np.abs(sort_cdf - 0.05))] - ratio_stats[0]
    ratio_stats[2] = ratio[sort_idx][np.argmin(np.abs(sort_cdf - 0.95))] - ratio_stats[0]
    return (ratio_stats)

