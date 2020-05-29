# ---------------------------------------------------------
# Time series and linear trends for each of the basins
# Save data for global time series and trends in GMT format
# ---------------------------------------------------------
import os
import numpy as np
import mod_gentools as gentools
def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['dir_gmt'] = os.getenv('HOME')+'/Scripts/GMT/Papers/Budget_20c/Basin_tseries/'
    settings['fn_stats'] = settings['dir_budget'] + 'results/stats_basin_global.npy'
    settings['years'] = np.arange(1900,2019)
    stats = np.load(settings['fn_stats'],allow_pickle=True).all()
    save_tseries(stats, settings)
    save_linear_trends(stats, settings)
    return

def save_tseries(stats,settings):
    # Global time series of mean and uncertainty
    for basin in range(6):
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'obs_'+   str(basin)+'_tseries', settings['years'], stats['basin'][basin]['obs']['tseries'])
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'alt_'+   str(basin)+'_tseries', settings['years'], stats['basin'][basin]['altimetry']['tseries'])
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'steric_'+str(basin)+'_tseries', settings['years'], stats['basin'][basin]['steric']['tseries'])
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'grd_'+   str(basin)+'_tseries', settings['years'], stats['basin'][basin]['grd_total']['tseries'])
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'gia_'+   str(basin)+'_tseries', settings['years'], stats['basin'][basin]['gia']['tseries'])
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'budget_'+str(basin)+'_tseries', settings['years'], stats['basin'][basin]['budget']['tseries'])
    return

def save_linear_trends(stats,settings):
    for basin in range(6):
        save_trend_periods('obs_'+   str(basin), stats['basin'][basin]['obs'], [0,4,9], settings)
        save_trend_periods('budget_'+str(basin), stats['basin'][basin]['budget'], [-1,5,11], settings)
        save_trend_periods('steric_'+str(basin), stats['basin'][basin]['steric'], [-1,7,13], settings)
        save_trend_periods('grd_'+   str(basin), stats['basin'][basin]['grd_total'], [2,6,12], settings)
        save_trend_periods('alt_'+   str(basin), stats['basin'][basin]['altimetry'], [-1,-1,10], settings)
        save_trend_periods('gia_'+   str(basin), stats['basin'][basin]['gia'], [1,-1,-1], settings)
    return

def save_trend_periods(fn,term,pos_arr,settings):
    gentools.gmt_save_trend(settings['dir_gmt']+fn+'_trend_1900.txt',pos_arr[0],term['trend_1900_2018'])
    gentools.gmt_save_trend(settings['dir_gmt']+fn+'_trend_1957.txt',pos_arr[1],term['trend_1957_2018'])
    gentools.gmt_save_trend(settings['dir_gmt']+fn+'_trend_1993.txt',pos_arr[2],term['trend_1993_2018'])
    return