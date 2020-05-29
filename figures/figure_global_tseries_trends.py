# ---------------------------------------------------------
# Save data for global time series and trends in GMT format
# ---------------------------------------------------------
import os
import numpy as np
import mod_gentools as gentools
def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['dir_gmt'] = os.getenv('HOME')+'/Scripts/GMT/Papers/Budget_20c/Global_tseries_trends/'
    settings['fn_stats'] = settings['dir_budget'] + 'results/stats_basin_global.npy'
    settings['years'] = np.arange(1900,2019)
    stats = np.load(settings['fn_stats'],allow_pickle=True).all()
    save_tseries(stats, settings)
    save_sliding_trends(stats, settings)
    save_linear_trends(stats, settings)
    return

def save_tseries(stats,settings):
    # Global time series of mean and uncertainty
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'obs_glb_tseries', settings['years'], stats['global']['obs']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'alt_glb_tseries', settings['years'], stats['global']['altimetry']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'steric_glb_tseries', settings['years'], stats['global']['steric']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'grd_glb_tseries', settings['years'], stats['global']['grd_total']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'glac_glb_tseries', settings['years'], stats['global']['grd_glac']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'GrIS_glb_tseries', settings['years'], stats['global']['grd_GrIS']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'AIS_glb_tseries', settings['years'], stats['global']['grd_AIS']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'tws_glb_tseries', settings['years'], stats['global']['grd_tws']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'dam_glb_tseries', settings['years'], stats['global']['grd_tws_dam']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'gwd_glb_tseries', settings['years'], stats['global']['grd_tws_gwd']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'nat_glb_tseries', settings['years'], stats['global']['grd_tws_natural']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'budget_tseries', settings['years'], stats['global']['budget']['tseries'])

    return

def save_sliding_trends(stats,settings):
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'obs_glb_sliding', settings['years'], stats['global']['obs']['sliding_trend'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'steric_glb_sliding', settings['years'], stats['global']['steric']['sliding_trend'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'grd_glb_sliding', settings['years'], stats['global']['grd_total']['sliding_trend'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'glac_glb_sliding', settings['years'], stats['global']['grd_glac']['sliding_trend'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'GrIS_glb_sliding', settings['years'], stats['global']['grd_GrIS']['sliding_trend'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'AIS_glb_sliding', settings['years'], stats['global']['grd_AIS']['sliding_trend'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'tws_glb_sliding', settings['years'], stats['global']['grd_tws']['sliding_trend'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'budget_sliding', settings['years'], stats['global']['budget']['sliding_trend'])
    return

def save_linear_trends(stats,settings):
    save_trend_periods('obs', stats['global']['obs'], [0,5,10], settings)
    save_trend_periods('budget', stats['global']['budget'], [1,6,11], settings)
    save_trend_periods('steric', stats['global']['steric'], [2,7,13], settings)
    save_trend_periods('grd', stats['global']['grd_total'], [3,8,14], settings)
    save_trend_periods('alt', stats['global']['altimetry'], [-1,-1,12], settings)

    save_trend_periods('grd_tot', stats['global']['grd_total'], [0,6,12], settings)
    save_trend_periods('glac', stats['global']['grd_glac'], [1,7,13], settings)
    save_trend_periods('GrIS', stats['global']['grd_GrIS'], [2,8,14], settings)
    save_trend_periods('AIS', stats['global']['grd_AIS'], [3,9,15], settings)
    save_trend_periods('tws', stats['global']['grd_tws'], [4,10,16], settings)
    return

def save_trend_periods(fn,term,pos_arr,settings):
    gentools.gmt_save_trend(settings['dir_gmt']+fn+'_glb_trend_1900.txt',pos_arr[0],term['trend_1900_2018'])
    gentools.gmt_save_trend(settings['dir_gmt']+fn+'_glb_trend_1957.txt',pos_arr[1],term['trend_1957_2018'])
    gentools.gmt_save_trend(settings['dir_gmt']+fn+'_glb_trend_1993.txt',pos_arr[2],term['trend_1993_2018'])
    return

def latexprint(trend):
    print(str(np.around(trend[0],2))[:5])
    print(str(np.around(trend[0]+trend[1],2))[:5])
    print(str(np.around(trend[0]+trend[2],2))[:5])
    return
