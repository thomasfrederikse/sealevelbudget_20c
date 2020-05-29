# ---------------------------------------------------------------
# Compare the budget and sea level between LC's model and ICE6G_D
# Plot figures for response letter
# ---------------------------------------------------------------
import numpy as np
import os
import matplotlib.pyplot as plt
import mod_gentools as gentools
def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_stats_ice6g'] = settings['dir_budget'] + 'results/stats_basin_global_ice6g.npy'
    settings['fn_stats_caron'] = settings['dir_budget'] + 'results/stats_basin_global.npy'
    settings['dir_gmt'] = os.getenv('HOME')+'/Scripts/GMT/Papers/Budget_20c/comp_ICE6G/'

    settings['years'] = np.arange(1900, 2019)

    stats_caron = np.load(settings['fn_stats_caron'],allow_pickle=True).all()
    stats_ice6g = np.load(settings['fn_stats_ice6g'],allow_pickle=True).all()

    # Save basin data
    for basin in range(6):
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'obs_caron_'+   str(basin)+'_tseries', settings['years'], stats_caron['basin'][basin]['obs']['tseries'])
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'bdg_caron_'+   str(basin)+'_tseries', settings['years'], stats_caron['basin'][basin]['budget']['tseries'])
        gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'gia_caron_'+   str(basin)+'_tseries', settings['years'], stats_caron['basin'][basin]['gia']['tseries'])

        gentools.gmt_save_tseries(settings['dir_gmt'], 'obs_ice6g_'+   str(basin)+'_tseries', settings['years'], stats_ice6g['basin'][basin]['obs']['tseries'][:,1])
        gentools.gmt_save_tseries(settings['dir_gmt'], 'bdg_ice6g_'+   str(basin)+'_tseries', settings['years'], stats_ice6g['basin'][basin]['budget']['tseries'][:,1])
        gentools.gmt_save_tseries(settings['dir_gmt'], 'gia_ice6g_'+   str(basin)+'_tseries', settings['years'], stats_ice6g['basin'][basin]['gia']['tseries'][:,1])

    # Global
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'obs_caron_glb_tseries', settings['years'], stats_caron['global']['obs']['tseries'])
    gentools.gmt_save_tseries_ci(settings['dir_gmt'], 'bdg_caron_glb_tseries', settings['years'], stats_caron['global']['budget']['tseries'])
    gentools.gmt_save_tseries(settings['dir_gmt'], 'obs_ice6g_glb_tseries', settings['years'], stats_ice6g['global']['obs']['tseries'][:, 1])
    gentools.gmt_save_tseries(settings['dir_gmt'], 'bdg_ice6g_glb_tseries', settings['years'], stats_ice6g['global']['budget']['tseries'][:, 1])

    # Plot global
    plt.figure(figsize=(5,4))
    plt.fill_between(settings['years'],stats_caron['global']['obs']['tseries'][:,0],stats_caron['global']['obs']['tseries'][:,2],color='C0',alpha=0.4)
    plt.fill_between(settings['years'],stats_caron['global']['budget']['tseries'][:,0],stats_caron['global']['budget']['tseries'][:,2],color='C1',alpha=0.4)
    plt.plot(settings['years'],stats_caron['global']['obs']['tseries'][:,1],'C0',linewidth=2,label='Observed sea level (This study)')
    plt.plot(settings['years'],stats_ice6g['global']['obs']['tseries'][:,1],'black',linestyle='--',linewidth=2,label='Observed sea level (ICE6G_D)')
    plt.plot(settings['years'],stats_caron['global']['budget']['tseries'][:,1],'C1',linewidth=2,label='Sum of processes (This study)')
    plt.plot(settings['years'],stats_ice6g['global']['budget']['tseries'][:,1],'C3',linestyle='--',linewidth=2,label='Sum of processes (ICE6G_D)')

    plt.xlim([1900, 2018])
    plt.ylim([-220, 50])
    plt.ylabel('Height (mm)',fontsize=9)
    plt.legend(fontsize=9)
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.title('Global',fontsize=9)
    plt.grid()
    plt.tight_layout()
    plt.savefig('Compare_ICE6G_D_global.png')

    # Plot basin
    region_names = ['Subpolar North Atlantic', 'Indian Ocean-South Pacific', 'Subtropical North Atlantic', 'East Pacific', 'South Atlantic', 'Northwest Pacific']
    plt.figure(figsize=(8,8))
    for basin in range(6):
        plt.subplot(3, 2, basin+1)
        plt.fill_between(settings['years'], stats_caron['basin'][basin]['obs']['tseries'][:, 0], stats_caron['basin'][basin]['obs']['tseries'][:, 2], color='C0', alpha=0.4)
        plt.fill_between(settings['years'], stats_caron['basin'][basin]['budget']['tseries'][:, 0], stats_caron['basin'][basin]['budget']['tseries'][:, 2], color='C1', alpha=0.4)
        plt.plot(settings['years'], stats_caron['basin'][basin]['obs']['tseries'][:, 1], 'C0', linewidth=2, label='Observed sea level (This study)')
        plt.plot(settings['years'], stats_ice6g['basin'][basin]['obs']['tseries'][:, 1], 'black', linestyle='--', linewidth=2, label='Observed sea level (ICE6G_D)')
        plt.plot(settings['years'], stats_caron['basin'][basin]['budget']['tseries'][:, 1], 'C1', linewidth=2, label='Sum of processes (This study)')
        plt.plot(settings['years'], stats_ice6g['basin'][basin]['budget']['tseries'][:, 1], 'C3', linestyle='--', linewidth=2, label='Sum of processes (ICE6G_D)')
        plt.xlim([1900, 2018])
        plt.ylim([-400, 100])
        plt.ylabel('Height (mm)', fontsize=9)
        if basin==0: plt.legend(fontsize=9)
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
        plt.title(region_names[basin], fontsize=9)
        plt.grid()
        plt.tight_layout()
    plt.savefig('Compare_ICE6G_D_basin.png')
    return

