# --------------------------------------------------------
# Compute virtual station for each region from region list
# 1. Define ensemble members
# 2. Determine region order
# 3. Merge for each ensemble member
# 4. Compute global estimate from basin-mean data
# --------------------------------------------------------
import numpy as np
import mod_gentools as gentools
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def main():
    settings = {}
    settings['test_run_ICE6G_D'] = True
    settings['dir_data']    = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_basin'] = settings['dir_data']+'Basins/ocean_regions_thompson.grd'
    settings['years'] = np.arange(1900,2019)
    settings['min_years_merge'] = 10
    settings['num_ens'] = 100
    if settings['test_run_ICE6G_D']:
        settings['fn_region_ensembles'] = settings['dir_budget']+'region_data/region_ensembles_ice6g.npy'
        settings['fn_obs_ensembles']    = settings['dir_budget']+'results/obs_basin_global_ens_ice6g.npy'
    else:
        settings['fn_region_ensembles'] = settings['dir_budget']+'region_data/region_ensembles.npy'
        settings['fn_obs_ensembles']    = settings['dir_budget']+'results/obs_basin_global_ens.npy'
    settings['fn_basin_terms'] = settings['dir_data'] + 'Budget_20c/results/basin_terms.npy'
    settings['virstat_visual'] = False # Show plots of merging process
    region_ensembles = np.load(settings['fn_region_ensembles'],allow_pickle=True)
    station_order,merge_info = determine_station_order(region_ensembles, settings)
    obs_ensembles = {}
    obs_ensembles['basin'] = compute_basin_ensemble(region_ensembles, station_order, merge_info, settings)
    obs_ensembles['global'] = compute_global_ensembles(obs_ensembles['basin'], settings)
    np.save(settings['fn_obs_ensembles'],obs_ensembles)
    return

def compute_global_ensembles(basin_ensembles,settings):
    # Compute weights of each basin
    # Take into account that not all basins are complete
    basin_mask = read_basin_mask(settings)
    area = gentools.grid_area(basin_mask['lat'],basin_mask['lon'])
    slm = np.isfinite(basin_mask['num'])
    total_area = np.sum(area*slm)
    basin_weight = np.zeros(len(basin_ensembles))
    for basin in range(len(basin_ensembles)): basin_weight[basin] = np.sum(area*(basin_mask['num']==basin)) / total_area
    reg_acc = np.zeros([len(basin_ensembles),len(settings['years'])],dtype=bool)
    for basin in range(len(basin_ensembles)):  reg_acc[basin,:] = np.isfinite(basin_ensembles[basin]['tg_no_corrections'][0, :])
    weight_idx = reg_acc * basin_weight[:,np.newaxis]
    weight_idx = weight_idx / weight_idx.sum(axis=0)
    global_ensembles = {}
    for method in basin_ensembles[0].keys():
        global_ensembles[method] = np.zeros([settings['num_ens'],len(settings['years'])])
        for basin in range(len(basin_ensembles)):
            basin_lcl = basin_ensembles[basin][method].copy()
            basin_lcl[np.isnan(basin_lcl)] = 0
            global_ensembles[method] = global_ensembles[method] + weight_idx[basin,:] * basin_lcl
    return(global_ensembles)

def compute_basin_ensemble(region_ensembles, station_order, merge_info, settings):
    print('Compute basin estimates:')
    corr_runs = ['tg_no_corrections', 'tg_gia', 'tg_gia_grd', 'tg_resvlm', 'tg_full']
    basin_ensembles = np.zeros(len(region_ensembles),dtype=object)
    # if settings['virstat_visual']: basin_terms = np.load(settings['fn_basin_terms'],allow_pickle=True).all()
    for basin in range(len(region_ensembles)):
        print('   Basin '+str(basin))
        basin_ensembles[basin] = {}
        for run in corr_runs:
            region_lcl = region_ensembles[basin][run].copy()
            for region_vs in range(len(station_order[basin])):
                combine_height = np.zeros([2,settings['num_ens'],len(settings['years'])])
                combine_height[0,...] = region_lcl[:,station_order[basin][region_vs][0],:]
                combine_height[1,...] = region_lcl[:,station_order[basin][region_vs][1],:]
                num_obs = np.sum(np.isfinite(combine_height[:,0,:]),axis=0).astype(float)
                num_obs[num_obs == 0] = np.nan
                combine_height[0,:,:] = combine_height[0,:,:] - np.mean(combine_height[0,:, num_obs == 2],axis=0)[:,np.newaxis]
                combine_height[1,:,:] = combine_height[1,:,:] - np.mean(combine_height[1,:, num_obs == 2],axis=0)[:,np.newaxis]
                new_height = np.nansum(combine_height, axis=0) / num_obs[np.newaxis,:]
                region_lcl = np.append(region_lcl,new_height[:,np.newaxis,:],axis=1)
                region_lcl = np.delete(region_lcl,station_order[basin][region_vs],axis=1)
                if settings['virstat_visual'] & (run == 'tg_full'):
                    print(merge_info[basin]['station_names'][region_vs])
                    fig = plt.figure(figsize=(6, 9))
                    ax1 = fig.add_subplot(3, 1, 1, projection=ccrs.PlateCarree())
                    ax1.coastlines()
                    if merge_info[basin]['untouched_stats'][region_vs].shape[0]>0: ax1.plot(merge_info[basin]['untouched_stats'][region_vs][:,1],merge_info[basin]['untouched_stats'][region_vs][:,0],'x',transform=ccrs.PlateCarree())
                    ax1.plot(merge_info[basin]['merged_stats'][region_vs][:,1],merge_info[basin]['merged_stats'][region_vs][:,0],'o',transform=ccrs.PlateCarree())
                    ax1.plot(merge_info[basin]['virstat'][region_vs][1],merge_info[basin]['virstat'][region_vs][0],'o',color='red',transform=ccrs.PlateCarree())
                    ax1.set_extent([0, 359, -88, 88], ccrs.PlateCarree())
                    ax2 = fig.add_subplot(3, 1, 2)
                    ax2.plot(settings['years'],combine_height[0,:,:].mean(axis=0))
                    ax2.plot(settings['years'],combine_height[1,:,:].mean(axis=0))
                    ax2.plot(settings['years'],new_height.mean(axis=0))
                    # if (np.isfinite(new_height.mean(axis=0)) & np.isfinite(basin_terms['basin'][basin]['total']['tseries'][:,1])).sum() > 5:
                    #     acc_idx = (np.isfinite(new_height.mean(axis=0)) & np.isfinite(basin_terms['basin'][basin]['total']['tseries'][:,1]))
                    #     ax2.plot(settings['years'],basin_terms['basin'][basin]['total']['tseries'][:,1]-np.nanmean(basin_terms['basin'][basin]['total']['tseries'][acc_idx,1]) + new_height.mean(axis=0)[acc_idx].mean())
                    #
                    # else:
                    #     ax2.plot(settings['years'],basin_terms['basin'][basin]['total']['tseries'][:,1]-np.nanmean(basin_terms['basin'][basin]['total']['tseries'][:,1]))
                    plt.grid()
                    ax3 = fig.add_subplot(3, 1, 3)
                    compute_trend_smooth(combine_height[0,:,:].mean(axis=0), settings)
                    ax3.plot(settings['years'],compute_trend_smooth(combine_height[0,:,:].mean(axis=0), settings))
                    ax3.plot(settings['years'],compute_trend_smooth(combine_height[1,:,:].mean(axis=0), settings))
                    ax3.plot(settings['years'],compute_trend_smooth(new_height.mean(axis=0), settings))
                    # ax3.plot(settings['years'],compute_trend_smooth(basin_terms['basin'][basin]['total']['tseries'][:,1]-np.nanmean(basin_terms['basin'][basin]['total']['tseries'][:,1]), settings))

                    plt.draw()
                    plt.pause(0.001)
                    input("Press Enter to continue...")
                    plt.close()
            basin_ensembles[basin][run] = region_lcl.squeeze()
    return(basin_ensembles)

def determine_station_order(region_ensembles,settings):
    # ----------------------------------------------------
    # Determine the order in which stations must be merged
    # ----------------------------------------------------
    print('Determining station order:')
    station_order = np.zeros(len(region_ensembles),dtype=object)
    merge_info    = np.zeros(len(region_ensembles),dtype=object)
    for basin in range(len(region_ensembles)):
        print('   Basin '+str(basin))
        # Rebuild station names
        station_names = region_ensembles[basin]['station_names'][:]
        lat = region_ensembles[basin]['coords'][:,0]
        lon = region_ensembles[basin]['coords'][:,1]
        time_acc = np.isfinite(region_ensembles[basin]['height'])
        merge_order = []
        merge_info[basin] = {}
        merge_info[basin]['untouched_stats'] = []
        merge_info[basin]['merged_stats']    = []
        merge_info[basin]['virstat']       = []
        merge_info[basin]['station_names'] = []

        run_flag    = True
        while run_flag:
            latmid,lonmid,merged_stats,virstat_time_acc = compute_midpoint_min_overlap(lat, lon, time_acc, settings)
            merge_info[basin]['merged_stats'].append(np.array([lat[merged_stats],lon[merged_stats]]).T)
            # Station names from merged stats
            virstat_names = station_names[merged_stats][0] + station_names[merged_stats][1]
            merge_info[basin]['station_names'].append(virstat_names)
            lat = np.delete(lat,merged_stats)
            lon = np.delete(lon,merged_stats)
            station_names = np.delete(station_names,merged_stats)
            station_new = np.zeros(len(station_names)+1,dtype=object)
            station_new[:-1] = station_names
            station_new[-1] = virstat_names
            station_names = station_new
            merge_info[basin]['untouched_stats'].append(np.array([lat,lon]).T)
            merge_info[basin]['virstat'].append(np.array([latmid,lonmid]))
            merge_order.append(merged_stats)
            # Delete merged stations
            time_acc = np.delete(time_acc,merged_stats,axis=0)
            # Append virtual station
            lat = np.append(lat,latmid)
            lon = np.append(lon,lonmid)
            time_acc = np.append(time_acc,virstat_time_acc[np.newaxis,:],axis=0)
            if len(lat)==1: run_flag=False
        station_order[basin] = np.array(merge_order)
    print('Done')
    return(station_order,merge_info)

## HELPER FUNCTIONS ##
def compute_midpoint_min_overlap(lat,lon,time_acc,settings):
    lonmat,latmat = np.meshgrid(lon,lat)
    distmat = np.arcsin(np.sqrt(np.sin(np.deg2rad(0.5 * (latmat - latmat.T))) ** 2 + np.cos(np.deg2rad(latmat.T)) * np.cos(np.deg2rad(latmat)) * np.sin(np.deg2rad(0.5 * (lonmat - lonmat.T))) ** 2))
    distmat = distmat + np.tri(N=len(distmat),k=0)*1e10
    ovl_notfound = True
    while ovl_notfound:
        merged_stats = np.array(np.unravel_index(np.argmin(distmat), distmat.shape))
        n_ovl = np.sum(time_acc[merged_stats[0], :] & time_acc[merged_stats[1], :])
        if n_ovl >= settings['min_years_merge']: ovl_notfound=False
        else: distmat[merged_stats[0], merged_stats[1]] = 1e10
     # Midpoint computation
    Bx = np.cos(np.deg2rad(lat[merged_stats[1]])) * np.cos(np.deg2rad(lon[merged_stats[1]] - lon[merged_stats[0]]))
    By = np.cos(np.deg2rad(lat[merged_stats[1]])) * np.sin(np.deg2rad(lon[merged_stats[1]] - lon[merged_stats[0]]))

    latmid = np.rad2deg(np.arctan2(np.sin(np.deg2rad(lat[merged_stats[0]])) + np.sin(np.deg2rad(lat[merged_stats[1]])), np.sqrt((np.cos(np.deg2rad(lat[merged_stats[0]])) + Bx) ** 2 + By ** 2)))
    lonmid = np.rad2deg(np.deg2rad(lon[merged_stats[0]]) + np.arctan2(By, np.cos(np.deg2rad(lat[merged_stats[0]])) + Bx))
    # Accepted points of new virtual station
    virstat_time_acc = time_acc[merged_stats[0], :] | time_acc[merged_stats[1], :]
    return(latmid,lonmid,merged_stats,virstat_time_acc)

def read_basin_mask(settings):
    basin_mask = {}
    file_handle = Dataset(settings['fn_basin'], 'r')
    file_handle.set_auto_mask(False)
    basin_mask['lat'] = file_handle.variables['y'][:]
    basin_mask['lon'] = file_handle.variables['x'][:]
    basin_mask['num'] = file_handle.variables['z'][:]
    basin_mask['regions'] = np.arange(0, 6)
    file_handle.close()
    basin_mask['num'][basin_mask['num'] == 0] = np.nan
    basin_mask['num'] = basin_mask['num']-1
    return (basin_mask)

def compare_methods(obs_ensembles,settings):
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.figure(figsize=(9,8))
    for basin in range(6):
        ax = plt.subplot(321+basin)
        for idx, prod in enumerate(obs_ensembles['basin'][basin]):
            plt.plot(settings['years'],obs_ensembles['basin'][basin][prod].mean(axis=0),label=prod)
    plt.grid()
    plt.legend(fontsize=8)
    plt.axis([2005, 2017, -10, 40])
    plt.ylabel('Height (mm)',fontsize=8)
    ax2 = plt.subplot(312)
    return

def compute_trend_smooth(tseries,settings):
    filt_width=30
    fhalf = np.int(filt_width/2)
    trend_smooth = np.zeros(len(settings['years']))*np.nan
    for idx, yr in enumerate(settings['years']):
        idx_start = np.max([0,idx-fhalf])
        idx_stop  = np.min([len(settings['years']),idx+fhalf])
        yr_acc=settings['years'][idx_start:idx_stop]
        if np.isfinite(tseries[idx_start:idx_stop]).sum()>24:
            tg_acc = np.isfinite(tseries[idx_start:idx_stop])
            amat = np.ones([len(yr_acc), 2])
            amat[:, 1] = yr_acc
            amat = amat[tg_acc,:]
            trend_smooth[idx] = np.linalg.lstsq(amat, tseries[idx_start:idx_stop][tg_acc], rcond=None)[0][1]
    return(trend_smooth)

if __name__ == '__main__':
    main()
