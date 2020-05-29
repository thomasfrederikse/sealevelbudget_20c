# ------------------------------------------------
# Read TG data as a preliminary data set used to
# manually check all regions for outliers
# ------------------------------------------------
# - Merge nearby records
# - Merge nearby stations into regional estimates
# - Remove meteorological forcing and nodal cycle
# - Compute basin estimates for region selection
# ------------------------------------------------
import numpy as np
import os
from netCDF4 import Dataset
import mod_gentools as gentools
import scipy.stats as scistats

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'

    settings['fn_nodal']   = settings['dir_data']+'Nodal/Nodal.nc'
    settings['fn_ERA']     = settings['dir_budget']+'tg/ERA.nc'
    settings['fn_ERA_slm'] = settings['dir_budget']+'tg/ERA_slm.nc'
    settings['fn_basin'] = settings['dir_data']+'Basins/ocean_regions_thompson.grd'
    settings['fn_gia_rsl'] = settings['dir_data']+'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['fn_station_data'] = settings['dir_budget']+'tg/station_data.npy'
    settings['fn_regions_for_selection'] = settings['dir_budget']+'tg/regions_for_selection.npy'

    settings['years'] = np.arange(1900,2019)
    settings['num_ens'] = 100
    settings['min_ovl']   = 5
    settings['merge_dist'] = 20000

    station_data = np.load(settings['fn_station_data'],allow_pickle=True).all()
    regions_for_selection = merge_nearby_stations(station_data, settings)
    regions_for_selection = attribute_station_to_basin(regions_for_selection, settings)
    regions_for_selection = compute_basin_estimate(regions_for_selection, settings)
    save_data(regions_for_selection, settings)
    return

def compute_basin_estimate(regions_for_selection,settings):
    print('Computing basin estimates:')
    # --------------------------------------------
    # Read GIA and GRD and compute basin estimates
    # --------------------------------------------
    basin_mask = read_basin_mask(settings)

    # Sample locations
    regions_for_selection['grd_sample_points'] = np.zeros([len(regions_for_selection['id']),2],dtype=int)
    lat = np.arange(-89.75,90.25,0.5)
    lon = np.arange(0.25,360.25,0.5)
    regions_for_selection['grd_sample_points'][:,0] = np.argmin(np.abs(regions_for_selection['coords'][:, 0][np.newaxis, :] - lat[:, np.newaxis]), axis=0)
    regions_for_selection['grd_sample_points'][:,1] = np.argmin(np.abs(regions_for_selection['coords'][:, 1][np.newaxis, :] - lon[:, np.newaxis]), axis=0)

    # Sample GIA at region locations
    GIA = read_GIA_rsl(settings)

    regions_for_selection['rsl_gia_mean'] = np.zeros(len(regions_for_selection['id']))
    regions_for_selection['rsl_gia_dev'] = np.zeros(len(regions_for_selection['id'])) # Deviation of region from basin
    for region in range(len(regions_for_selection['id'])):
        regions_for_selection['rsl_gia_mean'][region] = (GIA['probability'] * GIA['rsl'][:, regions_for_selection['grd_sample_points'][region, 0], regions_for_selection['grd_sample_points'][region, 1]]).sum()
    GIA_basin = np.zeros(len(basin_mask['basins']))  # Basin-mean GIA
    area = gentools.grid_area(GIA['lat'],GIA['lon'])
    for basin in range(len(basin_mask['basins'])):
        GIA_basin[basin] = (GIA['probability']*(((area * (basin_mask['num'] == basin))[np.newaxis,:,:] * GIA['rsl']).sum(axis=(1,2)) / (area * (basin_mask['num'] == basin)).sum())).sum()
    for region in range(len(regions_for_selection['id'])):
        regions_for_selection['rsl_gia_dev'][region] = regions_for_selection['rsl_gia_mean'][region] - GIA_basin[regions_for_selection['basin_num'][region]]

    # Sample GRD at region locations
    regions_for_selection['rsl_grd_mean'] = np.zeros([len(regions_for_selection['id']),len(settings['years'])])
    regions_for_selection['rsl_grd_dev']  = np.zeros([len(regions_for_selection['id']),len(settings['years'])]) # Deviation of region from basin
    grd_region_ens = np.zeros([settings['num_ens'],len(regions_for_selection['id']),len(settings['years'])])
    grd_basin_ens = np.zeros([settings['num_ens'],len(basin_mask['basins']),len(settings['years'])])
    for ens in range(settings['num_ens']):
        print('   Ensemble '+str(ens))
        GRD_rsl_ens = read_GRD_rsl_ens(ens, settings)
        for region in range(len(regions_for_selection['id'])):
            grd_region_ens[ens,region,:] = GRD_rsl_ens[:, regions_for_selection['grd_sample_points'][region, 0], regions_for_selection['grd_sample_points'][region, 1]]
        for basin in range(len(basin_mask['basins'])):
            grd_basin_ens[ens,basin,:] = ((area * (basin_mask['num'] == basin))[np.newaxis, :, :] * GRD_rsl_ens).sum(axis=(1, 2)) / (area * (basin_mask['num'] == basin)).sum()
    grd_basin = (GIA['probability'][:,np.newaxis,np.newaxis]*grd_basin_ens).sum(axis=0)
    grd_region = (GIA['probability'][:,np.newaxis,np.newaxis]*grd_region_ens).sum(axis=0)
    for region in range(len(regions_for_selection['id'])):
        regions_for_selection['rsl_grd_mean'][region,:] = grd_region[region,:]
        regions_for_selection['rsl_grd_dev'][region] = grd_region[region,:] - grd_basin[regions_for_selection['basin_num'][region],:]
    return(regions_for_selection)

def merge_nearby_stations(station_data,settings):
    print('Merge nearby stations into regional estimates:')
    # -------------------------------------------
    # Merge nearby stations into region estimates
    # Store results in list
    # -------------------------------------------
    untouched_stations = np.ones(len(station_data['id']),dtype=bool)
    regions_for_selection = {}
    regions_for_selection['id']     = []
    regions_for_selection['name']   = []
    regions_for_selection['coords'] = []
    regions_for_selection['height_corr'] = []
    while untouched_stations.sum()>0:
        workstat   = np.where(untouched_stations)[0][0] # Index of first untouched station:
        dist_array = 2*6371000*np.arcsin(np.sqrt(np.sin(np.deg2rad(0.5*(station_data['coords'][:,0]-station_data['coords'][workstat,0])))**2+np.cos(np.deg2rad(station_data['coords'][workstat,0]))*np.cos(np.deg2rad(station_data['coords'][:,0]))*np.sin(np.deg2rad(0.5*(station_data['coords'][:,1]-station_data['coords'][workstat,1])))**2))
        acc_merge  = (dist_array<settings['merge_dist']) & (untouched_stations)
        if acc_merge.sum()==1: # No stations to merge
            untouched_stations[workstat]=False
            regions_for_selection['id'].append([station_data['id'][workstat]])
            regions_for_selection['name'].append(station_data['name'][workstat])
            regions_for_selection['coords'].append(station_data['coords'][workstat])
            regions_for_selection['height_corr'].append(station_data['height_corr'][workstat])
        else:
            merge_idx   = np.where(dist_array<settings['merge_dist'])[0]
            merge_array = np.zeros([len(settings['years']), len(merge_idx)])
            for idx,num in enumerate(merge_idx): merge_array[:, idx] = station_data['height_corr'][num]
            ts_merged, is_station_merged = merge_lcl_stats(merge_array.copy(), settings)
            if is_station_merged.sum()>0: # Stations have been merged: new entry in regions_for_selection
                merged_stats = merge_idx[is_station_merged]
                id_list = []
                name_list = []
                for idx,statnum in enumerate(merged_stats):
                    id_list.append(station_data['id'][statnum])
                    name_list.append(station_data['name'][statnum])
                    untouched_stations[statnum] = False
                regions_for_selection['id'].append(id_list)
                regions_for_selection['name'].append(name_list)
                regions_for_selection['coords'].append(station_data['coords'][merged_stats[0]])
                regions_for_selection['height_corr'].append(ts_merged)
            else: # No stations have been merged
                untouched_stations[workstat] = False
                regions_for_selection['id'].append([station_data['id'][workstat]])
                regions_for_selection['name'].append(station_data['name'][workstat])
                regions_for_selection['coords'].append(station_data['coords'][workstat])
                regions_for_selection['height_corr'].append(station_data['height_corr'][workstat])

    regions_for_selection['id']     = np.array(regions_for_selection['id'])
    regions_for_selection['name']   = np.array(regions_for_selection['name'])
    regions_for_selection['coords'] = np.array(regions_for_selection['coords'])
    regions_for_selection['height_corr'] = np.array(regions_for_selection['height_corr'])
    return(regions_for_selection)


def attribute_station_to_basin(regions_for_selection,settings):
    print('Attribute station to basin:')
    regions_for_selection['basin_num'] = np.zeros(len(regions_for_selection['id']),dtype=int)
    # Determine basin to which each station belongs
    basin_mask = read_basin_mask(settings)
    for region in range(len(regions_for_selection['id'])):
        if (regions_for_selection['coords'][region,0]<55) & (regions_for_selection['coords'][region,0]>35.79)& (regions_for_selection['coords'][region,1]>280)&(regions_for_selection['coords'][region,1]<310):
            regions_for_selection['basin_num'][region] = 0
        else:
            distance = gentools.point_grid_distance(regions_for_selection['coords'][region,0],regions_for_selection['coords'][region,1],basin_mask['lat'],basin_mask['lon'])
            distance[np.isnan(basin_mask['num'])] = 1e9
            if basin_mask['num'][distance < 250000].size==0: regions_for_selection['basin_num'][region] = basin_mask['num'].flatten()[distance.argmin()].astype(int)
            else: regions_for_selection['basin_num'][region] = scistats.mode(basin_mask['num'][distance < 250000].astype(int))[0][0]
    return(regions_for_selection)

def save_data(regions_for_selection,settings):
    print('Saving data:')
    np.save(settings['fn_regions_for_selection'],regions_for_selection)
    return


## HELPER FUNCTIONS ##

def merge_lcl_stats(merge_array,settings):
    # Merge nearby stations in merge_array into single record.
    # Only merge if enough overlap
    is_station_merged = np.zeros(merge_array.shape[1],dtype=bool)
    is_station_touched = np.zeros(merge_array.shape[1],dtype=bool)
    is_station_touched[0] = True
    while is_station_touched.sum()<(merge_array.shape[1]):
        n_ovl = np.zeros(merge_array.shape[1]-1,dtype=int) # Determine number of overlaps
        for stat in range(len(n_ovl)):
            n_ovl[stat] = (np.isfinite(merge_array[:, 0] * merge_array[:,stat+1])).sum()
        if n_ovl.max()>=settings['min_ovl']:
            merge_array_lcl = np.zeros([len(settings['years']),2])
            merge_array_lcl[:,0] = merge_array[:,0]
            merge_array_lcl[:,1] = merge_array[:,n_ovl.argmax()+1]
            merge_array[:,n_ovl.argmax()+1] = np.nan
            is_station_merged[n_ovl.argmax()+1]=True
            is_station_merged[0]=True
            merge_array[:,0] = gentools.merge_common_mean(merge_array_lcl)
        else: # No single station has enough overlap:
            is_station_touched[:] = True
        is_station_touched[np.where(n_ovl==n_ovl.max())[0]+1]=True
    merged_tseries = merge_array[:,0]
    return(merged_tseries,is_station_merged)

def read_basin_mask(settings):
    basin_mask = {}
    file_handle = Dataset(settings['fn_basin'], 'r')
    file_handle.set_auto_mask(False)
    basin_mask['lat'] = file_handle.variables['y'][:]
    basin_mask['lon'] = file_handle.variables['x'][:]
    basin_mask['num'] = file_handle.variables['z'][:]
    basin_mask['basins'] = np.arange(0, 6)
    file_handle.close()
    basin_mask['num'][basin_mask['num'] == 0] = np.nan
    basin_mask['num'] = basin_mask['num'] - 1
    return (basin_mask)

def read_GIA_rsl(settings):
    GIA = {}
    file_handle = Dataset(settings['fn_gia_rsl'], 'r')
    file_handle.set_auto_mask(False)
    GIA['lat'] =  file_handle.variables['y'][:]
    GIA['lon'] =  file_handle.variables['x'][:]
    GIA['probability'] = file_handle.variables['probability'][:settings['num_ens']]
    GIA['probability'] = GIA['probability']/GIA['probability'].sum()
    GIA['rsl'] = file_handle.variables['rsl'][:settings['num_ens'],:,:]
    file_handle.close()
    return(GIA)

def read_GRD_rsl_ens(ens,settings):
    file_handle = Dataset(settings['dir_data']+'Budget_20c/grd/grd_'+str(ens)+'.nc', 'r')
    file_handle.set_auto_mask(False)
    GRD_rsl_ens = file_handle.variables['rsl'][:]
    file_handle.close()
    return(GRD_rsl_ens)