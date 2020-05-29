# Save data for region map in GMT
import numpy as np
import os

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_basin'] = settings['dir_data'] + 'Basins/ocean_regions_thompson.grd'

    settings['fn_station_data'] = settings['dir_budget']+'tg/station_data.npy'
    settings['fn_region_list'] = settings['dir_budget']+'region_data/region_list_beta_march_9.npy'
    settings['years']     = np.arange(1900,2019)
    settings['dir_gmt'] = os.getenv('HOME') + '/Scripts/GMT/Papers/Budget_20c/Station_map/'

    station_data = np.load(settings['fn_station_data'],allow_pickle=True).all()
    region_list  = np.load(settings['fn_region_list'],allow_pickle=True)


    # Determine total number of tide gauge stations
    tg_stations = []
    for basin in range(len(region_list)):
        for region in range(len(region_list[basin]['list'])):
            for stat in region_list[basin]['list'][region]['id']:
                tg_stations.append(stat)
    tg_stations = np.unique(np.array(tg_stations))
    stations_per_basin = np.zeros(len(region_list))
    for basin in range(len(region_list)):
        novlm = []
        vlm_gps = []
        vlm_alttg = []
        vlm_all = []
        num_tg = np.zeros(len(settings['years']))
        for region in range(len(region_list[basin]['list'])):
            stations_per_basin[basin]+=len(region_list[basin]['list'][region]['id'])
            stat_data = np.ones(4)
            acc_tg = np.zeros(len(settings['years']),dtype=bool)
            for stat in region_list[basin]['list'][region]['id']:
                stat_idx = station_data['id']==stat
                acc_tg[np.isfinite(station_data['height'][stat_idx].flatten())] = True
            stat_data[:2]  = np.fliplr(station_data['coords'][stat_idx])
            num_tg += acc_tg*1.0
            stat_data[3] = acc_tg.sum()
            if len(region_list[basin]['list'][region]['vlm_id'])==0:
                novlm.append(stat_data)
            elif 'ALTTG' in region_list[basin]['list'][region]['vlm_id']:
                if len(region_list[basin]['list'][region]['vlm_id'])>1:
                    vlm_all.append(stat_data)
                else:
                    vlm_alttg.append(stat_data)
            else:
                vlm_gps.append(stat_data)
        novlm = np.array(novlm)
        novlm=novlm[np.argsort(novlm[:,3]),:]

        vlm_gps = np.array(vlm_gps)
        vlm_gps=vlm_gps[np.argsort(vlm_gps[:,3]),:]

        vlm_alttg = np.array(vlm_alttg)
        if len(vlm_alttg>0): vlm_alttg=vlm_alttg[np.argsort(vlm_alttg[:,3]),:]

        vlm_all = np.array(vlm_all)
        if len(vlm_all>0): vlm_all=vlm_all[np.argsort(vlm_all[:,3]),:]

        novlm[:,3] = 0.75*(novlm[:,3])**0.6
        vlm_gps[:,3] = 0.75*(vlm_gps[:,3])**0.6
        np.savetxt(settings['dir_gmt']+'novlm_'+str(basin)+'.txt',novlm)
        np.savetxt(settings['dir_gmt']+'vlm_gps_'+str(basin)+'.txt',vlm_gps)
        if len(vlm_alttg>0):
            vlm_alttg[:,3] =  0.75*(vlm_alttg[:,3])**0.6
            np.savetxt(settings['dir_gmt'] + 'vlm_alttg_' + str(basin) + '.txt', vlm_alttg)
        if len(vlm_all>0):
            vlm_all[:,3] =  0.75*(vlm_all[:,3])**0.6
            np.savetxt(settings['dir_gmt']+'vlm_all_'+str(basin)+'.txt',vlm_all)
        np.savetxt(settings['dir_gmt']+'num_tg'+str(basin)+'.txt',np.vstack([settings['years'],num_tg]).astype(int).T,fmt='%3u; %3u')
    return