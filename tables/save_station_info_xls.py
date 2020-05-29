# -------------------------------------------------
# Create an Excel file with all region information
# 1. Station names
# 2. Station IDs
# 3. Latitude
# 4. Longitude
# 5. First data year
# 6. Last data year
# 7. Number of years
# 8. Residual vlm lower
# 9. Residual vlm mean
# 10. Residual vlm upper
# 11. Residual vlm sources
# -------------------------------------------------
import numpy as np
import os
import pandas as pd
from netCDF4 import Dataset
def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_station_data'] = settings['dir_budget'] + 'tg/station_data.npy'
    settings['fn_region_ensembles']   = settings['dir_budget']+'region_data/region_ensembles.npy'
    settings['fn_output'] = settings['dir_budget'] + 'data_supplement/region_info.xlsx'
    settings['years'] = np.arange(1900, 2019)
    settings['num_ens'] = 5000
    settings['fn_gia_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'
    settings['fn_region_list'] = settings['dir_budget']+'region_data/region_list_final.npy'
    settings['probability'] = Dataset(settings['fn_gia_rsl'], 'r').variables['probability'][:settings['num_ens']]._get_data()
    settings['probability'] = settings['probability'] / settings['probability'].sum()
    region_ensembles = np.load(settings['fn_region_ensembles'],allow_pickle=True)

    station_data = np.load(settings['fn_station_data'], allow_pickle=True).all()
    region_list = np.load(settings['fn_region_list'], allow_pickle=True)
    write_array = pd.ExcelWriter(settings['fn_output'])
    basin_name = ['Subpolar North Atlantic', 'Indian Ocean - South Pacific', 'Subtropical North Atlantic', 'East Pacific', 'South Atlantic', 'Northwest Pacific']
    for basin in range(len(region_list)):
        data_basin = pd.DataFrame(index=np.arange(len(region_list[basin]['list']))+1, columns=['Station names', 'PSMSL IDs', 'Latitude','Longitude','First year','Last year','Number of annual observations','VLM sources','Residual VLM [lower]','Residual VLM [mean]','Residual VLM [upper]'])
        # Mean and STD of region ensembles
        resvlm_trend = np.zeros([len(region_list[basin]['list']), 3]) * np.nan
        resvlm_trend[:,1] = np.sum(settings['probability'] * region_ensembles[basin]['resvlm_ens'], axis=1)
        for region in range(len(resvlm_trend)):
            sort_idx = np.argsort(region_ensembles[basin]['resvlm_ens'][region, :])
            sort_cdf = np.cumsum(settings['probability'][sort_idx])
            resvlm_trend[region, 0] = region_ensembles[basin]['resvlm_ens'][region,sort_idx][np.argmin(np.abs(sort_cdf - 0.05))]
            resvlm_trend[region, 2] = region_ensembles[basin]['resvlm_ens'][region,sort_idx][np.argmin(np.abs(sort_cdf - 0.95))]

        for region in range(len(region_list[basin]['list'])):
            psmsl_name = []
            psmsl_id = []
            acc_tg = np.zeros(len(settings['years']),dtype=bool)
            for stat in region_list[basin]['list'][region]['id']:
                stat_idx = station_data['id']==stat
                psmsl_name.append(station_data['name'][stat_idx][0])
                psmsl_id.append(station_data['id'][stat_idx][0])
                coords = station_data['coords'][stat_idx][0]
                acc_tg[np.isfinite(station_data['height'][stat_idx].flatten())] = True
            data_basin.at[region+1, 'Station names'] = psmsl_name
            data_basin.at[region+1, 'PSMSL IDs'] = psmsl_id
            data_basin.at[region+1, 'Latitude'] = coords[0]
            data_basin.at[region+1, 'Longitude'] = coords[1]
            data_basin.at[region+1, 'Number of annual observations'] = acc_tg.sum()
            data_basin.at[region+1, 'First year'] = settings['years'][acc_tg][0]
            data_basin.at[region+1, 'Last year'] = settings['years'][acc_tg][-1]
            vlm_stations = region_list[basin]['list'][region]['vlm_id']
            if len(vlm_stations) == 0:
                vlm_stations = ['None']
                data_basin.at[region+1, 'Residual VLM [lower]'] = np.nan
                data_basin.at[region+1, 'Residual VLM [mean]'] = np.nan
                data_basin.at[region+1, 'Residual VLM [upper]'] = np.nan
            else:
                data_basin.at[region+1, 'Residual VLM [lower]'] = resvlm_trend[region,0]
                data_basin.at[region+1, 'Residual VLM [mean]'] = resvlm_trend[region,1]
                data_basin.at[region+1, 'Residual VLM [upper]'] = resvlm_trend[region,2]
            data_basin.at[region+1, 'VLM sources'] = vlm_stations
        data_basin['Station names'] = data_basin['Station names'].str.join('; ')
        data_basin['PSMSL IDs'] = data_basin['PSMSL IDs'].astype(str).str[1:-1]
        data_basin['VLM sources'] = data_basin['VLM sources'].str.join('; ')
        data_basin.to_excel(write_array, sheet_name=basin_name[basin])
    write_array.close()
    return

if __name__ == '__main__':
    main()

