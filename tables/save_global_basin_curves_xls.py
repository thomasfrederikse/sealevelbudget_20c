# ---------------------------------------------------------
# Save global and basin-mean sea-level curves as Excel file
# Single file per basin
# ---------------------------------------------------------
import numpy as np
import os
import pandas as pd

def main():
    settings = {}
    settings['dir_data']   = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_stats'] = settings['dir_budget'] + 'results/stats_basin_global.npy'
    settings['fn_region_ensembles']   = settings['dir_budget']+'region_data/region_ensembles.npy'
    settings['fn_output'] = settings['dir_budget'] + 'data_supplement/global_basin_timeseries.xlsx'

    settings['years'] = np.arange(1900,2019)
    stats = np.load(settings['fn_stats'],allow_pickle=True).all()

    write_array = pd.ExcelWriter(settings['fn_output'])

    data_global = pd.DataFrame(data=stats['global']['obs']['tseries'], index=settings['years'], columns=['Observed GMSL [lower]','Observed GMSL [mean]','Observed GMSL [upper]'])
    data_global['Sum of contributors [lower]'] = stats['global']['budget']['tseries'][:,0]
    data_global['Sum of contributors [mean]'] = stats['global']['budget']['tseries'][:,1]
    data_global['Sum of contributors [upper]'] = stats['global']['budget']['tseries'][:,2]
    data_global['Steric [lower]'] = stats['global']['steric']['tseries'][:,0]
    data_global['Steric [mean]'] = stats['global']['steric']['tseries'][:,1]
    data_global['Steric [upper]'] = stats['global']['steric']['tseries'][:,2]
    data_global['Glaciers [lower]'] = stats['global']['grd_glac']['tseries'][:,0]
    data_global['Glaciers [mean]'] = stats['global']['grd_glac']['tseries'][:,1]
    data_global['Glaciers [upper]'] = stats['global']['grd_glac']['tseries'][:,2]
    data_global['Greenland Ice Sheet [lower]'] = stats['global']['grd_GrIS']['tseries'][:,0]
    data_global['Greenland Ice Sheet [mean]'] = stats['global']['grd_GrIS']['tseries'][:,1]
    data_global['Greenland Ice Sheet [upper]'] = stats['global']['grd_GrIS']['tseries'][:,2]
    data_global['Antarctic Ice Sheet [lower]'] = stats['global']['grd_AIS']['tseries'][:,0]
    data_global['Antarctic Ice Sheet [mean]'] = stats['global']['grd_AIS']['tseries'][:,1]
    data_global['Antarctic Ice Sheet [upper]'] = stats['global']['grd_AIS']['tseries'][:,2]
    data_global['Terrestrial Water Storage [lower]'] = stats['global']['grd_tws']['tseries'][:,0]
    data_global['Terrestrial Water Storage [mean]'] = stats['global']['grd_tws']['tseries'][:,1]
    data_global['Terrestrial Water Storage [upper]'] = stats['global']['grd_tws']['tseries'][:,2]
    data_global['Reservoir impoundment [lower]'] = stats['global']['grd_tws_dam']['tseries'][:,0]
    data_global['Reservoir impoundment [mean]'] = stats['global']['grd_tws_dam']['tseries'][:,1]
    data_global['Reservoir impoundment [upper]'] = stats['global']['grd_tws_dam']['tseries'][:,2]
    data_global['Groundwater depletion [lower]'] = stats['global']['grd_tws_gwd']['tseries'][:,0]
    data_global['Groundwater depletion [mean]'] = stats['global']['grd_tws_gwd']['tseries'][:,1]
    data_global['Groundwater depletion [upper]'] = stats['global']['grd_tws_gwd']['tseries'][:,2]
    data_global['Natural TWS [lower]'] = stats['global']['grd_tws_natural']['tseries'][:,0]
    data_global['Natural TWS [mean]'] = stats['global']['grd_tws_natural']['tseries'][:,1]
    data_global['Natural TWS [upper]'] = stats['global']['grd_tws_natural']['tseries'][:,2]
    data_global['Altimetry [lower]'] = stats['global']['altimetry']['tseries'][:,0]
    data_global['Altimetry [mean]'] = stats['global']['altimetry']['tseries'][:,1]
    data_global['Altimetry [upper]'] = stats['global']['altimetry']['tseries'][:,2]
    data_global.to_excel(write_array, sheet_name='Global')

    basin_name = ['Subpolar North Atlantic', 'Indian Ocean - South Pacific', 'Subtropical North Atlantic', 'East Pacific', 'South Atlantic', 'Northwest Pacific']
    for basin in range(len(stats['basin'])):
        data_basin = pd.DataFrame(data=stats['basin'][basin]['obs']['tseries'], index=settings['years'], columns=['Observed basin-mean sea level [lower]', 'Observed basin-mean sea level [mean]', 'Observed basin-mean sea level [upper]'])
        data_basin['Sum of contributors [lower]'] = stats['basin'][basin]['budget']['tseries'][:, 0]
        data_basin['Sum of contributors [mean]'] = stats['basin'][basin]['budget']['tseries'][:, 1]
        data_basin['Sum of contributors [upper]'] = stats['basin'][basin]['budget']['tseries'][:, 2]
        data_basin['Steric [lower]'] = stats['basin'][basin]['steric']['tseries'][:, 0]
        data_basin['Steric [mean]'] = stats['basin'][basin]['steric']['tseries'][:, 1]
        data_basin['Steric [upper]'] = stats['basin'][basin]['steric']['tseries'][:, 2]
        data_basin['Glaciers [lower]'] = stats['basin'][basin]['grd_glac']['tseries'][:, 0]
        data_basin['Glaciers [mean]'] = stats['basin'][basin]['grd_glac']['tseries'][:, 1]
        data_basin['Glaciers [upper]'] = stats['basin'][basin]['grd_glac']['tseries'][:, 2]
        data_basin['Greenland Ice Sheet [lower]'] = stats['basin'][basin]['grd_GrIS']['tseries'][:, 0]
        data_basin['Greenland Ice Sheet [mean]'] = stats['basin'][basin]['grd_GrIS']['tseries'][:, 1]
        data_basin['Greenland Ice Sheet [upper]'] = stats['basin'][basin]['grd_GrIS']['tseries'][:, 2]
        data_basin['Antarctic Ice Sheet [lower]'] = stats['basin'][basin]['grd_AIS']['tseries'][:, 0]
        data_basin['Antarctic Ice Sheet [mean]'] = stats['basin'][basin]['grd_AIS']['tseries'][:, 1]
        data_basin['Antarctic Ice Sheet [upper]'] = stats['basin'][basin]['grd_AIS']['tseries'][:, 2]
        data_basin['Terrestrial Water Storage [lower]'] = stats['basin'][basin]['grd_tws']['tseries'][:, 0]
        data_basin['Terrestrial Water Storage [mean]'] = stats['basin'][basin]['grd_tws']['tseries'][:, 1]
        data_basin['Terrestrial Water Storage [upper]'] = stats['basin'][basin]['grd_tws']['tseries'][:, 2]
        data_basin['GIA [lower]'] = stats['basin'][basin]['gia']['tseries'][:, 0]
        data_basin['GIA [mean]'] = stats['basin'][basin]['gia']['tseries'][:, 1]
        data_basin['GIA [upper]'] = stats['basin'][basin]['gia']['tseries'][:, 2]
        data_basin['Altimetry [lower]'] = stats['basin'][basin]['altimetry']['tseries'][:, 0]
        data_basin['Altimetry [mean]'] = stats['basin'][basin]['altimetry']['tseries'][:, 1]
        data_basin['Altimetry [upper]'] = stats['basin'][basin]['altimetry']['tseries'][:, 2]
        data_basin.to_excel(write_array, sheet_name=basin_name[basin])
    write_array.close()
    return


