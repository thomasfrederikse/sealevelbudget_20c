# Save GMSL and its components in a MEAsUREs netcdf file
import numpy as np
import os
from netCDF4 import Dataset, date2num
import datetime as dt

def main():
    settings = {}
    settings['dir_data']   = os.getenv('HOME') + '/Data/'
    settings['dir_budget'] = settings['dir_data'] + 'Budget_20c/'
    settings['fn_stats'] = settings['dir_budget'] + 'results/stats_basin_global.npy'
    settings['fn_region_ensembles']   = settings['dir_budget']+'region_data/region_ensembles.npy'
    settings['fn_output'] = settings['dir_budget'] + 'supplement/global_timeseries_measures.nc'

    settings['years'] = np.arange(1900,2019)
    stats = np.load(settings['fn_stats'],allow_pickle=True).all()

    date_array = []
    for year in settings['years']:
            date_array.append(dt.datetime(year, 6, 15))

    fh = Dataset(settings['fn_output'],'w')
    fh.createDimension('time', len(settings['years']))
    # Global attributes
    attributeDictionary = {'title': 'Global sea-level changes and contributing processes over 1900-2018',
                           'summary': 'This file contains reconstructed global-mean sea level evolution and the estimated contributing processes over 1900-2018. Reconstructed sea level is based on annual-mean tide-gauge observations and uses the virtual-station method to aggregrate the individual observations into a global estimate. The contributing processes consist of thermosteric changes, glacier mass changes, mass changes of the Greenland and Antarctic Ice Sheet, and terrestrial water storage changes. The glacier, ice sheet, and terrestrial water storage are estimated by combining GRACE observations (2003-2018) with long-term estimates from in-situ observations and models. Steric estimates are based on in-situ temperature profiles. The upper- and lower bound represent the 5 and 95 percent confidence level. The numbers are equal to the ones presented in Frederikse et al. The causes of sea-level rise since 1900, Nature, 2020, reformatted to meet the specifications of the JPL PO.DAAC',
                           'Conventions': 'CF-1.6, ACDD-1.3',
                           'history': 'This version provides the data as presented in Frederikse et al. 2020.',
                           'source': 'Frederikse et al. The causes of sea-level rise since 1900, Nature, 2020 https://doi.org/10.1038/s41586-020-2591-3',
                           'standard_name_vocabulary': 'NetCDF Climate and Forecast (CF) Metadata Convention',
                           'acknowledgement': 'This research was carried out by the Jet Propulsion Laboratory, managed by the California Institute of Technology under a contract with the National Aeronautics and Space Administration.',
                           'license': 'Public Domain',
                           'product_version': '1.0',
                           'creator_name': 'Thomas Frederikse',
                           'creator_email': 'thomas.frederikse@jpl.nasa.gov',
                           'creator_institution': 'NASA Jet Propulsion Laboratory (JPL)',
                           'institution': 'NASA Jet Propulsion Laboratory (JPL)',
                           'project': 'NASA sea-level change science team (N-SLCT)',
                           'program': 'NASA sea-level change science team (N-SLCT)',
                           'contributor_name': 'Thomas Frederikse, Felix Landerer, Lambert Caron, Surendra Adhikari, David Parkes, Vincent Humphrey, Sönke Dangendorf, Peter Hogarth, Laure Zanna, Lijing Cheng, Yun-Hao Wu',
                           'contributor_role': 'main author,author,author,author,author,author,author,author,author,author,author',
                           'publisher_name': 'Physical Oceanography Distributed Active Archive Center (PO.DAAC)',
                           'publisher_email': 'podaac@podaac.jpl.nasa.gov',
                           'publisher_type': 'group',
                           'publisher_institution': 'NASA Jet Propulsion Laboratory (JPL)',
                           'time_coverage_start': '1900-01-01',
                           'time_coverage_end': '2019-01-01',
                           'time_coverage_duration': '119 years',
                           'time_coverage_resolution': 'Annual',
                           'date_created': 'July 7, 2020',
                           'date_modified': '',
                           'date_issued': 'PO.DAAC'}
    fh.setncatts(attributeDictionary)

    # Time
    fh_time = fh.createVariable('time', 'i4', ('time'), zlib=True, complevel=4)
    fh_time[:] = date2num(date_array, 'days since 1900-01-01 00:00:00 UTC')
    fh_time.standard_name='time'
    fh_time.valid_min = -1e6
    fh_time.valid_max = 1e6
    fh_time.missing_value = -1e7
    fh_time.fill_value = -1e7
    fh_time.units = 'days since 1900-01-01 00:00:00 UTC'
    fh_time.long_name = 'time'
    # Observed
    fh_slm_mean = fh.createVariable('global_average_sea_level_change', 'f4', ('time'), zlib=True, complevel=4)
    fh_slm_mean[:] = stats['global']['obs']['tseries'][:,1]
    fh_slm_mean.units = 'mm'
    fh_slm_mean.long_name = 'Observed global-average sea level (mean value)'
    fh_slm_mean.valid_min = -1e6
    fh_slm_mean.valid_max = 1e6
    fh_slm_mean.comment = "Value is relative to 2002-2019 baseline"
    fh_slm_mean.standard_name='global_average_sea_level_change'
    fh_slm_mean.missing_value = -1e7
    fh_slm_mean.fill_value = -1e7
    fh_slm_mean.coverage_content_type = 'physicalMeasurement'

    fh_slm_upper = fh.createVariable('global_average_sea_level_change_upper', 'f4', ('time'), zlib=True, complevel=4)
    fh_slm_upper[:] = stats['global']['obs']['tseries'][:,2]
    fh_slm_upper.units = 'mm'
    fh_slm_upper.long_name = 'Observed global-average sea level (upper bound)'
    fh_slm_upper.standard_name='global_average_sea_level_change'
    fh_slm_upper.valid_min = -1e6
    fh_slm_upper.valid_max = 1e6
    fh_slm_upper.comment = "Value is relative to 2002-2019 baseline"
    fh_slm_upper.missing_value = -1e7
    fh_slm_upper.fill_value = -1e7
    fh_slm_upper.coverage_content_type = 'physicalMeasurement'

    fh_slm_lower = fh.createVariable('global_average_sea_level_change_lower', 'f4', ('time'), zlib=True, complevel=4)
    fh_slm_lower[:] = stats['global']['obs']['tseries'][:,1]
    fh_slm_lower.units = 'mm'
    fh_slm_lower.long_name = 'Observed global-mean sea level (lower bound)'
    fh_slm_lower.standard_name='global_average_sea_level_change'
    fh_slm_lower.valid_min = -1e6
    fh_slm_lower.valid_max = 1e6
    fh_slm_lower.comment = "Value is relative to 2002-2019 baseline"
    fh_slm_lower.missing_value = -1e7
    fh_slm_lower.fill_value = -1e7
    fh_slm_lower.coverage_content_type = 'physicalMeasurement'

    # Glaciers
    fh_glac_mean = fh.createVariable('glac_mean', 'f4', ('time'), zlib=True, complevel=4)
    fh_glac_mean[:] = stats['global']['grd_glac']['tseries'][:,1]
    fh_glac_mean.units = 'mm'
    fh_glac_mean.long_name = 'Glacier contribution (mean value)'
    fh_glac_mean.valid_min = -1e6
    fh_glac_mean.valid_max = 1e6
    fh_glac_mean.comment = "Value is relative to 2002-2019 baseline. The glacier term does not include peripheral glaciers of both ice sheets"
    fh_glac_mean.missing_value = -1e7
    fh_glac_mean.fill_value = -1e7
    fh_glac_mean.standard_name='global_average_sea_level_change'
    fh_glac_mean.coverage_content_type = 'physicalMeasurement'

    fh_glac_upper = fh.createVariable('glac_upper', 'f4', ('time'), zlib=True, complevel=4)
    fh_glac_upper[:] = stats['global']['grd_glac']['tseries'][:,2]
    fh_glac_upper.units = 'mm'
    fh_glac_upper.long_name = 'Glacier contribution (upper bound)'
    fh_glac_upper.valid_min = -1e6
    fh_glac_upper.valid_max = 1e6
    fh_glac_upper.comment = "Value is relative to 2002-2019 baseline. The glacier term does not include peripheral glaciers of both ice sheets"
    fh_glac_upper.missing_value = -1e7
    fh_glac_upper.fill_value = -1e7
    fh_glac_upper.standard_name='global_average_sea_level_change'
    fh_glac_upper.coverage_content_type = 'physicalMeasurement'

    fh_glac_lower = fh.createVariable('glac_lower', 'f4', ('time'), zlib=True, complevel=4)
    fh_glac_lower[:] = stats['global']['grd_glac']['tseries'][:,0]
    fh_glac_lower.units = 'mm'
    fh_glac_lower.long_name = 'Glacier contribution (lower bound)'
    fh_glac_lower.valid_min = -1e6
    fh_glac_lower.valid_max = 1e6
    fh_glac_lower.comment = "Value is relative to 2002-2019 baseline. The glacier term does not include peripheral glaciers of both ice sheets"
    fh_glac_lower.missing_value = -1e7
    fh_glac_lower.fill_value = -1e7
    fh_glac_lower.standard_name='global_average_sea_level_change'
    fh_glac_lower.coverage_content_type = 'physicalMeasurement'

    # Greenland
    fh_GrIS_mean = fh.createVariable('GrIS_mean', 'f4', ('time'), zlib=True, complevel=4)
    fh_GrIS_mean[:] = stats['global']['grd_GrIS']['tseries'][:,1]
    fh_GrIS_mean.units = 'mm'
    fh_GrIS_mean.long_name = 'Greenland Ice Sheet contribution (mean value)'
    #fh_GrIS_mean.coordinates = "time"
    fh_GrIS_mean.valid_min = -1e6
    fh_GrIS_mean.valid_max = 1e6
    fh_GrIS_mean.comment= "Value is relative to 2002-2019 baseline. This term includes glaciers in the Greenland periphery"
    fh_GrIS_mean.missing_value = -1e7
    fh_GrIS_mean.fill_value = -1e7
    fh_GrIS_mean.standard_name='global_average_sea_level_change'
    fh_GrIS_mean.coverage_content_type = 'physicalMeasurement'

    fh_GrIS_upper = fh.createVariable('GrIS_upper', 'f4', ('time'), zlib=True, complevel=4)
    fh_GrIS_upper[:] = stats['global']['grd_GrIS']['tseries'][:,2]
    fh_GrIS_upper.units = 'mm'
    fh_GrIS_upper.long_name = 'Greenland Ice Sheet contribution (upper bound)'
    #fh_GrIS_upper.coordinates = "time"
    fh_GrIS_upper.valid_min = -1e6
    fh_GrIS_upper.valid_max = 1e6
    fh_GrIS_upper.comment = "Value is relative to 2002-2019 baseline. This term includes glaciers in the Greenland periphery"
    fh_GrIS_upper.missing_value = -1e7
    fh_GrIS_upper.fill_value = -1e7
    fh_GrIS_upper.standard_name='global_average_sea_level_change'
    fh_GrIS_upper.coverage_content_type = 'physicalMeasurement'

    fh_GrIS_lower = fh.createVariable('GrIS_lower', 'f4', ('time'), zlib=True, complevel=4)
    fh_GrIS_lower[:] = stats['global']['grd_GrIS']['tseries'][:,0]
    fh_GrIS_lower.units = 'mm'
    fh_GrIS_lower.long_name = 'Greenland Ice Sheet contribution (lower bound)'
    fh_GrIS_lower.valid_min = -1e6
    fh_GrIS_lower.valid_max = 1e6
    fh_GrIS_lower.comment = "Value is relative to 2002-2019 baseline. This term includes glaciers in the Greenland periphery"
    fh_GrIS_lower.missing_value = -1e7
    fh_GrIS_lower.fill_value = -1e7
    fh_GrIS_lower.standard_name='global_average_sea_level_change'
    fh_GrIS_lower.coverage_content_type = 'physicalMeasurement'

    # Antarctica
    fh_AIS_mean = fh.createVariable('AIS_mean', 'f4', ('time'), zlib=True, complevel=4)
    fh_AIS_mean[:] = stats['global']['grd_AIS']['tseries'][:,1]
    fh_AIS_mean.units = 'mm'
    fh_AIS_mean.long_name = 'Antarctic Ice Sheet contribution (mean value)'
    fh_AIS_mean.valid_min = -1e6
    fh_AIS_mean.valid_max = 1e6
    fh_AIS_mean.comment= "Value is relative to 2002-2019 baseline."
    fh_AIS_mean.missing_value = -1e7
    fh_AIS_mean.fill_value = -1e7
    fh_AIS_mean.standard_name='global_average_sea_level_change'
    fh_AIS_mean.coverage_content_type = 'physicalMeasurement'

    fh_AIS_upper = fh.createVariable('AIS_upper', 'f4', ('time'), zlib=True, complevel=4)
    fh_AIS_upper[:] = stats['global']['grd_AIS']['tseries'][:,2]
    fh_AIS_upper.units = 'mm'
    fh_AIS_upper.long_name = 'Antarctic Ice Sheet contribution (upper bound)'
    fh_AIS_upper.valid_min = -1e6
    fh_AIS_upper.valid_max = 1e6
    fh_AIS_upper.comment = "Value is relative to 2002-2019 baseline."
    fh_AIS_upper.missing_value = -1e7
    fh_AIS_upper.fill_value = -1e7
    fh_AIS_upper.standard_name='global_average_sea_level_change'
    fh_AIS_upper.coverage_content_type = 'physicalMeasurement'

    fh_AIS_lower = fh.createVariable('AIS_lower', 'f4', ('time'), zlib=True, complevel=4)
    fh_AIS_lower[:] = stats['global']['grd_AIS']['tseries'][:,0]
    fh_AIS_lower.units = 'mm'
    fh_AIS_lower.long_name = 'Antarctic Ice Sheet contribution (lower bound)'
    fh_AIS_lower.valid_min = -1e6
    fh_AIS_lower.valid_max = 1e6
    fh_AIS_lower.comment = "Value is relative to 2002-2019 baseline."
    fh_AIS_lower.missing_value = -1e7
    fh_AIS_lower.fill_value = -1e7
    fh_AIS_lower.standard_name='global_average_sea_level_change'
    fh_AIS_lower.coverage_content_type = 'physicalMeasurement'

    # TWS
    fh_tws_mean = fh.createVariable('tws_mean', 'f4', ('time'), zlib=True, complevel=4)
    fh_tws_mean[:] = stats['global']['grd_tws']['tseries'][:,1]
    fh_tws_mean.units = 'mm'
    fh_tws_mean.long_name = 'Terrestrial water storage contribution (mean value)'
    fh_tws_mean.valid_min = -1e6
    fh_tws_mean.valid_max = 1e6
    fh_tws_mean.comment= "Value is relative to 2002-2019 baseline. This term includes the effects for groundwater depletion, reservoir impoundment and non-anthropogenic tws changes"
    fh_tws_mean.missing_value = -1e7
    fh_tws_mean.fill_value = -1e7
    fh_tws_mean.standard_name='global_average_sea_level_change'
    fh_tws_mean.coverage_content_type = 'physicalMeasurement'

    fh_tws_upper = fh.createVariable('tws_upper', 'f4', ('time'), zlib=True, complevel=4)
    fh_tws_upper[:] = stats['global']['grd_tws']['tseries'][:,2]
    fh_tws_upper.units = 'mm'
    fh_tws_upper.long_name = 'Terrestrial water storage contribution (upper bound)'
    fh_tws_upper.valid_min = -1e6
    fh_tws_upper.valid_max = 1e6
    fh_tws_upper.comment= "Value is relative to 2002-2019 baseline. This term includes the effects for groundwater depletion, reservoir impoundment and non-anthropogenic tws changes"
    fh_tws_upper.missing_value = -1e7
    fh_tws_upper.fill_value = -1e7
    fh_tws_upper.standard_name='global_average_sea_level_change'
    fh_tws_upper.coverage_content_type = 'physicalMeasurement'

    fh_tws_lower = fh.createVariable('tws_lower', 'f4', ('time'), zlib=True, complevel=4)
    fh_tws_lower[:] = stats['global']['grd_tws']['tseries'][:,0]
    fh_tws_lower.units = 'mm'
    fh_tws_lower.long_name = 'Terrestrial water storage contribution (lower bound)'
    fh_tws_lower.valid_min = -1e6
    fh_tws_lower.valid_max = 1e6
    fh_tws_lower.comment= "Value is relative to 2002-2019 baseline. This term includes the effects for groundwater depletion, reservoir impoundment and non-anthropogenic tws changes"
    fh_tws_lower.missing_value = -1e7
    fh_tws_lower.fill_value = -1e7
    fh_tws_lower.standard_name='global_average_sea_level_change'
    fh_tws_lower.coverage_content_type = 'physicalMeasurement'

    # Steric
    fh_steric_mean = fh.createVariable('global_average_thermosteric_sea_level_change', 'f4', ('time'), zlib=True, complevel=4)
    fh_steric_mean[:] = stats['global']['steric']['tseries'][:,1]
    fh_steric_mean.units = 'mm'
    fh_steric_mean.long_name = 'Global thermosteric contribution (mean value)'
    fh_steric_mean.valid_min = -1e6
    fh_steric_mean.valid_max = 1e6
    fh_steric_mean.comment= "Value is relative to 2002-2019 baseline"
    fh_steric_mean.missing_value = -1e7
    fh_steric_mean.fill_value = -1e7
    fh_steric_mean.standard_name='global_average_thermosteric_sea_level_change'
    fh_steric_mean.coverage_content_type = 'physicalMeasurement'

    fh_steric_upper = fh.createVariable('global_average_thermosteric_sea_level_change_upper', 'f4', ('time'), zlib=True, complevel=4)
    fh_steric_upper[:] = stats['global']['steric']['tseries'][:,2]
    fh_steric_upper.units = 'mm'
    fh_steric_upper.long_name = 'Global thermosteric contribution (upper bound)'
    fh_steric_upper.valid_min = -1e6
    fh_steric_upper.valid_max = 1e6
    fh_steric_upper.comment= "Value is relative to 2002-2019 baseline"
    fh_steric_upper.missing_value = -1e7
    fh_steric_upper.fill_value = -1e7
    fh_steric_upper.standard_name='global_average_thermosteric_sea_level_change'
    fh_steric_upper.coverage_content_type = 'physicalMeasurement'

    fh_steric_lower = fh.createVariable('global_average_thermosteric_sea_level_change_lower', 'f4', ('time'), zlib=True, complevel=4)
    fh_steric_lower[:] = stats['global']['steric']['tseries'][:,0]
    fh_steric_lower.units = 'mm'
    fh_steric_lower.long_name = 'Global thermosteric contribution (lower bound)'
    fh_steric_lower.valid_min = -1e6
    fh_steric_lower.valid_max = 1e6
    fh_steric_lower.comment= "Value is relative to 2002-2019 baseline"
    fh_steric_lower.missing_value = -1e7
    fh_steric_lower.fill_value = -1e7
    fh_steric_lower.standard_name='global_average_thermosteric_sea_level_change'
    fh_steric_lower.coverage_content_type = 'physicalMeasurement'

    # Sum of processes
    fh_budget_mean = fh.createVariable('sum_of_contrib_processes_mean', 'f4', ('time'), zlib=True, complevel=4)
    fh_budget_mean[:] = stats['global']['budget']['tseries'][:,1]
    fh_budget_mean.units = 'mm'
    fh_budget_mean.long_name = 'Sum of contributing processes (mean value)'
    fh_budget_mean.valid_min = -1e6
    fh_budget_mean.valid_max = 1e6
    fh_budget_mean.comment= "Value is relative to 2002-2019 baseline"
    fh_budget_mean.missing_value = -1e7
    fh_budget_mean.fill_value = -1e7
    fh_budget_mean.coverage_content_type = 'physicalMeasurement'

    fh_budget_upper = fh.createVariable('sum_of_contrib_processes_upper', 'f4', ('time'), zlib=True, complevel=4)
    fh_budget_upper[:] = stats['global']['budget']['tseries'][:,2]
    fh_budget_upper.units = 'mm'
    fh_budget_upper.long_name = 'Sum of contributing processes (upper bound)'
    fh_budget_upper.valid_min = -1e6
    fh_budget_upper.valid_max = 1e6
    fh_budget_upper.comment= "Value is relative to 2002-2019 baseline"
    fh_budget_upper.missing_value = -1e7
    fh_budget_upper.fill_value = -1e7
    fh_budget_upper.coverage_content_type = 'physicalMeasurement'

    fh_budget_lower = fh.createVariable('sum_of_contrib_processes_lower', 'f4', ('time'), zlib=True, complevel=4)
    fh_budget_lower[:] = stats['global']['budget']['tseries'][:,0]
    fh_budget_lower.units = 'mm'
    fh_budget_lower.long_name = 'Sum of contributing processes (lower bound)'
    fh_budget_lower.valid_min = -1e6
    fh_budget_lower.valid_max = 1e6
    fh_budget_lower.comment= "Value is relative to 2002-2019 baseline"
    fh_budget_lower.missing_value = -1e7
    fh_budget_lower.fill_value = -1e7
    fh_budget_lower.coverage_content_type = 'physicalMeasurement'


    fh.close()

    # write_array = pd.ExcelWriter(settings['fn_output'])
    #
    # data_global = pd.DataFrame(data=stats['global']['obs']['tseries'], index=settings['years'], columns=['Observed GMSL [lower]','Observed GMSL [mean]','Observed GMSL [upper]'])
    # data_global['Sum of contributors [lower]'] = stats['global']['budget']['tseries'][:,0]
    # data_global['Sum of contributors [mean]'] = stats['global']['budget']['tseries'][:,1]
    # data_global['Sum of contributors [upper]'] = stats['global']['budget']['tseries'][:,2]
    # data_global['Steric [lower]'] = stats['global']['steric']['tseries'][:,0]
    # data_global['Steric [mean]'] = stats['global']['steric']['tseries'][:,1]
    # data_global['Steric [upper]'] = stats['global']['steric']['tseries'][:,2]
    # data_global['Glaciers [lower]'] = stats['global']['grd_glac']['tseries'][:,0]
    # data_global['Glaciers [mean]'] = stats['global']['grd_glac']['tseries'][:,1]
    # data_global['Glaciers [upper]'] = stats['global']['grd_glac']['tseries'][:,2]
    # data_global['Greenland Ice Sheet [lower]'] = stats['global']['grd_GrIS']['tseries'][:,0]
    # data_global['Greenland Ice Sheet [mean]'] = stats['global']['grd_GrIS']['tseries'][:,1]
    # data_global['Greenland Ice Sheet [upper]'] = stats['global']['grd_GrIS']['tseries'][:,2]
    # data_global['Antarctic Ice Sheet [lower]'] = stats['global']['grd_AIS']['tseries'][:,0]
    # data_global['Antarctic Ice Sheet [mean]'] = stats['global']['grd_AIS']['tseries'][:,1]
    # data_global['Antarctic Ice Sheet [upper]'] = stats['global']['grd_AIS']['tseries'][:,2]
    # data_global['Terrestrial Water Storage [lower]'] = stats['global']['grd_tws']['tseries'][:,0]
    # data_global['Terrestrial Water Storage [mean]'] = stats['global']['grd_tws']['tseries'][:,1]
    # data_global['Terrestrial Water Storage [upper]'] = stats['global']['grd_tws']['tseries'][:,2]
    # data_global['Reservoir impoundment [lower]'] = stats['global']['grd_tws_dam']['tseries'][:,0]
    # data_global['Reservoir impoundment [mean]'] = stats['global']['grd_tws_dam']['tseries'][:,1]
    # data_global['Reservoir impoundment [upper]'] = stats['global']['grd_tws_dam']['tseries'][:,2]
    # data_global['Groundwater depletion [lower]'] = stats['global']['grd_tws_gwd']['tseries'][:,0]
    # data_global['Groundwater depletion [mean]'] = stats['global']['grd_tws_gwd']['tseries'][:,1]
    # data_global['Groundwater depletion [upper]'] = stats['global']['grd_tws_gwd']['tseries'][:,2]
    # data_global['Natural TWS [lower]'] = stats['global']['grd_tws_natural']['tseries'][:,0]
    # data_global['Natural TWS [mean]'] = stats['global']['grd_tws_natural']['tseries'][:,1]
    # data_global['Natural TWS [upper]'] = stats['global']['grd_tws_natural']['tseries'][:,2]
    # data_global['Altimetry [lower]'] = stats['global']['altimetry']['tseries'][:,0]
    # data_global['Altimetry [mean]'] = stats['global']['altimetry']['tseries'][:,1]
    # data_global['Altimetry [upper]'] = stats['global']['altimetry']['tseries'][:,2]
    # data_global.to_excel(write_array, sheet_name='Global')



    return