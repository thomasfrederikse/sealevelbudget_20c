# -------------------------------------------------
# Download ERA 20c and ERA 5 data for TG regression
# -------------------------------------------------
from ecmwfapi import ECMWFDataServer
import numpy as np
import cdsapi

def main():
    # 1. ERA 20C 1900-1979
    create_date = ''
    years = np.arange(1900,1980)
    mnth = np.arange(1,13)
    for yr in range(len(years)):
        for mn in range(len(mnth)):
            create_date = create_date+(str(years[yr]) + str(mnth[mn]).zfill(2)+'01/')
    create_date = create_date[:-1]

    server = ECMWFDataServer()
    server.retrieve({
        "class": "e2",
        "dataset": "era20c",
        "date": create_date,
        "expver": "1",
        "levtype": "sfc",
        "param": "151.128/165.128/166.128",
        "stream": "moda",
        "type": "an",
        "target": "output",
        'grid': "1/1",

        'format' : "netcdf"
    })

    # ERA 5 1980-2018
    c = cdsapi.Client()

    data = c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'product_type':'monthly_averaged_reanalysis',
            'variable':[
                '10m_u_component_of_wind','10m_v_component_of_wind','mean_sea_level_pressure'
            ],
            'year':[
                '1979','1980','1981',
                '1982','1983','1984',
                '1985','1986','1987',
                '1988','1989','1990',
                '1991','1992','1993',
                '1994','1995','1996',
                '1997','1998','1999',
                '2000','2001','2002',
                '2003','2004','2005',
                '2006','2007','2008',
                '2009','2010','2011',
                '2012','2013','2014',
                '2015','2016','2017',
                '2018'
            ],
            'month':[
                '01','02','03',
                '04','05','06',
                '07','08','09',
                '10','11','12'
            ],
            'time':'00:00',
            'format':'netcdf',
            'grid': "1/1",
    },"/Users/tfrederi/Downloads/ERA5.nc")
    return
