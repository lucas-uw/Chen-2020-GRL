#!/people/chen423/sw/anaconda3/bin/python

import numpy as np
import xarray as xr
import netCDF4 as nc
from calendar import monthrange
import sys

model = sys.argv[1]
ARmodel_tag = sys.argv[2] # abs or p85

def crt_filenames(model, year, month):
    
    ARfile_NARRbase = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/AR_tagged/Gershunov/SERDP6km/WRF_ARtag.%s.Gershunov.%d.%d.AR%s.nc' % (model, model, year, month, ARmodel_tag)
    
    WRFdir = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/by_month_SERDP6km/' % model
        
    WRFfile = WRFdir + 'WRF_IWV_uvIVT.6hr.%d.%d.nc' % (year, month)
    
    return ARfile_NARRbase, WRFfile

def get_AR_thre_maps():
    file_WRF_IWV = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/HIST/ref/WRF_IWV.1981-2015.6hr.pctl85.nc'
    file_WRF_IVT = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/HIST/ref/WRF_uvIVT.1981-2015.6hr.pctl85.nc'

    fix_IWV = xr.open_dataset(file_WRF_IWV).IWV.values[0]
    fix_IVT = xr.open_dataset(file_WRF_IVT).uvIVT.values[0]

    return fix_IWV, fix_IVT



def gen_mod_ARtag(model, year, month, ARmodel_tag):
    
    ARfile_NARRbase, WRFfile = crt_filenames(model, year, month)
    
    ARtag = xr.open_dataset(ARfile_NARRbase).AR_tag.values
    
    WRF_IVT = xr.open_dataset(WRFfile).uvIVT.values
    WRF_IWV = xr.open_dataset(WRFfile).IWV.values
        
    ARtag_adj = np.zeros(ARtag.shape)

    if ARmodel_tag=='abs':
        valid_AR_index = (ARtag>0)*(WRF_IVT>250)*(WRF_IWV>15)
        ARtag_adj[valid_AR_index] = 1
    elif ARmodel_tag=='p85':
       for t in np.arange(ARtag.shape[0]):
            ARtag_adj[t] = (ARtag[t]>0)*(WRF_IVT[t]>map_IVT_thre)*(WRF_IWV[t]>map_IWV_thre)
    
    return ARtag_adj

def save_to_nc(outfile, stime, indata):
    outgroup = nc.Dataset(outfile, 'w', format='NETCDF4')

    # dimension
    tdim = outgroup.createDimension('time', None)
    xdim = outgroup.createDimension('x', 450)
    ydim = outgroup.createDimension('y', 450)

    # vars
    tvar = outgroup.createVariable('time', 'i4', ('time'))
    tvar.units = 'Hours since %s-01 00:00' % stime
    tvar[:] = (np.arange(indata.shape[0])*6).astype('int32')
    xvar = outgroup.createVariable('x', 'i4', ('x'))
    xvar[:] = np.arange(450)
    yvar = outgroup.createVariable('y', 'i4', ('y'))
    yvar[:] = np.arange(450)

    var_lat = outgroup.createVariable('lat', 'f4', ('y', 'x'))
    var_lat.long_name = 'Latitude'
    var_lat.standard_name = 'latitude'
    var_lat.units = 'degrees_north'
    var_lon = outgroup.createVariable('lon', 'f4', ('y', 'x'))
    var_lon.long_name = 'Longitude'
    var_lon.standard_name = 'longitude'
    var_lon.units = 'degrees_east'
    var_tag = outgroup.createVariable('AR_tag', 'i4', ('time','y', 'x'))
    var_tag.long_name = 'AR tag based on NARR-Gershunov, adjusted using WRF data'
    var_tag.coordinates = 'lon lat'
    
    var_lat[:] = wrf_lat
    var_lon[:] = wrf_lon
    var_tag[:] = indata

    outgroup.close()

reffile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/ref/SERDP6km_latlon.nc'
wrf_lat = xr.open_dataset(reffile).XLAT_M.values
wrf_lon = xr.open_dataset(reffile).XLONG_M.values

map_IWV_thre, map_IVT_thre = get_AR_thre_maps()

year = 2003
for month in np.arange(10,13):
    print(year, month)
    ARtag_adj = gen_mod_ARtag(model, year, month, ARmodel_tag)
    outfile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/AR_tagged/Gershunov/SERDP6km_adj/WRF_ARtag_adj.%s.Gershunov.%d.%d.AR%s.nc' % (model, model, year, month, ARmodel_tag)
            
    time_string = '%d-%02d' % (year, month)
    save_to_nc(outfile, time_string, ARtag_adj)


for year in np.arange(2004, 2016):
#for year in np.arange(2013, 2016):
    for month in np.arange(1,13):
        print(year, month)
        ARtag_adj = gen_mod_ARtag(model, year, month, ARmodel_tag)
        outfile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/AR_tagged/Gershunov/SERDP6km_adj/WRF_ARtag_adj.%s.Gershunov.%d.%d.AR%s.nc' % (model, model, year, month, ARmodel_tag)
            
        time_string = '%d-%02d' % (year, month)
        save_to_nc(outfile, time_string, ARtag_adj)
