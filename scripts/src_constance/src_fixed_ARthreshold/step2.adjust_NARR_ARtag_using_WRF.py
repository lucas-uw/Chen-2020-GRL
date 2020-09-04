#!/people/chen423/sw/anaconda3/bin/python

import numpy as np
import xarray as xr
import netCDF4 as nc
from calendar import monthrange
import sys

model = sys.argv[1]

def crt_filenames(model, year, month):
    
    ARfile_NARRbase = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/AR_tagged/Gershunov/SERDP6km/WRF_ARtag.%s.Gershunov.%d.%d.nc' % (model, model, year, month)
    
    WRFdir = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/SERDP6km_by_month/' % model
        
    WRFfile_u = WRFdir + 'WRF_%s.6hr.diag.uIVT.%d.%d.nc' % (model, year, month)
    WRFfile_v = WRFdir + 'WRF_%s.6hr.diag.vIVT.%d.%d.nc' % (model, year, month)
    WRFfile_IWV = WRFdir + 'WRF_%s.6hr.diag.IWV.%d.%d.nc' % (model, year, month)
    
    return ARfile_NARRbase, WRFfile_u, WRFfile_v, WRFfile_IWV

def gen_mod_ARtag(model, year, month):
    
    ARfile_NARRbase, WRFfile_u, WRFfile_v, WRFfile_IWV = crt_filenames(model, year, month)
    
    ARtag = xr.open_dataset(ARfile_NARRbase).AR_tag.values
    
    WRF_IVT = np.sqrt(xr.open_dataset(WRFfile_u).uIVT.values**2 + xr.open_dataset(WRFfile_v).vIVT.values**2)
    WRF_IWV = xr.open_dataset(WRFfile_IWV).IWV.values
        
    ARtag_adj = np.zeros(ARtag.shape)
    valid_AR_index = (ARtag>0)*(WRF_IVT>250)*(WRF_IWV>15)
    ARtag_adj[valid_AR_index] = 1
    
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

year = 2003
for month in np.arange(10,13):
    ARtag_adj = gen_mod_ARtag(model, year, month)
    outfile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/AR_tagged/Gershunov/SERDP6km_adj/WRF_ARtag_adj.%s.Gershunov.%d.%d.nc' % (model, model, year, month)
            
    time_string = '%d-%02d' % (year, month)
    save_to_nc(outfile, time_string, ARtag_adj)


for year in np.arange(2004, 2016):
    for month in np.arange(1,13):
        ARtag_adj = gen_mod_ARtag(model, year, month)
        outfile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/AR_tagged/Gershunov/SERDP6km_adj/WRF_ARtag_adj.%s.Gershunov.%d.%d.nc' % (model, model, year, month)
            
        time_string = '%d-%02d' % (year, month)
        save_to_nc(outfile, time_string, ARtag_adj)
