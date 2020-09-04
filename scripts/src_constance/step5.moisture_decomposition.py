#!/people/chen423/sw/anaconda3/bin/python

import numpy as np
import xarray as xr
import scipy.io as sio
import sys

scenario = 'HIST'
year = int(sys.argv[1])
month = int(sys.argv[2])
para_b = int(10)


def compute_moisture_intensity(in_ARtag, in_uIVT, in_ET, ref_mask):
    uIVT_total = in_uIVT[:,0][in_ARtag[:,0]==1].sum()*6000*86400
    ET_total = in_ET[(in_ARtag==1)&(ref_mask==0)].sum()*6000*6000
   
    out_ratio = -9999
    sub_ocean_grids = ((in_ARtag==1)&(ref_mask==0)).sum()
    
    if (ET_total+uIVT_total)!=0:
        out_ratio = ET_total/(ET_total+uIVT_total)
        
    return sub_ocean_grids, ET_total, uIVT_total

reffile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/ref/WRF_latlon.nc'
landmask = xr.open_dataset(reffile).LANDMASK.values[para_b:(450-para_b),para_b:(450-para_b)]

ETdir = '/pic/projects/next_gen_idf/chen423/data/WRF_hist/raw/ET_by_year/'
uIVTdir = '/pic/projects/next_gen_idf/chen423/data/WRF_hist/diagnosic_vars.correct_SST/organized/monthly/'
ARdir = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/AR_tagged/Gershunov/SERDP6km_adj/' % (scenario)

ETfile = ETdir + 'WRF_NARR.%s.SFCEVP.%d.%d.nc' % (scenario, year, month)
uIVTfile = uIVTdir + 'WRF_NARR.%s.uIVT.%d.%d.nc' % (scenario, year, month)

ARfile = ARdir + 'WRF_ARtag_adj.%s.Gershunov.%d.%d.ARp85.nc' % (scenario, year, month)

ETdata = xr.open_dataset(ETfile).SFCEVP.values[:,para_b:(450-para_b),para_b:(450-para_b)]
uIVTdata = xr.open_dataset(uIVTfile).uIVT.values[:,para_b:(450-para_b),para_b:(450-para_b)]
ARtag = xr.open_dataset(ARfile).AR_tag.values[:,para_b:(450-para_b),para_b:(450-para_b)]

nt = ARtag.shape[0]
array_grids = np.zeros(nt)
array_ET = np.zeros(nt)
array_uIVT = np.zeros(nt)
for t in np.arange(nt):
    array_grids[t], array_ET[t], array_uIVT[t] = compute_moisture_intensity(ARtag[t], uIVTdata[t], ETdata[int(np.floor(t/4))], landmask)
    
outfile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/HIST/moisture/ETratio.%s.ARp85.%d.%d.mat' % (scenario, year, month)
sio.savemat(outfile, {'array_grids':array_grids, 'array_ET':array_ET, 'array_uIVT':array_uIVT})
