#!/people/chen423/sw/anaconda3/bin/python

import numpy as np
import xarray as xr
import scipy.io as sio
import sys

scenario = 'HIST'
year = int(sys.argv[1])
month = int(sys.argv[2])
para_b = int(10)


def compute_moisture_intensity(in_uIVT, in_vIVT, in_ET, ref_mask):
    # refmask is landmak, land is 1, ocean is 0
    
    # note:
    #      uIVT: into the domain
    #      vIVT: bottom: into the domain;   top: away from the domain
    
    uIVT_total = in_uIVT[:,0].sum()*6000*86400
    vIVT_sub_bottom = in_vIVT[0,:][ref_mask[0,:]==0].sum()*6000*86400
    vIVT_sub_top = in_vIVT[(450-2*para_b-1),:][ref_mask[(450-2*para_b-1),:]==0].sum()*6000*86400
    ET_total = in_ET[(ref_mask==0)].sum()*6000*6000

    return ET_total, uIVT_total, vIVT_sub_bottom, vIVT_sub_top

reffile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/ref/WRF_latlon.nc'
landmask = xr.open_dataset(reffile).LANDMASK.values[para_b:(450-para_b),para_b:(450-para_b)]

ETdir = '/pic/projects/next_gen_idf/chen423/data/WRF_hist/raw/ET_by_year/'
uIVTdir = '/pic/projects/next_gen_idf/chen423/data/WRF_hist/diagnosic_vars.correct_SST/organized/monthly/'
vIVTdir = '/pic/projects/next_gen_idf/chen423/data/WRF_hist/diagnosic_vars.correct_SST/organized/monthly/'

ETfile = ETdir + 'WRF_NARR.%s.SFCEVP.%d.%d.nc' % (scenario, year, month)
uIVTfile = uIVTdir + 'WRF_NARR.%s.uIVT.%d.%d.nc' % (scenario, year, month)
vIVTfile = vIVTdir + 'WRF_NARR.%s.vIVT.%d.%d.nc' % (scenario, year, month)


ETdata = xr.open_dataset(ETfile).SFCEVP.values[:,para_b:(450-para_b),para_b:(450-para_b)]
uIVTdata = xr.open_dataset(uIVTfile).uIVT.values[:,para_b:(450-para_b),para_b:(450-para_b)]
vIVTdata = xr.open_dataset(vIVTfile).vIVT.values[:,para_b:(450-para_b),para_b:(450-para_b)]

nt = uIVTdata.shape[0]
array_ET = np.zeros(nt)
array_uIVT = np.zeros(nt)
array_vIVT_bottom = np.zeros(nt)
array_vIVT_top = np.zeros(nt)
for t in np.arange(nt):
    array_ET[t], array_uIVT[t], array_vIVT_bottom[t], array_vIVT_top[t] = compute_moisture_intensity(uIVTdata[t], vIVTdata[t], ETdata[int(np.floor(t/4))], landmask)
    
outfile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/HIST/moisture/ETratio.%s.full_ocean.%d.%d.mat' % (scenario, year, month)
sio.savemat(outfile, {'array_ET':array_ET, 'array_uIVT':array_uIVT,
                      'array_vIVT_bottom':array_vIVT_bottom, 'array_vIVT_top':array_vIVT_top})
