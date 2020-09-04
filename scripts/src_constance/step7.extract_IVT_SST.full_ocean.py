#!/people/chen423/sw/anaconda3/bin/python
#######!/share/apps/python/anaconda3.6/bin/python

import numpy as np
import xarray as xr
import scipy.io as sio
import pandas as pd
import calendar

rootdir = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/'

def crt_filenames(model, year, month):
    
    WRFdir = rootdir + '%s/by_month_SERDP6km/' % model
    WRFfile = WRFdir + 'WRF_IWV_uvIVT.6hr.%d.%d.nc' % (year, month)
    
    return WRFfile

def get_intensity_data(model, year, month):
    
    WRFfile = crt_filenames(model, year, month)
    
    WRF_IVT = xr.open_dataset(WRFfile).uvIVT.values
    WRF_IWV = xr.open_dataset(WRFfile).IWV.values
    
    return WRF_IVT, WRF_IWV

def compute_6hrly_SST(in_SST):
    out_SSTmean = in_SST[ocean_mask==1].mean()
    
    return out_SSTmean

def compute_6hrly_intensity(in_intensity):
    out_intensity = in_intensity[ocean_mask==1].mean()
    
    return out_intensity


def compute_stats(year, month):
    
    file_SSTHIST = rootdir + 'HIST/SST/NARR_TS.SERDP6km.6hourly.%d.%d.nc' % (year, month)
    
    file_SSTfix = rootdir + 'HIST/SST/NARR_TS.SERDP6km.2000.10.01.00.nc'
    
    
    SST_HIST = xr.open_dataset(file_SSTHIST).var11.values
    IVT_HIST, IWV_HIST = get_intensity_data('HIST', year, month)

    SST_fSST = xr.open_dataset(file_SSTfix).var11.values[0]
    IVT_fSST, IWV_fSST = get_intensity_data('fSST', year, month)

    
    # compute various stats
    nt = SST_HIST.shape[0]
    stat_SSTmean = np.zeros((2,nt))-9999
    stat_IVT = np.zeros((2,nt))-9999
    stat_IWV = np.zeros((2,nt))-9999
    stat_IVTf = np.zeros((2,nt))-9999
    stat_IWVf = np.zeros((2,nt))-9999

    for t in np.arange(nt):
        #print('t= ', t)
        stat_SSTmean[0,t] = compute_6hrly_SST(SST_HIST[t])
        stat_IVT[0,t] = compute_6hrly_intensity(IVT_HIST[t])
        stat_IWV[0,t] = compute_6hrly_intensity(IWV_HIST[t])
        stat_IWVf[0,t] = IWV_HIST[t].mean()
        stat_IVTf[0,t] = IVT_HIST[t].mean()

        stat_SSTmean[1,t] = compute_6hrly_SST(SST_fSST)
        stat_IVT[1,t] = compute_6hrly_intensity(IVT_fSST[t])
        stat_IWV[1,t] = compute_6hrly_intensity(IWV_fSST[t])
        stat_IWVf[1,t] = IWV_fSST[t].mean()
        stat_IVTf[1,t] = IVT_fSST[t].mean()
            
            
    return stat_SSTmean, stat_IVT, stat_IWV, stat_IVTf, stat_IWVf



reffile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/ref/SERDP6km.dist_to_coastal.nc'
dist_to_coast = xr.open_dataset(reffile).dist_to_coast.values
dist_to_coast[dist_to_coast==9999] = 0
ocean_mask = np.zeros((450,450))
ocean_mask[dist_to_coast==0] = 1


## part 1
stats_SSTmean = np.zeros((2,17532))-9999
stats_IVT = np.zeros((2,17532))-9999 # over ocean
stats_IWV = np.zeros((2,17532))-9999 # over ocean
stats_IVTf = np.zeros((2,17532))-9999 # all
stats_IWVf = np.zeros((2,17532))-9999 # all
bg_year = np.zeros(17532)-9999
bg_month = np.zeros(17532)-9999

sindex = -31*4
eindex = 0

year = 2003
print('working on ', year)
for month in np.arange(10,13):
    tmp_SSTmean, tmp_IVT, tmp_IWV, tmp_IVTf, tmp_IWVf = compute_stats(year, month)
    
    sindex = eindex
    eindex = eindex + calendar.monthrange(year, month)[1]*4
    stats_SSTmean[:, sindex:eindex] = tmp_SSTmean
    stats_IVT[:, sindex:eindex] = tmp_IVT
    stats_IWV[:, sindex:eindex] = tmp_IWV
    stats_IVTf[:, sindex:eindex] = tmp_IVTf
    stats_IWVf[:, sindex:eindex] = tmp_IWVf
    bg_year[sindex:eindex] = np.ones(tmp_IWV.shape[1])*year
    bg_month[sindex:eindex] = np.ones(tmp_IWV.shape[1])*month
    
    
for year in np.arange(2004,2015):
    print('working on ', year)
    for month in np.arange(1,13):
        tmp_SSTmean, tmp_IVT, tmp_IWV, tmp_IVTf, tmp_IWVf = compute_stats(year, month)
    
        sindex = eindex
        eindex = eindex + calendar.monthrange(year, month)[1]*4
        stats_SSTmean[:, sindex:eindex] = tmp_SSTmean
        stats_IVT[:, sindex:eindex] = tmp_IVT
        stats_IWV[:, sindex:eindex] = tmp_IWV
        stats_IVTf[:, sindex:eindex] = tmp_IVTf
        stats_IWVf[:, sindex:eindex] = tmp_IWVf
        bg_year[sindex:eindex] = np.ones(tmp_IWV.shape[1])*year
        bg_month[sindex:eindex] = np.ones(tmp_IWV.shape[1])*month
        
year = 2015
print('working on ', year)
for month in np.arange(1,10):
    tmp_SSTmean, tmp_IVT, tmp_IWV, tmp_IVTf, tmp_IWVf = compute_stats(year, month)
    
    sindex = eindex
    eindex = eindex + calendar.monthrange(year, month)[1]*4
    stats_SSTmean[:, sindex:eindex] = tmp_SSTmean
    stats_IVT[:, sindex:eindex] = tmp_IVT
    stats_IWV[:, sindex:eindex] = tmp_IWV
    stats_IVTf[:, sindex:eindex] = tmp_IVTf
    stats_IWVf[:, sindex:eindex] = tmp_IWVf
    bg_year[sindex:eindex] = np.ones(tmp_IWV.shape[1])*year
    bg_month[sindex:eindex] = np.ones(tmp_IWV.shape[1])*month


tmpfile =  rootdir + 'intermediate_data/full_ocean_stats_separate.mat'
sio.savemat(tmpfile, {'stats_SSTmean':stats_SSTmean, 
                      'stats_IVT':stats_IVT, 'stats_IWV':stats_IWV,
                      'stats_IVTf':stats_IVTf, 'stats_IWVf':stats_IWVf,
                      'bg_year':bg_year, 'bg_month':bg_month})

