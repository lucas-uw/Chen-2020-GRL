#!/share/apps/python/anaconda3.6/bin/python
######!/people/chen423/sw/anaconda3/bin/python

import numpy as np
import xarray as xr
import scipy.io as sio
import pandas as pd
import calendar
import sys

rootdir = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/'

def crt_filenames(model, year, month):
    
    WRFdir = rootdir + '%s/by_month_SERDP6km/' % model
    WRFfile = WRFdir + 'WRF_IWV_uvIVT.6hr.%d.%d.nc' % (year, month)
    
    return WRFfile

def get_AR_intensity_data(model, year, month):
    
    WRFfile = crt_filenames(model, year, month)
    
    WRF_IVT = xr.open_dataset(WRFfile).uvIVT.values
    WRF_IWV = xr.open_dataset(WRFfile).IWV.values
    
    return WRF_IVT, WRF_IWV

def compute_6hrly_AR_SST(in_ARtag, in_SST):
    ocean_AR_union = (ocean_mask==1)*(in_ARtag==1)
    out_SSTmean = in_SST[ocean_AR_union==1].mean()
    
    return out_SSTmean

def compute_6hrly_AR_intensity(in_ARtag, in_intensity):
    land_AR_union = (ocean_mask==0)*(in_ARtag==1)
    out_intensity = in_intensity[land_AR_union==1].mean()
    
    return out_intensity

def compute_6hrly_AR_intrusion(in_ARtag):
    out_dist_max = (dist_to_coast[in_ARtag==1]).max()
    
    return out_dist_max

def compute_AR_stats_separateAR(year, month, ARtag='p85', flag_area=-9999, flag_USstate=-9999, flag_post_adj=-9999):
    if flag_post_adj==1:
        file_ARHIST = rootdir + 'HIST/AR_tagged/Gershunov/SERDP6km_adj/WRF_ARtag_adj.HIST.Gershunov.%d.%d.AR%s.nc' % (year, month, ARtag)
        file_ARfSST = rootdir + 'fSST/AR_tagged/Gershunov/SERDP6km_adj/WRF_ARtag_adj.fSST.Gershunov.%d.%d.AR%s.nc' % (year, month, ARtag)
    elif flag_post_adj==0:
        file_ARHIST = rootdir + 'HIST/AR_tagged/Gershunov/SERDP6km/WRF_ARtag.HIST.Gershunov.%d.%d.AR%s.nc' % (year, month, ARtag)
        file_ARfSST = rootdir + 'fSST/AR_tagged/Gershunov/SERDP6km/WRF_ARtag.fSST.Gershunov.%d.%d.AR%s.nc' % (year, month, ARtag)
    
    
    file_SSTHIST = rootdir + 'HIST/SST/NARR_TS.SERDP6km.6hourly.%d.%d.nc' % (year, month)
    
    file_SSTfix = rootdir + 'HIST/SST/NARR_TS.SERDP6km.2000.10.01.00.nc'
    
    
    ARtag_HIST = xr.open_dataset(file_ARHIST).AR_tag.values
    SST_HIST = xr.open_dataset(file_SSTHIST).var11.values
    IVT_HIST, IWV_HIST = get_AR_intensity_data('HIST', year, month)

    ARtag_fSST = xr.open_dataset(file_ARfSST).AR_tag.values
    SST_fSST = xr.open_dataset(file_SSTfix).var11.values[0]
    IVT_fSST, IWV_fSST = get_AR_intensity_data('fSST', year, month)

    
    # compute various stats
    nt = ARtag_HIST.shape[0]
    stat_AR_SSTmean = np.zeros((2,nt))-9999
    stat_AR_dist = np.zeros((2,nt))-9999
    stat_AR_landarea = np.zeros((2,nt))-9999
    stat_AR_IVT = np.zeros((2,nt))-9999
    stat_AR_IWV = np.zeros((2,nt))-9999
    valid_index = np.zeros((2,nt))
    common_AR = np.zeros(nt)

    for t in np.arange(nt):
        if flag_USstate==1:
            sig1 = ((ARtag_HIST[t]==1)*(ocean_mask==0)*(USstate==0)).sum() # land
            sig3 = ((ARtag_fSST[t]==1)*(ocean_mask==0)*(USstate==0)).sum() # land
        elif flag_USstate==0:
            sig1 = ((ARtag_HIST[t]==1)*(ocean_mask==0)).sum() # land
            sig3 = ((ARtag_fSST[t]==1)*(ocean_mask==0)).sum() # land
        sig2 = ((ARtag_HIST[t]==1)*(ocean_mask==1)).sum()  # ocean
        sig4 = ((ARtag_fSST[t]==1)*(ocean_mask==1)).sum() # ocean
        sig5 = (ARtag_HIST[t]*ARtag_fSST[t]*(ocean_mask==1)).sum()
        #print(t, sig1, sig2, sig3, sig4)
        if sig1>flag_area and sig2>flag_area:
            valid_index[0,t] = 1
            stat_AR_SSTmean[0,t] = compute_6hrly_AR_SST(ARtag_HIST[t], SST_HIST[t])
            stat_AR_dist[0,t] = compute_6hrly_AR_intrusion(ARtag_HIST[t])
            stat_AR_landarea[0,t] = sig1
            stat_AR_IVT[0,t] = compute_6hrly_AR_intensity(ARtag_HIST[t], IVT_HIST[t])
            stat_AR_IWV[0,t] = compute_6hrly_AR_intensity(ARtag_HIST[t], IWV_HIST[t])
            
        if sig3>flag_area and sig4>flag_area:
            valid_index[1,t] = 1
            stat_AR_SSTmean[1,t] = compute_6hrly_AR_SST(ARtag_fSST[t], SST_fSST)
            stat_AR_dist[1,t] = compute_6hrly_AR_intrusion(ARtag_fSST[t])
            stat_AR_landarea[1,t] = sig3
            stat_AR_IVT[1,t] = compute_6hrly_AR_intensity(ARtag_fSST[t], IVT_fSST[t])
            stat_AR_IWV[1,t] = compute_6hrly_AR_intensity(ARtag_fSST[t], IWV_fSST[t])
            
        if sig1>flag_area and sig2>flag_area and sig3>flag_area and sig4>flag_area and sig5>commonAR_thre:
            common_AR[t] = 1
            
    return stat_AR_SSTmean, stat_AR_dist, stat_AR_landarea, stat_AR_IVT, stat_AR_IWV, valid_index, common_AR



ARtag = sys.argv[1]
flag_area = int(sys.argv[2])   # minimum size of patches (over land and over ocean, both)
flag_USstate = int(sys.argv[3])   # whether to use US west coast 5 states along with land mask. 1 is to use, 0 is to skip
flag_post_adj = int(sys.argv[4])   # WRF further adjusted, or not (i.e., directly from modified NARR). 1 is further adjusted, 0 for raw
commonAR_thre = int(sys.argv[5])

version_tag = 'AR%s_s%d_state%d_post%d_c%d' % (ARtag, flag_area, flag_USstate, flag_post_adj, commonAR_thre)
print(version_tag)

reffile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/ref/SERDP6km.dist_to_coastal.nc'
dist_to_coast = xr.open_dataset(reffile).dist_to_coast.values
dist_to_coast[dist_to_coast==9999] = 0
ocean_mask = np.zeros((450,450))
ocean_mask[dist_to_coast==0] = 1

reffile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/ref/US_state.nc'
USstate = 1-xr.open_dataset(reffile).state_mask.values[0:5].sum(axis=0)


## part 1
stats_AR_SSTmean = np.zeros((2,17532))-9999
stats_AR_dist = np.zeros((2,17532))-9999
stats_AR_landarea = np.zeros((2,17532))-9999
stats_AR_IVT = np.zeros((2,17532))-9999 # over land
stats_AR_IWV = np.zeros((2,17532))-9999 # over land
bg_year = np.zeros(17532)-9999
bg_month = np.zeros(17532)-9999
ARday_index = np.zeros((2,17532))-9999
commonAR = np.zeros(17532)-9999

sindex = -31*4
eindex = 0

year = 2003
print('working on ', year)
for month in np.arange(10,13):
    tmp_SSTmean, tmp_dist, tmp_landarea, tmp_IVT, tmp_IWV, tmp_vindex, tmp_c = compute_AR_stats_separateAR(year, month,
                                                                                                           ARtag=ARtag,
                                                                                                           flag_area=flag_area,
                                                                                                           flag_USstate=flag_USstate,
                                                                                                           flag_post_adj=flag_post_adj)
    
    sindex = eindex
    eindex = eindex + calendar.monthrange(year, month)[1]*4
    stats_AR_SSTmean[:, sindex:eindex] = tmp_SSTmean
    stats_AR_dist[:, sindex:eindex] = tmp_dist
    stats_AR_landarea[:, sindex:eindex] = tmp_landarea
    stats_AR_IVT[:, sindex:eindex] = tmp_IVT
    stats_AR_IWV[:, sindex:eindex] = tmp_IWV
    ARday_index[:, sindex:eindex] = tmp_vindex
    bg_year[sindex:eindex] = np.ones(tmp_vindex.shape[1])*year
    bg_month[sindex:eindex] = np.ones(tmp_vindex.shape[1])*month
    commonAR[sindex:eindex] = tmp_c
    
    
for year in np.arange(2004,2015):
    print('working on ', year)
    for month in np.arange(1,13):
        tmp_SSTmean, tmp_dist, tmp_landarea, tmp_IVT, tmp_IWV, tmp_vindex, tmp_c = compute_AR_stats_separateAR(year, month,
                                                                                                               ARtag=ARtag,
                                                                                                               flag_area=flag_area,
                                                                                                               flag_USstate=flag_USstate,
                                                                                                               flag_post_adj=flag_post_adj)
    
        sindex = eindex
        eindex = eindex + calendar.monthrange(year, month)[1]*4
        stats_AR_SSTmean[:, sindex:eindex] = tmp_SSTmean
        stats_AR_dist[:, sindex:eindex] = tmp_dist
        stats_AR_landarea[:, sindex:eindex] = tmp_landarea
        stats_AR_IVT[:, sindex:eindex] = tmp_IVT
        stats_AR_IWV[:, sindex:eindex] = tmp_IWV
        ARday_index[:, sindex:eindex] = tmp_vindex
        bg_year[sindex:eindex] = np.ones(tmp_vindex.shape[1])*year
        bg_month[sindex:eindex] = np.ones(tmp_vindex.shape[1])*month
        commonAR[sindex:eindex] = tmp_c
        
year = 2015
print('working on ', year)
for month in np.arange(1,10):
    tmp_SSTmean, tmp_dist, tmp_landarea, tmp_IVT, tmp_IWV, tmp_vindex, tmp_c = compute_AR_stats_separateAR(year, month,
                                                                                                           ARtag=ARtag,
                                                                                                           flag_area=flag_area,
                                                                                                           flag_USstate=flag_USstate,
                                                                                                           flag_post_adj=flag_post_adj)
    
    sindex = eindex
    eindex = eindex + calendar.monthrange(year, month)[1]*4
    stats_AR_SSTmean[:, sindex:eindex] = tmp_SSTmean
    stats_AR_dist[:, sindex:eindex] = tmp_dist
    stats_AR_landarea[:, sindex:eindex] = tmp_landarea
    stats_AR_IVT[:, sindex:eindex] = tmp_IVT
    stats_AR_IWV[:, sindex:eindex] = tmp_IWV
    ARday_index[:, sindex:eindex] = tmp_vindex
    bg_year[sindex:eindex] = np.ones(tmp_vindex.shape[1])*year
    bg_month[sindex:eindex] = np.ones(tmp_vindex.shape[1])*month
    commonAR[sindex:eindex] = tmp_c


tmpfile = rootdir + 'intermediate_data/AR_stats_separate.%s.mat' % (version_tag)
sio.savemat(tmpfile, {'stats_AR_SSTmean':stats_AR_SSTmean, 'stats_AR_dist':stats_AR_dist,
                      'stats_AR_landarea':stats_AR_landarea, 'stats_AR_IVT':stats_AR_IVT,
                      'stats_AR_IWV':stats_AR_IWV, 'ARday_index':ARday_index,
                      'bg_year':bg_year, 'bg_month':bg_month, 'commonAR':commonAR})




# part 2
def get_AR_maps(model, year, month):
    
    ARfile = rootdir + '%s/AR_tagged/Gershunov/SERDP6km_adj/WRF_ARtag_adj.%s.Gershunov.%d.%d.AR%s.nc' % (model, model, year, month, ARtag)
    AR_maps = xr.open_dataset(ARfile).AR_tag.values
    
    return AR_maps

ts_full = pd.period_range(start='2003-10-01-00', end='2015-09-30-18', freq='6H')

def compute_monthly_stats(year, month):
    
    # need to use ts_full and ARday_index
    
    ARtag_HIST = get_AR_maps('HIST', year, month)
    ARtag_fSST = get_AR_maps('fSST', year, month)
    
    ARindex_clip = ARday_index[:, ((ts_full.year==year)*(ts_full.month==month))==1]
    nt = ARindex_clip.shape[1]
    
    sum_HIST1 = np.zeros((450,450)) # common
    sum_HIST2 = np.zeros((450,450)) # only HIST
    sum_fSST1 = np.zeros((450,450)) # common
    sum_fSST2 = np.zeros((450,450)) # only fSST
    
    for t in np.arange(nt):
        if ARindex_clip[0,t]==1 and ARindex_clip[1,t]==1:  # common days
            sum_HIST1 = sum_HIST1 + ARtag_HIST[t]
            sum_fSST1 = sum_fSST1 + ARtag_fSST[t]
        if ARindex_clip[0,t]==1 and ARindex_clip[1,t]==0: # only HIST
            sum_HIST2 = sum_HIST2 + ARtag_HIST[t]
        if ARindex_clip[0,t]==0 and ARindex_clip[1,t]==1: # only fSST
            sum_fSST2 = sum_fSST2 + ARtag_fSST[t]
        
    return sum_HIST1, sum_HIST2, sum_fSST1, sum_fSST2, nt
    
frac_HIST1 = np.zeros((144,450,450))
frac_HIST2 = np.zeros((144,450,450))
frac_fSST1 = np.zeros((144,450,450))
frac_fSST2 = np.zeros((144,450,450))

count_HIST1 = np.zeros((144,450,450))
count_HIST2 = np.zeros((144,450,450))
count_fSST1 = np.zeros((144,450,450))
count_fSST2 = np.zeros((144,450,450))

count = 0
year = 2003
print(year)
for month in np.arange(10,13):
    subdata1, subdata2, subdata3, subdata4, nt = compute_monthly_stats(year, month)

    frac_HIST1[count] = subdata1/nt
    frac_HIST2[count] = subdata2/nt
    frac_fSST1[count] = subdata3/nt
    frac_fSST2[count] = subdata4/nt
    count_HIST1[count] = subdata1
    count_HIST2[count] = subdata2
    count_fSST1[count] = subdata3
    count_fSST2[count] = subdata4
    count = count + 1
    
for year in np.arange(2004,2015):
    print(year)
    for month in np.arange(1,13):
        subdata1, subdata2, subdata3, subdata4, nt = compute_monthly_stats(year, month)

        frac_HIST1[count] = subdata1/nt
        frac_HIST2[count] = subdata2/nt
        frac_fSST1[count] = subdata3/nt
        frac_fSST2[count] = subdata4/nt
        count_HIST1[count] = subdata1
        count_HIST2[count] = subdata2
        count_fSST1[count] = subdata3
        count_fSST2[count] = subdata4
        count = count + 1
        
year = 2015
print(year)
for month in np.arange(1,10):
    subdata1, subdata2, subdata3, subdata4, nt = compute_monthly_stats(year, month)

    frac_HIST1[count] = subdata1/nt
    frac_HIST2[count] = subdata2/nt
    frac_fSST1[count] = subdata3/nt
    frac_fSST2[count] = subdata4/nt
    count_HIST1[count] = subdata1
    count_HIST2[count] = subdata2
    count_fSST1[count] = subdata3
    count_fSST2[count] = subdata4
    count = count + 1
    
    
tmpfile = rootdir + 'intermediate_data/ARstats.monthly_frac.%s.mat' % version_tag
sio.savemat(tmpfile, {'frac_HIST1':frac_HIST1, 'frac_HIST2':frac_HIST2,
                      'frac_fSST1':frac_fSST1, 'frac_fSST2':frac_fSST2})

    
tmpfile = rootdir + 'intermediate_data/ARstats.monthly_count.%s.mat' % version_tag
sio.savemat(tmpfile, {'count_HIST1':count_HIST1, 'count_HIST2':count_HIST2,
                      'count_fSST1':count_fSST1, 'count_fSST2':count_fSST2})
