#!/people/chen423/sw/anaconda3/bin/python

import numpy as np
import xarray as xr
import netCDF4 as nc
from skimage.draw import polygon
from skimage import measure
import cv2 as cv
from calendar import monthrange
import sys

model = sys.argv[1]

NARR_mask = np.zeros((279,351))
NARR_coast = np.zeros((279,351))

reffile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/ref/NARR_landmask_coastline.nc'
NARR_mask[1:278,1:350] = xr.open_dataset(reffile).landmask.values
NARR_coast[1:278,1:350] = xr.open_dataset(reffile).coastline.values
reffile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/ref/NARR_full_latlon.nc'
NARR_lats = xr.open_dataset(reffile).lat.values
NARR_lons = xr.open_dataset(reffile).lon.values

NARR_gridsize = 32 # km

# misc. functions around IO
def crt_filenames(model, year, month):
    NARRfile = '/pic/projects/next_gen_idf/chen423/data/temp/narr_AR/narr_AR.%d.%d.nc'%(year, month)

    WRFdir = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/NARR_by_month/' % model        
    WRFfile_u = WRFdir + 'WRF_NARR_%s.NARRgrid.6hr.uIVT.%d.%d.nc' % (model, year, month)
    WRFfile_v = WRFdir + 'WRF_NARR_%s.NARRgrid.6hr.vIVT.%d.%d.nc' % (model, year, month)
    WRFfile_IWV = WRFdir + 'WRF_NARR_%s.NARRgrid.6hr.IWV.%d.%d.nc' % (model, year, month)
    
    return NARRfile, WRFfile_u, WRFfile_v, WRFfile_IWV

def get_NARRraw_data(infile, var):
    # 0:248:2 is used to extract 6hourly data
    tmpgroup = xr.open_dataset(infile)
    if var=='uvIVT':
        tmpdata1 = tmpgroup['uIVT'].values[0:248:2]
        tmpdata2 = tmpgroup['vIVT'].values[0:248:2]
        outdata = np.sqrt(tmpdata1*tmpdata1+tmpdata2*tmpdata2)
    else:
        outdata = tmpgroup[var].values[0:248:2]
    tmpgroup.close()
    return outdata

def gen_mod_NARRdata(model, year, month):
    
    NARRfile, WRFfile_u, WRFfile_v, WRFfile_IWV = crt_filenames(model, year, month)
    
    data_IVT = get_NARRraw_data(NARRfile, 'uvIVT')
    WRF_IVT = np.sqrt(xr.open_dataset(WRFfile_u).uIVT.values**2 + xr.open_dataset(WRFfile_v).vIVT.values**2)
    data_IVT[(np.isinf(data_IVT))|(np.isnan(data_IVT))] = 0
    
    data_IWV = get_NARRraw_data(NARRfile, 'IWV')
    WRF_IWV = xr.open_dataset(WRFfile_IWV).IWV.values
    data_IWV[(np.isinf(data_IWV))|(np.isnan(data_IWV))] = 0
        
    data_IVT[np.where(np.isfinite(WRF_IVT)==True)] = WRF_IVT[np.where(np.isfinite(WRF_IVT)==True)]
    data_IWV[np.where(np.isfinite(WRF_IWV)==True)] = WRF_IWV[np.where(np.isfinite(WRF_IWV)==True)]
    
    return data_IVT, data_IWV

def concatenate_file(model, year, month):
    yearp = year
    monthp = month-1
    if monthp==0:
        yearp = yearp - 1
        monthp = 12
    yearn = year
    monthn = month + 1
    if monthn==13:
        yearn = yearn + 1
        monthn = 1
    
    if year==2003 and month==10:
        #print('here is the case')
        IVTdata, IWVdata = gen_mod_NARRdata(model, year, month)
        IVTdatan, IWVdatan = gen_mod_NARRdata(model, yearn, monthn)
        out_IVT = np.concatenate((IVTdata, IVTdatan), axis=0)
        out_IWV = np.concatenate((IWVdata, IWVdatan), axis=0)
        sindex = 0
        eindex = monthrange(year, month)[1]*4
    elif year==2015 and month==12:
        IVTdata, IWVdata = gen_mod_NARRdata(model, year, month)
        IVTdatap, IWVdatap = gen_mod_NARRdata(model, yearp, monthp)
        out_IVT = np.concatenate((IVTdatap, IVTdata), axis=0)
        out_IWV = np.concatenate((IWVdatap, IWVdata), axis=0)
        sindex = monthrange(yearp, monthp)[1]*4
        eindex = monthrange(yearp, monthp)[1]*4 + monthrange(year, month)[1]*4
    else:
        IVTdata, IWVdata = gen_mod_NARRdata(model, year, month)
        IVTdatap, IWVdatap = gen_mod_NARRdata(model, yearp, monthp)
        IVTdatan, IWVdatan = gen_mod_NARRdata(model, yearn, monthn)
        out_IVT = np.concatenate((IVTdatap, IVTdata, IVTdatan), axis=0)
        out_IWV = np.concatenate((IWVdatap, IWVdata, IWVdatan), axis=0)
        sindex = monthrange(yearp, monthp)[1]*4
        eindex = monthrange(yearp, monthp)[1]*4 + monthrange(year, month)[1]*4
        
    return out_IVT, out_IWV, sindex, eindex


def save_to_nc(outfile, stime, indata):
    outgroup = nc.Dataset(outfile, 'w', format='NETCDF4')

    # dimension
    tdim = outgroup.createDimension('time', None)
    xdim = outgroup.createDimension('x', 349)
    ydim = outgroup.createDimension('y', 277)

    # vars
    tvar = outgroup.createVariable('time', 'i4', ('time'))
    tvar.units = 'Hours since %s-01 00:00' % stime
    tvar[:] = (np.arange(indata.shape[0])*6).astype('int32')
    xvar = outgroup.createVariable('x', 'i4', ('x'))
    xvar[:] = np.arange(349)
    yvar = outgroup.createVariable('y', 'i4', ('y'))
    yvar[:] = np.arange(277)

    var_lat = outgroup.createVariable('lat', 'f4', ('y', 'x'))
    var_lat.long_name = 'Latitude'
    var_lat.standard_name = 'latitude'
    var_lat.units = 'degrees_north'
    var_lon = outgroup.createVariable('lon', 'f4', ('y', 'x'))
    var_lon.long_name = 'Longitude'
    var_lon.standard_name = 'longitude'
    var_lon.units = 'degrees_east'
    var_tag = outgroup.createVariable('AR_tag', 'f4', ('time','y', 'x'))
    var_tag.long_name = 'AR tag based on Gershunov approach and 6-hourly data'
    var_tag.coordinates = 'lon lat'
    
    var_lat[:] = NARR_lats
    var_lon[:] = NARR_lons
    var_tag[:] = indata

    outgroup.close()


# my AR identification toolkit
def delinerate_potential_AR_cells(in_IVT, in_IWV):
    # threshold taken from Gershunov algorithm
    out_tag = np.logical_and(in_IVT>250, in_IWV>15)
    return out_tag

def crt_high_intense_contours(AR_tag):

    # the boundary should be 1, but when 0.99999 here is accompanied with round_ function,
    #  the correct grids are still identified
    contours = measure.find_contours(AR_tag, 0.99999)
    
    return contours

def find_contour_minEncCir(in_contour):
    # give polygon, return the perimeter
    
    cnt = np.array(np.round_(in_contour), dtype=int)
    (x,y),radius = cv.minEnclosingCircle(cnt) 
    
    return radius*2

def analyze_individual_contour(contour_single, IVT_map):
    contour_single_round = np.round_(contour_single)
    conclusion = False
    center_lat = 0
    length = 0
    out_rr = -9999
    out_cc = -9999
    if contour_single_round.shape[0]>20:
        tmp_patch_matrix = np.zeros((279,351))
        rr, cc = polygon(contour_single_round[:,0], contour_single_round[:,1], (279,351))
        #print(rr, cc)
        tmp_patch_matrix[rr,cc] = 1
        
        # 1. check the intersection with coastal line
        signal = (tmp_patch_matrix*NARR_coast).sum()
        if signal>1:
            tmp_intensity_dist = IVT_map*tmp_patch_matrix*NARR_coast
            #print('debug1: ', np.isinf(IVT_map).sum(), IVT_map.mean())
            #print('debug2: ', np.isinf(tmp_patch_matrix).sum(), tmp_patch_matrix.mean())
            #print('debug3: ', np.isinf(NARR_coast).sum(), NARR_coast.mean())
            #print('debug: max dist: ', tmp_intensity_dist.max())
            max_intensity = tmp_intensity_dist.max()
            location = np.where(tmp_intensity_dist==max_intensity)
            #print('max_intensity: ', max_intensity)
            center_lat = NARR_lats[location]
            center_lon = NARR_lons[location]
            #print(max_intensity, location)
            
            # 2. check the length of patch, it needs to be longer than 1500km
            MECperi = find_contour_minEncCir(contour_single_round)
            if MECperi*NARR_gridsize>=1500:
                #print(MECperi*NARR_gridsize, max_intensity, location)
                conclusion = True
                out_rr = rr
                out_cc = cc
                length = MECperi
            
    return conclusion, center_lat, length, out_rr, out_cc

def idenfity_AR(model, year, month):
    
    # get data for previous, current, and next month
    print(year, month, 'reading input ...')
    IVT_raw, IWV_raw, sindex, eindex = concatenate_file(model,year, month)

    nt = IVT_raw.shape[0]
    nx = IVT_raw.shape[1]
    ny = IVT_raw.shape[2]

    IVT_full = np.zeros((nt, nx+2, ny+2))
    IVT_full[:,1:(nx+1),1:(ny+1)] = IVT_raw

    IWV_full = np.zeros((nt, nx+2, ny+2))
    IWV_full[:,1:(nx+1),1:(ny+1)] = IWV_raw
        
    #print(IVT_full.shape, IWV_full.shape, sindex, eindex)
    
    # compute high intense region
    print('  computing high intense region ...')
    high_intens_tag = delinerate_potential_AR_cells(IVT_full, IWV_full)
    high_intens_tag[:,:,220:349] = np.zeros((nt,IVT_full.shape[1],129))
    
    # orgrnize into individual maps (one cluster in one map)
    print('  organizing into individual maps ...')
    AR_location_time_matrix = np.zeros((nt,10))
    AR_size_matrix = np.zeros((nt,10))
    AR_tag_matrix_stage1 = np.zeros((nt,10,IVT_full.shape[1], IVT_full.shape[2]))

    for i in np.arange(nt):
        IVT_step = IVT_full[i]
        raw_contours = crt_high_intense_contours(high_intens_tag[i,:,:])
        nARs = len(raw_contours)
        #print(nARs)
        count = 0
        for n in np.arange(nARs):
            conclusion, location, AR_length, rr, cc = analyze_individual_contour(raw_contours[n], IVT_step)
            if conclusion==True:
                #print(i, n, location)
                AR_location_time_matrix[i,count] = location
                AR_size_matrix[i,count] = AR_length
                AR_tag_matrix_stage1[i,count,:,:][rr,cc] = 1
                count = count + 1
                
    # considering 18-hour requirement and 2.5-degree spatial shift requirement
    print('  cleaning clusters ...')
    AR_tag_matrix_stage2 = np.zeros((nt, 10, IVT_full.shape[1], IVT_full.shape[2]))

    for i in np.arange(nt):  # nt
        #print(i)
        for j in np.arange(10):
            location_current = AR_location_time_matrix[i,j]
            if location_current>0:
                sig_backward = np.zeros(10) # meaning no succedding occurrence of same event
                sig_forward = np.zeros(10)
                loc_backward = np.zeros(10)
                loc_forward = np.zeros(10)
                # forward
                loc_forward[0] = location_current
                sig_forward[0] = 1
                for t in np.arange(1,6): # loop over backward/forward time step
                    if i+t<=(nt-1):
                        for k in np.arange(10): # loop over each cluster at each time step
                            if AR_location_time_matrix[i+t,k]>0:
                                #print(i+t,k,AR_location_time_matrix[i+t,k], loc_forward[t-1])
                                if abs(AR_location_time_matrix[i+t,k]-loc_forward[t-1])<=5: # 5 in Gershonove which is based on 6-hr data
                                    sig_forward[t] = 1
                                    loc_forward[t] = AR_location_time_matrix[i+t, k]
                # backward
                loc_backward[0] = location_current
                sig_backward[0] = 1
                for t in np.arange(1,6):
                    if i-t>=0:
                        for k in np.arange(10):
                            if AR_location_time_matrix[i-t,k]>0:
                                if abs(AR_location_time_matrix[i-t,k]-loc_backward[t-1])<=5: # see note above
                                    sig_backward[t] = 1
                                    loc_backward[t] = AR_location_time_matrix[i-t,k]
                #print(i,j, sig_backward.sum()+sig_forward.sum()-1)              
                if sig_backward.sum()+sig_forward.sum()>=4: # do some mannual derivation, basically 6 snaps, but 0th is the same
                    AR_tag_matrix_stage2[i,j] = AR_tag_matrix_stage1[i,j]
                    
    # 3. compositve maps for each time step
    AR_tag_matrix_stage3 = AR_tag_matrix_stage2.sum(axis=1)
    
    # write results
    print('  writing ...')
    outfile = '/pic/projects/hyperion/chen423/data/papers/AR-SST/data/%s/AR_tagged/Gershunov/NARR/NARR_ARtag.%s.Gershunov.%d.%d.nc' % (model, model, year, month)
           
    time_string = '%d-%02d' % (year, month)
    save_to_nc(outfile, time_string, AR_tag_matrix_stage3[sindex:eindex,1:278,1:350])

for month in np.arange(10,13):
	print(model, 2003, month)
	idenfity_AR(model, 2003, month)


for year in np.arange(2004, 2016):
	for month in np.arange(1,13):
		print(model, year, month)
		idenfity_AR(model, year, month)
