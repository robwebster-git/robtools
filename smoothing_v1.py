#!/usr/bin/env python3

import rasterio as rio
from rasterio.mask import mask
from rasterio import features
import numpy as np
from scipy.ndimage import median_filter
from skimage.morphology import area_opening, area_closing, dilation, erosion, disk, square
from tqdm import tqdm
import geopandas as gpd
from shapely import geometry
import json
import geojson
from rasterstats import zonal_stats
from glob import glob
import os
import re
import click
from subprocess import call

dbh_sw_bins = [4, 12, 22, 37]
dbh_hw_bins = [4, 12, 29, 37]
cc_bins = [20, 35, 60, 80]

def rasterize(gdf, meta, column='STRAT_ZON', fill=-9999):
    """
    Convert geodataframe to raster based on a column's values
    """
    gdf = gdf.to_crs(meta['crs'])
    geoms = ((g,v) for g, v in zip(gdf.geometry, gdf[column]))
    rasterized = features.rasterize(geoms, out_shape=(meta['height'], meta['width']), transform=meta['transform'], all_touched=False, fill=fill, dtype=np.int16)
    return np.ma.masked_array(rasterized, rasterized==-9999, fill_value=-9999)

def polygonize(in_data, ref_tif, colname='STRATA'):
    """
    Not used. rio shapes is called instead. Somehow it is faster.
    """
    gdf = gpd.GeoDataFrame.from_dict(list(features.shapes(in_data, transform=ref_tif.transform)))
    gdf['geometry'] = gdf[0].apply(lambda x: geometry.shape(geojson.loads(json.dumps(x))))
    gdf.crs = ref_tif.crs
    gdf = gdf.drop(0, axis=1)
    gdf = gdf.rename({1: colname}, axis=1)
    return gdf

def smooth_tif(data, mask, area = 41):
    """
    Run sieve on a raster in a step-wise manner. I did not like what sieve does to smaller shapes when run in one-shot
    so I made sure that it goes up to the value in steps. It's slower but it works better.
    """
    smooth = None
    #steps = np.linspace(0.01, 1, 6)
    steps = [1]
    for a in steps:
        if smooth is None:
            smooth = features.sieve(data.astype(np.int16), size=int(area*a), mask=mask)
        smooth = features.sieve(smooth, size=int(area*a), mask=mask)
#    strat_smooth[strat_smooth == 0] = -9999
    return np.ma.masked_array(smooth, smooth == -9999, fill_value=-9999)

def prep_multi_cc(multi_cc_file, outdir):
    
    """Load CC data, bin and iteratively filter/smooth"""
    #multi_cc_file = '../inputs/inputs_crop/nd_crop_cmp_merged_PredMultiLayer_properties_clipped.tif' # example!!!!
    
    with rio.open(multi_cc_file) as src:
        profile = src.profile.copy()    
        mult_cc_data = src.read(masked = True)     
    masked_multi_cc_med = np.ma.masked_array(mult_cc_data, mult_cc_data == -9999, fill_data=-9999)

    multi_cc_bins = [200]
    print(f'--CC BINS = {multi_cc_bins}--')
    multi_cc_bin = np.digitize(masked_multi_cc_med, multi_cc_bins)
    multi_cc_bin[mult_cc_data < 1] = -9999
    ma_multi_cc_bin = np.ma.masked_array(multi_cc_bin, mask = mult_cc_data.mask, fill_data=-9999)
    save_raster(ma_multi_cc_bin, profile, f'{outdir}/multi_cc_bin.tif')


    multi_cc_med = median_filter(ma_multi_cc_bin, 3)
    multi_cc_med = np.ma.masked_array(multi_cc_med, multi_cc_med == -9999, fill_data=-9999)
    
    save_raster(multi_cc_med, profile, f'{outdir}/multi_cc_bin_med.tif')
    
    # how many bins actually are there? May not need any of this looping. 
    for multi_cc_group in [0, 1]:

        tmp = np.ma.masked_array(multi_cc_med, mask=multi_cc_med.mask | (~np.isin(multi_cc_med, multi_cc_group)))
        multi_cc_smooth = smooth_tif(tmp, mask=~tmp.mask[0], area=ac_to_px(5))

    profile['dtype'] = np.int16
    profile['nodata'] = -9999
    with rio.open(f'{outdir}/multi_cc_bin_med_smooth.tif' , 'w', **profile) as dst:
        dst.write(multi_cc_smooth.astype(np.int16), 1)
    
    return multi_cc_smooth


def simple_smooth(data, harvest, area):
    
    smooth = np.ma.masked_array(data, mask=(data == -9999), dtype=np.int32)
    
    tmp = np.ma.masked_array(smooth, mask=((data.mask) & (harvest != 1))) # Sieve inside harvest areas
    smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
    tmp = np.ma.masked_array(smooth, mask=((data.mask) & (harvest == 1))) # Sieve outside harvest areas
    smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)

    return smooth


def dbh_isolation_sieve(data, supertype, harvest, area=404):
    """
    Perform a sieve masking out areas of DBH bins. That is, low-valued classes get treated separate
    from the high-valued classes, both within and outwith harvest areas.
    """
    smooth = np.ma.masked_array(data, mask=(data == -9999), dtype=np.int32)
    
    for dbh_grp in [[0, 1], [3,4], [2, 3, 4]]:
        tmp = np.ma.masked_array(smooth, mask=(((supertype.mask) | (data.mask)) | (~np.isin(smooth, dbh_grp))) & (harvest == 1)) # Sieve inside harvest areas
        smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
        tmp = np.ma.masked_array(smooth, mask=(((supertype.mask) | (data.mask)) | (~np.isin(smooth, dbh_grp))) & (harvest != 1)) # Sieve outside harvest areas
        smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)

#    smooth[smooth == 0] = -9999
    
    return np.ma.masked_array(smooth, smooth == -9999, fill_value=-9999)


def st_isolation_sieve(data, dbh, harvest, area=404):
    """
    Perform a sieve masking out areas of DBH bins and supertype classes. 
    That is, low-valued DBH classes get treated separate from the high-valued classes, both within and outwith harvest areas.
    Simultaneously we create groupings in the mixed, hardwood, and softwood areas of supertype.
    """
    smooth = np.ma.masked_array(data, mask=(data == -9999), dtype=np.int32)
    
    for dbh_grp in tqdm([[0, 1], [2, 3, 4]]):
        for st_grp in [[2, 3], [1, 2], [3, 4], list(range(1, 5))]: 
            tmp = np.ma.masked_array(data, mask=(((dbh.mask) | (data.mask)) | (~np.isin(dbh, dbh_grp)) | (~np.isin(data, st_grp)) & (harvest != 1)))
            smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
            tmp = np.ma.masked_array(data, mask=(((dbh.mask) | (data.mask)) | (~np.isin(dbh, dbh_grp)) | (~np.isin(data, st_grp)) & (harvest == 1)))
            smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)

    return np.ma.masked_array(smooth, (smooth == -9999) | (smooth == 0), fill_value=-9999)


def cc_isolation_sieve(data, dbh, supertype, harvest, area=404):
    """
    Perform a sieve masking out areas of DBH bins and supertype classes. 
    That is, low-valued DBH classes get treated separate from the high-valued classes, both within and outwith harvest areas.
    Simultaneously we create groupings in the mixed, hardwood, and softwood areas of supertype.
    """
    smooth = np.ma.masked_array(data, mask=(data == -9999), dtype=np.int32)
    
    for dbh_grp in tqdm([[0, 1], [2, 3, 4]]):
        for st_grp in [[2, 3], [1, 2], [3, 4], list(range(1, 5))]: 
            tmp = np.ma.masked_array(data, mask=(((dbh.mask) | (data.mask)) | (~np.isin(dbh, dbh_grp)) | (~np.isin(supertype, st_grp)) & (harvest != 1)))
            smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
            tmp = np.ma.masked_array(data, mask=(((dbh.mask) | (data.mask)) | (~np.isin(dbh, dbh_grp)) | (~np.isin(supertype, st_grp)) & (harvest == 1)))
            smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)

    return np.ma.masked_array(smooth, (smooth == -9999) | (smooth == 0), fill_value=-9999)


def multilayer_isolation_sieve(data, dbh, harvest, area=404):
    """
    Perform a sieve masking out areas of DBH bins. That is, low-valued classes get treated separate
    from the high-valued classes, both within and outwith harvest areas.
    """
    smooth = np.ma.masked_array(data, mask=(data == -9999), dtype=np.int32)
    
    for dbh_grp in [[0, 1], [2, 3, 4]]:
        for multi_grp in [0,1]:
            tmp = np.ma.masked_array(smooth, mask=(((dbh.mask) | (data.mask)) | (~np.isin(dbh, dbh_grp)) | (~np.isin(data, multi_grp ))) & (harvest != 1)) # Sieve inside harvest areas
            smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
            tmp = np.ma.masked_array(smooth, mask=(((dbh.mask) | (data.mask)) | (~np.isin(dbh, dbh_grp)) | (~np.isin(data, multi_grp ))) & (harvest == 1)) # Sieve outside harvest areas
            smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)

    return smooth


def dbh_st_isolation_sieve(data, area=404):
    """
    Perform a sieve masking out areas of DBH bins and supertype classes after they have been stratified.
    That is, low-valued DBH classes get treated separate from the high-valued classes, both within and outwith harvest areas.
    Simultaneously we create groupings in the mixed, hardwood, and softwood areas of supertype.

    Raster must have values where
    100 is harvest/not harvest
    10-40 is supertype
    0-5 is DBH
    """
    smooth = np.ma.masked_array(data, mask=(data == -9999), dtype=np.int32)
    
    dbh = smooth % 10
    supertype = ((smooth / 10).astype(np.int16) % 10).filled(-9999)
    supertype[supertype == 0] = -9999
    harvest = (smooth / 100).astype(np.int16) % 10

    for dbh_grp in tqdm([[0, 1], [2, 3, 4]]):
        for st_grp in [[2, 3], [1, 2], [3, 4], list(range(1, 5))]:

            tmp = np.ma.masked_array(smooth, mask=(data == -9999) | ((~np.isin(dbh, dbh_grp)) | (~np.isin(supertype, st_grp))) & (harvest != 1))
            smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
            tmp = np.ma.masked_array(smooth, mask=(data == -9999) | ((~np.isin(dbh, dbh_grp)) | (~np.isin(supertype, st_grp))) & (harvest == 1))
            smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)     
    
    return np.ma.masked_array(smooth, smooth == -9999, fill_value=-9999)


def dbh_st_cc_isolation_sieve(data, area=404):
    """
    Perform a sieve masking out areas of DBH bins and supertype classes after they have been stratified.
    That is, low-valued DBH classes get treated separate from the high-valued classes, both within and outwith harvest areas.
    Simultaneously we create groupings in the mixed, hardwood, and softwood areas of supertype.
    Canopy bins for under 60% are separate to above 60%

    Raster must have values where
    1000 is harvest/not harvest
    100-400 is supertype
    00-50 is DBH
    0-5 is CC
    """
    smooth = np.ma.masked_array(data, mask=(data == -9999), dtype=np.int32)
    
    cc = smooth % 10
    dbh = (smooth / 10).astype(np.int16) % 10
    supertype = ((smooth / 100).astype(np.int16) % 10).filled(-9999)
    supertype[supertype == 0] = -9999
    harvest = (smooth / 1000).astype(np.int16) % 10

    for dbh_grp in tqdm([[0, 1], [2, 3, 4]]):
        for st_grp in [[2, 3], [1, 2], [3, 4], list(range(1, 5))]: 
            for cc_group in [[0, 1, 2], [3, 4], list(range(5))]:
                tmp = np.ma.masked_array(smooth, mask=(data == -9999) | ((~np.isin(dbh, dbh_grp)) | (~np.isin(supertype, st_grp)) | (~np.isin(cc, cc_group))) & (harvest != 1))
                smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
                tmp = np.ma.masked_array(smooth, mask=(data == -9999) | ((~np.isin(dbh, dbh_grp)) | (~np.isin(supertype, st_grp)) | (~np.isin(cc, cc_group))) & (harvest == 1))
                smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)

#    tmp = np.ma.masked_array(smooth, mask=(data == -9999) | (harvest == 1))
#    smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
    
    return np.ma.masked_array(smooth, smooth == -9999, fill_value=-9999)


def dbh_st_cc_multi_isolation_sieve(data, area=404):
    """
    Perform a sieve masking out areas of DBH bins and supertype classes after they have been stratified.
    That is, low-valued DBH classes get treated separate from the high-valued classes, both within and outwith harvest areas.
    Simultaneously we create groupings in the mixed, hardwood, and softwood areas of supertype.
    Canopy bins for under 60% are separate to above 60%

    Raster must have values where
    10000 is multilayer/not multilayer
    1000 is harvest/not harvest
    100-400 is supertype
    00-50 is DBH
    0-5 is CC
    """
    smooth = np.ma.masked_array(data, mask=(data == -9999), dtype=np.int32)
    
    cc = smooth % 10
    dbh = (smooth / 10).astype(np.int16) % 10
    supertype = ((smooth / 100).astype(np.int16) % 10).filled(-9999)
    supertype[supertype == 0] = -9999
    harvest = (smooth / 1000).astype(np.int16) % 10
    multilayer = (smooth / 10000).astype(np.int16) % 10

    for dbh_grp in tqdm([[0, 1], [2, 3, 4]]):
        for st_grp in [[2, 3], [1, 2], [3, 4], list(range(1, 5))]: 
            for cc_group in [[0, 1, 2], [3, 4], list(range(5))]:
                for multilayer_group in [0, 1]:
                    tmp = np.ma.masked_array(smooth, mask=(data == -9999) | ((~np.isin(dbh, dbh_grp)) | (~np.isin(supertype, st_grp)) | (~np.isin(cc, cc_group)) | (~np.isin(multilayer, multilayer_group))) & (harvest != 1))
                    smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
                    tmp = np.ma.masked_array(smooth, mask=(data == -9999) | ((~np.isin(dbh, dbh_grp)) | (~np.isin(supertype, st_grp)) | (~np.isin(cc, cc_group)) | (~np.isin(multilayer, multilayer_group))) & (harvest == 1))
                    smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)

#    tmp = np.ma.masked_array(smooth, mask=(data == -9999) | (harvest == 1))
#    smooth = smooth_tif(tmp, mask=~tmp.mask, area=area)
    
    return np.ma.masked_array(smooth, smooth == -9999, fill_value=-9999)



def dbh_st_zonal(gdf, st_data, hv_data, dbh_data, transform, outdir, save=None):
    """
    Run zonal stats over level-1 stratification, replacing values of supertype, harvest, and DBH columns.
    """
    gdf['STRAT'] = gdf['val']

    gdf['DBH'] = gdf['STRAT'].astype(int) % 10
    gdf['ST'] = (gdf['STRAT'] / 10).astype(int) % 10
    gdf['HV'] = (gdf['STRAT'] / 100).astype(int)
    gdf['STRAT'] = gdf['STRAT'].astype(int)

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, st_data, affine=transform, geojson_out=True, stats='majority', nodata=-9999), crs=gdf.crs)
    gdf['ST_ZON'] = gdf['majority'].fillna(0).astype(np.int16)
    del gdf['majority']

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, hv_data, affine=transform, geojson_out=True, stats='majority', nodata=-9999), crs=gdf.crs)
    gdf['HV_ZON'] = gdf['majority'].astype(np.int16)
    del gdf['majority']

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, dbh_data, affine=transform, geojson_out=True, stats='mean'), crs=gdf.crs)
    gdf['DBH_ZON'] = gdf['mean']
    del gdf['mean']

    gdf = gdf[gdf['ST_ZON'] > 0]

    gdf.loc[gdf['ST_ZON'] > 2, 'DBH_B_ZON'] = np.digitize(gdf.loc[gdf['ST_ZON'] > 2, 'DBH_ZON'], dbh_sw_bins)
    gdf.loc[gdf['ST_ZON'] < 3, 'DBH_B_ZON'] = np.digitize(gdf.loc[gdf['ST_ZON'] < 3, 'DBH_ZON'], dbh_hw_bins)
    gdf.loc[gdf['DBH_ZON'] <= 0, 'DBH_B_ZON'] = 0

    gdf['STRAT_ZON'] = (gdf['HV_ZON'] * 100) + (gdf['ST_ZON'] * 10) + (gdf['DBH_B_ZON'] * 1)
    if save:
        gdf.to_file(f'{outdir}/temp_dbh_st_zonal.gpkg', driver='GPKG')
    return gdf


def dbh_st_cc_zonal(gdf, st_data, hv_data, dbh_data, cc_data, transform, outdir, save=None):
    """
    Run zonal stats over level-2 stratification, replacing values of supertype, harvest, CC, and DBH columns.
    """
    gdf['STRAT'] = gdf['val']
    gdf['CC'] = gdf['STRAT'].astype(int) % 10
    gdf['DBH'] = (gdf['STRAT'] / 10).astype(int) % 10
    gdf['ST'] = (gdf['STRAT'] / 100).astype(int) % 10
    gdf['HV'] = (gdf['STRAT'] / 1000).astype(int) % 10
    gdf['STRAT'] = gdf['STRAT'].astype(int)

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, st_data, affine=transform, geojson_out=True, stats='majority', nodata=-9999), crs=gdf.crs)
    gdf['ST_ZON'] = gdf['majority'].fillna(0).astype(np.int16)
    del gdf['majority']

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, hv_data, affine=transform, geojson_out=True, stats='majority', nodata=-9999), crs=gdf.crs)
    gdf['HV_ZON'] = gdf['majority'].fillna(0).astype(np.int16)
    del gdf['majority']

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, cc_data, affine=transform, geojson_out=True, stats='mean', nodata=-9999), crs=gdf.crs)
    gdf['CC_ZON'] = gdf['mean'].fillna(0).astype(np.int16)
    del gdf['mean']

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, dbh_data, affine=transform, geojson_out=True, stats='mean'), crs=gdf.crs)
    gdf['DBH_ZON'] = gdf['mean']
    del gdf['mean']

    gdf.loc[gdf['ST_ZON'] > 2, 'DBH_B_ZON'] = np.digitize(gdf.loc[gdf['ST_ZON'] > 2, 'DBH_ZON'], dbh_sw_bins)
    gdf.loc[gdf['ST_ZON'] < 3, 'DBH_B_ZON'] = np.digitize(gdf.loc[gdf['ST_ZON'] < 3, 'DBH_ZON'], dbh_hw_bins)
    gdf.loc[gdf['DBH_ZON'] <= 0, 'DBH_B_ZON'] = 0

    gdf['CC_B_ZON'] = np.digitize(gdf['CC_ZON'], cc_bins)

    gdf['STRAT_ZON'] = (gdf['HV_ZON'] * 1000) + (gdf['ST_ZON'] * 100) + (gdf['DBH_B_ZON'] * 10) + gdf['CC_B_ZON']

    if save:
        gdf.to_file(f'{outdir}/temp_dbh_st_cc_zonal.gpkg', driver='GPKG')

    return gdf


def dbh_st_cc_multi_zonal(gdf, st_data, hv_data, dbh_data, cc_data, multilayer_binned, transform, outdir, save=None):
    """
    Run zonal stats over level-2 stratification, replacing values of supertype, harvest, CC, and DBH columns.

    This function now includes the multilayer predictions

    """
    gdf['STRAT'] = gdf['val']
    gdf['CC'] = gdf['STRAT'].astype(int) % 10
    gdf['DBH'] = (gdf['STRAT'] / 10).astype(int) % 10
    gdf['ST'] = (gdf['STRAT'] / 100).astype(int) % 10
    gdf['HV'] = (gdf['STRAT'] / 1000).astype(int) % 10
    gdf['Multilayer'] = (gdf['STRAT'] / 10000).astype(int) % 10
    gdf['STRAT'] = gdf['STRAT'].astype(int)

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, st_data, affine=transform, geojson_out=True, stats='majority', nodata=-9999), crs=gdf.crs)
    gdf['ST_ZON'] = gdf['majority'].fillna(0).astype(np.int16)
    del gdf['majority']

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, hv_data, affine=transform, geojson_out=True, stats='majority', nodata=-9999), crs=gdf.crs)
    gdf['HV_ZON'] = gdf['majority'].fillna(0).astype(np.int16)
    del gdf['majority']

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, cc_data, affine=transform, geojson_out=True, stats='mean', nodata=-9999), crs=gdf.crs)
    gdf['CC_ZON'] = gdf['mean'].fillna(0).astype(np.int16)
    del gdf['mean']

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, dbh_data, affine=transform, geojson_out=True, stats='mean'), crs=gdf.crs)
    gdf['DBH_ZON'] = gdf['mean']
    del gdf['mean']

    gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, multilayer_binned, affine=transform, geojson_out=True, stats='majority'), crs=gdf.crs)
    gdf['ML_ZON'] = gdf['majority']
    del gdf['majority']


    gdf.loc[gdf['ST_ZON'] > 2, 'DBH_B_ZON'] = np.digitize(gdf.loc[gdf['ST_ZON'] > 2, 'DBH_ZON'], dbh_sw_bins)
    gdf.loc[gdf['ST_ZON'] < 3, 'DBH_B_ZON'] = np.digitize(gdf.loc[gdf['ST_ZON'] < 3, 'DBH_ZON'], dbh_hw_bins)
    gdf.loc[gdf['DBH_ZON'] <= 0, 'DBH_B_ZON'] = 0

    gdf['CC_B_ZON'] = np.digitize(gdf['CC_ZON'], cc_bins)
    gdf['STRAT_ZON'] = 0
    gdf['STRAT_ZON'] = (gdf['ML_ZON'] * 10000) + (gdf['HV_ZON'] * 1000) + (gdf['ST_ZON'] * 100) + (gdf['DBH_B_ZON'] * 10) + gdf['CC_B_ZON']

    if save:
        gdf.to_file(f'{outdir}/temp_dbh_st_cc_multi_zonal.gpkg', driver='GPKG')

    #gdf = gdf.dissolve(by='STRAT_ZON')
    # Explode to individual polygons again, and reset the index
    #gdf = gdf.reset_index().explode(drop=True, reset_index=True)

    return gdf





def spc_zonal(gdf, spc_dir, fia_file, save=None):
    """
    Run zonal stats over a list of species rasters. Make sure pc is in the name. 
    There must be a CSV file with a lookup between USFS names and LV/FIA codes.
    """
    df = gpd.pd.read_csv(fia_file)
    files = glob(spc_dir)
    for fi in tqdm(files):
        tif = rio.open(fi)
        data = tif.read(1, masked=True)
        spc = re.search('\w+(?=_pc)', os.path.basename(fi)).group(0)
        lv_code = df[df['USFS_Code'] == spc]['LV_CODE'].values[0]
        if isinstance(lv_code, str):
            gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, data, affine=tif.transform, geojson_out=True, stats='mean', nodata=-9999), crs=gdf.crs)
            gdf[lv_code] = gdf['mean'].fillna(0)
            del gdf['mean']
    gdf['StandID'] = gdf.index

    if save:
        gdf.to_file('merging/lv_strat_allzon_5ac_t6')
    return gdf


def load_data(path, props=None):
    """
    Read the data as either
        - a data set that has already been clipped
        - a full-tile image that has a property to clip by
    Can be used on an entire area level but use at your own risk.
    """
    tif = rio.open(path)
    profile = tif.profile.copy()
    if props:
        prop_df = gpd.read_file(props).to_crs(tif.crs)
        prop_df.geometry = prop_df.buffer(20)
        data, transform = mask(tif, prop_df.geometry, nodata=-9999, crop=True)
        data = np.ma.masked_array(data[0], mask=data[0] == -9999, fill_value=-9999)
    else:
        data = tif.read(1, masked=True)
        transform = tif.transform

    profile.update({
        "transform": transform,
        "height": data.shape[0],
        "width": data.shape[1],
        "nodata": -9999,
        "driver": "GTiff"        
    })
    return data, profile


def prep_supertype(data, filter_size=3, radius=1):
    data[data == 1] = 5
    st_clean = erosion(median_filter(data, filter_size), disk(radius))
    st_clean[st_clean == 5] = 1
    return np.ma.masked_array(st_clean, (st_clean == -9999) | (st_clean == 0), fill_value=-9999)


# the function below was used to test various ways of trying to improve DBH smoothing results but to no avail
'''
def prep_dbh(data, supertype):
    dbh_hwsw = np.zeros(data.shape, dtype=np.int16)
    dbh_hwsw[supertype < 3] = np.digitize(data[supertype < 3], dbh_hw_bins) # error here with wrong bins.
    dbh_hwsw[supertype > 2] = np.digitize(data[supertype > 2], dbh_sw_bins) # errro here with wrong bins. 
    dbh_hwsw = median_filter(dbh_hwsw, 5)
    dbh_hwsw = dilation(dbh_hwsw)
    #dbh_hwsw = median_filter(dbh_hwsw, 3)
    return np.ma.masked_array(dbh_hwsw, data.mask, fill_value=-9999)

'''
def prep_dbh(data, supertype):
    dbh_hwsw = np.zeros(data.shape, dtype=np.int16)
    dbh_hwsw[supertype < 3] = np.digitize(data[supertype < 3], dbh_hw_bins) # error here with wrong bins.
    dbh_hwsw[supertype > 2] = np.digitize(data[supertype > 2], dbh_sw_bins) # errro here with wrong bins. 
    dbh_hwsw = median_filter(dbh_hwsw, 3)
    return np.ma.masked_array(dbh_hwsw, data.mask, fill_value=-9999)


def ac_to_px(area):
    return round((area * 4046.86) / 100)


def load_harvest(hv_path, props=None, profile=None, column_name=None):
    if 'tif' in hv_path:
        return load_data(hv_path, props)
    else:
        left = profile['transform'][2]
        bottom = profile['transform'][5] - (profile['height']  * 10)
        right = profile['transform'][2] + (profile['width'] * 10)
        top = profile['transform'][5]

        hv_gdf = gpd.read_file(hv_path).to_crs(profile['crs'])
        #hv_gdf.geometry = hv_gdf.buffer(10) # add 10 metres so that when rasterized it covers enough ground
        hv_gdf['hv'] = 1
        hv_gdf = hv_gdf[hv_gdf.within(geometry.box(*[left, bottom, right, top]))]

        hv = rasterize(hv_gdf, profile, 'hv', fill=0)
        return hv, profile['transform']


def save_raster(data, profile, name):
    profile['dtype'] = np.int16
    profile['nodata'] = -9999
    with rio.open(name , 'w', **profile) as dst:
        dst.write(data.filled(-9999), 1)


def clean_ct(data, profile, outname):
    out_tif = outname.replace('gpkg', 'tif')
    ds_name = os.path.basename(outname).replace('.gpkg', '')
    med = median_filter(data, 3)
    smooth = smooth_tif(med, mask=~(data == -9999), area = 202)
    smooth = np.ma.masked_array(smooth, (smooth == -9999), fill_value=-9999)
    save_raster(smooth, profile, out_tif)
    call(f'gdal_polygonize.py {out_tif} {outname} -b 1 -f "GPKG" {ds_name} val', shell=True)
#    call(f'rio shapes {out_tif} -o {outname}', shell=True)


@click.command()
@click.option('--dbh', default=None, type=click.Path(exists=True), help='layer')
@click.option('--st', default=None, type=click.Path(exists=True), help='layer')
@click.option('--hv', default=None, type=click.Path(exists=True), help='layer')
@click.option('--props', default=None, type=click.Path(exists=True), help='shapefile of properties')
@click.option('--min_area', type=int, default=5, help='Minimum are size (in Acres).')
@click.option('--outdir', type=click.Path(), default='outputs')
def main(dbh, st, hv, props, min_area, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    tif = rio.open(st)

    dbh_data, profile = load_data(dbh, props)
    print('DBH data loaded')
    st_data, _ = load_data(st, props)
    st_data.mask = (st_data.mask) | (st_data == 0)
    print('Supertype data loaded')

    hv_data, _ = load_harvest(hv, props=props, profile=profile)
    print('Harvest data loaded')
    save_raster(hv_data, profile, f"{outdir}/hv_mask.tif")
        
    dbh_ma = prep_dbh(dbh_data, st_data)

    save_raster(dbh_ma, profile, f'{outdir}/dbh_prep.tif')

    dbh_smooth = dbh_isolation_sieve(dbh_ma.astype(np.int16), st_data, hv_data, ac_to_px(min_area))
    print('DBH data cleaned')

    st_clean = prep_supertype(st_data.copy())
    print('Supertype data prepared')

    save_raster(st_clean, profile, f'{outdir}/st_prep.tif')

    dbh_smooth = np.ma.masked_array(dbh_smooth, (st_clean.mask), fill_value=-9999)


    save_raster(dbh_smooth, profile, f'{outdir}/dbh_smooth.tif')

    print('done')

    '''
        hwsw_smooth = st_isolation_sieve(st_clean, dbh_smooth, hv_data, area=ac_to_px(2.5))
        print('Supertype data cleaned')

        hwsw_smooth = (hwsw_smooth * 10) + dbh_smooth
        hwsw_smooth[hv_data == 1] += (hv_data.data[hv_data == 1] * 100)

        save_raster(hwsw_smooth, profile, f'{outdir}/hwsw_dbh_smooth.tif')

        strat_smooth = dbh_st_isolation_sieve(hwsw_smooth, area=ac_to_px(2.5))
        print('DBH Supertype data cleaned')
        strat_smooth.mask = (strat_smooth.data == 0) | (strat_smooth.data == -9999)

        save_raster(strat_smooth, profile, f'{outdir}/strat_smooth.tif')

        p = call(f'gdal_polygonize.py {outdir}/strat_smooth.tif {outdir}/strat_l1.gpkg -b 1 -f "GPKG" strat_l1 val', shell=True)
        print('Level 1 delineations created')

    if not os.path.exists(f'{outdir}/strat_l2.gpkg'):
        gdf = gpd.read_file(f'{outdir}/strat_l1.gpkg').to_crs(tif.crs)
        gdf = dbh_st_zonal(gdf, st_data, hv_data, dbh_data, profile['transform'], outdir, save)
        print("Zonal stats calculated")

        strat_zon = rasterize(gdf, profile)
        print("Zonal stats rasterized")

        if save:
            save_raster(strat_zon, profile, f'{outdir}/strat_zon_l1.tif')
        print(f"Cleaning level 1 zonal delineation to 2.5 acres")

        smooth_zon = dbh_st_isolation_sieve(strat_zon, area=ac_to_px(2.5))
        if save:
            save_raster(smooth_zon, profile, f'{outdir}/smooth_zon.tif')
        
        cc_bin = np.digitize(cc_data, cc_bins)
        cc_med = median_filter(cc_bin, 3)
        cc_med = np.ma.masked_array(cc_med, cc_med == -9999, fill_data=-9999)

        dbh_smooth, _ = load_data(f'{outdir}/dbh_smooth.tif', props)

        print('CC data prepared')
        cc_smooth = st_isolation_sieve(cc_med.astype(np.int16), dbh_smooth, hv_data, ac_to_px(2.5))
        print("CC data cleaned to 5 acres")
        cc_ma = np.ma.masked_array(cc_smooth, mask=(smooth_zon.mask) | (cc_smooth == -9999), fill_value=-9999)

        if save:
            save_raster(cc_ma, profile, f'{outdir}/cc_smooth.tif')
        zon_cc = (smooth_zon * 10) + cc_ma
        print("Level 2 strat prepared, cleaning")

        if save:
            save_raster(zon_cc, profile, f'{outdir}/zon_cc.tif')

        smooth_zon_cc = dbh_st_cc_isolation_sieve(zon_cc, area=ac_to_px(2.5)) 
        #smooth_zon_cc = dbh_isolation_sieve(zon_cc, st_data, hv_data, area=ac_to_px(4))
        
        if save:
            save_raster(smooth_zon_cc, profile, f'{outdir}/smooth_zon_cc.tif')
        print(f"Level 2 strat cleaned to {min_area} acres")

        smooth_zon_cc, _ = load_data(f'{outdir}/smooth_zon_cc.tif', props)
        smooth_zon_cc[smooth_zon_cc < 0] = -9999
        #smooth = simple_smooth(smooth_zon_cc, hv_data, area=ac_to_px(5)) # changing 2.5 to 5ac doesn't seem to make any difference
        #smooth = smooth_tif(smooth_zon_cc, ~smooth_zon_cc.mask, area=ac_to_px(2.5)) # again, no difference changing from 2.5 to 5 ac
        #smooth = dbh_isolation_sieve(smooth_zon_cc, st_data, hv_data, area=ac_to_px(2.5))
        #smooth = np.ma.masked_array(smooth, smooth==-9999, fill_value=-9999)
        save_raster(smooth_zon_cc, profile, f'{outdir}/strat_cc_smooth.tif')

        p = call(f'gdal_polygonize.py {outdir}/strat_cc_smooth.tif {outdir}/strat_l2.gpkg -b 1 -f "GPKG" strat_l2 val', shell=True)
        print("Level 2 delineations created")

        gdf = gpd.read_file(f'{outdir}/strat_l2.gpkg').to_crs(tif.crs)
        gdf = dbh_st_cc_zonal(gdf, st_data, hv_data, dbh_data, cc_data, profile['transform'], outdir, save)
        strat_zon = rasterize(gdf, profile)
        if save:
            save_raster(strat_zon, profile, f'{outdir}/strat_zon_l2.tif')
        gdf['StandID'] = gdf.index
        gdf.to_file(f'{outdir}/strat_l3.gpkg', driver='GPKG')


        ####    ADDITION OF MULTILAYER STRATIFICATION   ####

        multilayer_bins = [200]

        multilayer, _= load_data(multilayer_path, props)
        multilayer_bin = np.digitize(multilayer, multilayer_bins)
        multilayer_med = median_filter(multilayer_bin, 3)
        multilayer_med = np.ma.masked_array(multilayer_med, multilayer_med == -9999, fill_data=-9999)
        multilayer_smooth = multilayer_isolation_sieve(multilayer_med.astype(np.int16), dbh_smooth, hv_data, ac_to_px(2.5))
        print("Multilayer data cleaned to 2.5 acres")
        multilayer_ma = np.ma.masked_array(multilayer_smooth, mask=((smooth_zon_cc.mask) | (multilayer_smooth == -9999)), fill_value=-9999)
        
        if save:
            save_raster(multilayer_ma, profile, f'{outdir}/multilayer_smooth.tif')
        
        #multilayer_cc_st_dbh = (strat_zon * 10) + multilayer_ma
        multilayer_cc_st_dbh = (10000 * multilayer_ma) + strat_zon

        smooth_multilayer_cc_st_dbh = dbh_st_cc_multi_isolation_sieve(multilayer_cc_st_dbh, area=ac_to_px(2.5))
        smooth_multilayer_cc_st_dbh[smooth_multilayer_cc_st_dbh < 0] = -9999


        # this dbh_blocks bit has been added as an attempt to get the DBH boundaries between classes 1/2 and 3/4/5 (i.e. 0/1 and 2/3/4 in the code)
        # to be respected in the final output
        dbh_blocks = dbh_smooth.copy()

        dbh_blocks[dbh_smooth < 2] = 0
        dbh_blocks[dbh_smooth > 1] = 10
        dbh_blocks = np.ma.masked_array(dbh_blocks, mask=dbh_smooth.mask, fill_value=0)

        smooth_multilayer_cc_st_dbh += dbh_smooth

        smooth = dbh_isolation_sieve(smooth_multilayer_cc_st_dbh, st_data, hv_data, area=ac_to_px(5))
        #smooth = simple_smooth(smooth_multilayer_cc_st_dbh, hv_data, area=ac_to_px(5))
        #smooth = smooth_tif(smooth_multilayer_cc_st_dbh, ~smooth_multilayer_cc_st_dbh.mask, area=ac_to_px(5))
        smooth_multilayer_cc_st_dbh = np.ma.masked_array(smooth, smooth==-9999, fill_value=-9999)

    
       # Write out raster of smoothed stratified version of DBH, CC, ST and Multilayer canopy data
        save_raster(smooth_multilayer_cc_st_dbh, profile, f'{outdir}/smooth_multilayer_cc_st_dbh.tif')

        # crop with a little bit extra around the edges so that final clipping to property boundaries works more cleanly

        p = call(f'gdalwarp -crop_to_cutline -cutline {props} -dstnodata -9999 -wo CUTLINE_ALL_TOUCHED=TRUE -wo SOURCE_EXTRA=2 -co COMPRESS=DEFLATE {outdir}/smooth_multilayer_cc_st_dbh.tif {outdir}/smooth_multilayer_cc_st_dbh_recrop.tif', shell=True)

        # Turn the raster above into a polygon layer ready for zonal statistics
        p = call(f'gdal_polygonize.py {outdir}/smooth_multilayer_cc_st_dbh_recrop.tif {outdir}/strat_l2_multilayer.gpkg -b 1 -f "GPKG" strat_l2_multilayer val', shell=True)
        print("Level 2 multilayer delineations created")

    #if not os.path.exists(f'{outdir}/strat_l3_multilayer.gpkg'):

        ####    ZONAL STATS SECTION FOR MULTILAYER CANOPY   ####
        print('Running zonal stats to create level-3 multilayer stratification')
        
        gdf = gpd.read_file(f'{outdir}/strat_l2_multilayer.gpkg').to_crs(tif.crs)
        gdf = dbh_st_cc_multi_zonal(gdf, st_data, hv_data, dbh_data, cc_data, multilayer_bin, profile['transform'], outdir, save)

        #strat_zonal_stats = rasterize(gdf, profile)
        #save_raster(strat_zonal_stats, profile, f'{outdir}/strat_zonal_l2_multilayer.tif')
        
        #p = call(f'gdal_polygonize.py {outdir}/strat_zonal_l2_multilayer.tif {outdir}/strat_l3_multilayer.gpkg -b 1 -f "GPKG" strat_l3_multilayer val', shell=True)
        #print("Level 3 multilayer delineations created")
        
        gdf['StandID'] = gdf.index
        gdf.to_file(f'{outdir}/strat_l3_multilayer.gpkg', driver='GPKG')

'''
'''
    if not os.path.exists(f'{outdir}/strat_l3_multilayer_ct.gpkg'):
        gdf = gpd.read_file(f'{outdir}/strat_l3_multilayer.gpkg')
        ct_gdf = gpd.read_file(f'{outdir}/ct_shapes.gpkg').to_crs(gdf.crs)
        print('Covertype shapes loaded')
        gdf.geometry = gdf.buffer(0)
        ct_gdf.geometry = ct_gdf.buffer(0)
        print('Shapes buffered.')
        ct_gdf['CT'] = ct_gdf['val']
        print('Combining Strat shaped with CT shapes (this may take a while...)')
        ov_gdf = gpd.overlay(gdf, ct_gdf, how='union')
        ov_gdf.to_file(f'{outdir}/strat_l3_multilayer_ct.gpkg', driver='GPKG')

    if not os.path.exists(f'{outdir}/strat_l4.gpkg'):
        gdf = gpd.read_file(f'{outdir}/strat_l3_ct.gpkg')
        clean_gdf = process(gdf)
        clean_gdf['CTStandID'] = clean_gdf.index
        clean_gdf.to_file(f'{outdir}/strat_l4.gpkg', driver='GPKG')
    
    if not os.path.exists(f'{outdir}/strat_l4_spc.gpkg', driver='GPKG'):
        gdf = gpd.read_file(f'{outdir}/strat_l4.gpkg')
        spc_gdf = gpd.read_file(f'{outdir}/spc_shapes.gpkg').to_crs(gdf.crs)
        print('Species shapes loaded')

        gdf.geometry = gdf.buffer(0)
        spc_gdf.geometry = spc_gdf.buffer(0)
        print('Shapes buffered.')

        spc_gdf['SPC'] = spc_gdf['val']
        print('Combining Strat shaped with SPC shapes (this may take a while...)')

        ov_gdf = gpd.overlay(gdf, spc_gdf, how='union')
        ov_gdf.to_file(f'{outdir}/strat_l4_spc.gpkg', driver='GPKG')

    if not os.path.exists(f'{outdir}/strat_l5.gpkg'):
        gdf = gpd.read_file(f'{outdir}/strat_l4_spc.gpkg')
        clean_gdf = process(gdf, group_col='CTStandID', elim_col='SPC')
        clean_gdf['SPCStandID'] = clean_gdf.index
        clean_gdf.to_file(f'{outdir}/strat_l5.gpkg', driver='GPKG')

    if not os.path.exists(f'{outdir}/strat_l6.gpkg'):
        pass
'''
if __name__ == "__main__":
    main()
