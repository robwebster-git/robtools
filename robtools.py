'''
A series of functions to perform common geospatial tasks
that may not be straightforward using existing libraries

Author      : Rob Webster
Date        : October 2020       

'''

import geopandas as gpd
import pandas as pd
import rasterio as rio
from rasterio.io import MemoryFile
from rasterio.mask import mask


def clip(raster_array, raster_profile, shapefile):
    '''
    Take a numpy array with a rasterio profile
    and clip this to using a shapefile

    Parameters
    ----------

    raster_array : ndarray
        An array representing a geospatial raster

    raster_profile : rasterio profile, which is a dict
        Information on how `raster_array` is mapped onto
        the real world

    shapefile : str
        a filename of a shapefile to use
        to clip the raster dataset

    Returns
    -------

    clipped_array : ndarray
        An array clipped to the shapefile

    raster_profile : rasterio profile object, which is a dict
        The information required to use rasterio
        to write a new clipped raster

    '''

    print(f'Cropping array {raster_array} with shapefile {shapefile}')

    #  Create a temporary raster file in memory only
    #  from the input array and profile
    with MemoryFile() as memfile:
        with memfile.open(**raster_profile) as dataset:
            dataset.write(raster_array)

            #  Read shapefile and convert to same CRS
            #  as the raster dataset  
            shp = gpd.read_file(shapefile)
            shp_in_raster_crs = shp.to_crs(raster_profile['crs'])
            
            #  Perform the clip / crop operation
            clipped_array, clipped_transform = mask(dataset, shp_in_raster_crs.geometry, crop=True)
    
    print(f'Shape of cropped data : {clipped_array.shape}')

    #  Update the profile to reflect cropped extents
    #  and transform
    raster_profile['height'] = clipped_array.shape[1]
    raster_profile['width'] = clipped_array.shape[2]
    raster_profile['transform'] = clipped_transform

    return clipped_array, raster_profile
