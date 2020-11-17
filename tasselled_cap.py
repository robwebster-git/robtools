import click
import os
import rasterio as rio
import numpy as np
import glob
import sys


def tasseled_cap_tm(rast, reflectance=True):
    '''
    Applies the Tasseled Cap transformation for S2 data. Assumes that the
    data are TM reflectance data (i.e., Landsat Surface Reflectance). The
    coefficients for reflectance factor data are taken from Crist (1985) in
    Remote Sensing of Environment 17:302.
    '''
    if reflectance:
        # Reflectance factor coefficients for S2 bands 1-12 and 8a
        # Rows top down = bands 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 8a
        # Columns L-R = Brightness (TC1), Greenness (TC2), and Wetness (TC3)
        '''
        https://www.researchgate.net/profile/R_Nedkov/publication/329183655_ORTHOGONAL_TRANSFORMATION_OF_SEGMENTED_IMAGES_FROM_THE_SATELLITE_SENTINEL-2_Roumen_Nedkov/links/5c978a9645851506d728189f/ORTHOGONAL-TRANSFORMATION-OF-SEGMENTED-IMAGES-FROM-THE-SATELLITE-SENTINEL-2-Roumen-Nedkov.pdf
        '''
        r = np.array([
            (0.0356, -0.0635, 0.0649),
            (0.0822, -0.1128, 0.1363),
            (0.1360, -0.1680, 0.2802),
            (0.2611, -0.3480, 0.3072),
            (0.2964, -0.3303, 0.5288),
            (0.3338, 0.0852, 0.1379),
            (0.3877, 0.3302, -0.0001),
            (0.3895, 0.3165, -0.0807),
            (0.0949, 0.0467, -0.0302),
            (0.0009, -0.0009, 0.0003),
            (0.3882, -0.4578, -0.4064),
            (0.1366, -0.4064, -0.5602),
            (0.4750, 0.3625, -0.1389)
        ], dtype=np.float32)

    else:
        raise NotImplemented('Only support for Sentinel 2 reflectance has been implemented')

    shp = rast.shape

    # Can accept either a gdal.Dataset or numpy.array instance
    if not isinstance(rast, np.ndarray):
        x = rast.ReadAsArray().reshape(shp[0], shp[1]*shp[2])

    else:
        x = rast.reshape(shp[0], shp[1]*shp[2])

    return np.dot(r, x).reshape(shp)


def calc_tc_components(rast, nodata=-9999):
    '''
    Calculates the first three components of the tasselled cap transformed
    input data - this is expected to be
    a tasseled cap-transformed raster. The NoData value is assumed to be
    negative (could never be the maximum value in a band).
    '''
    shp = rast.shape

    # Can accept either a gdal.Dataset or numpy.array instance
    if not isinstance(rast, np.ndarray):
        x = rast.ReadAsArray().reshape(shp[0], shp[1]*shp[2])

    else:
        x = rast.reshape(shp[0], shp[1]*shp[2])

    unit = np.ones((1, shp[1] * shp[2]))

    stack = []
    for i in range(0, 3):
        # Calculate the minimum values after excluding NoData values
        tcmin = np.setdiff1d(x[i, ...].ravel(), np.array([nodata])).min()

        # Calculate the normalized TC component TC_i for i in {1, 2, 3}
        stack.append(np.divide(np.subtract(x[i, ...], unit * tcmin),
            unit * (x[i, ...].max() - tcmin)))

    # Unpack the High-albedo, Vegetation, and Low-albedo components
    h, v, l = stack



    return h, v, l


# function also resamples bands from 20 to 10m 
def open_bands(s2_dir, write_resampled = False):
    # S2 data are scaled reflectances, scale_factor=10000
    # Indices formulae are based on reflectances
    # - Rescale S2 data and calulate indices from floats
    scale2reflect = np.float32(0.0001)
    # 10m bands
    path_10 = s2_dir + '/R10m/'
    band2_10m = rio.open(glob(path_10 + '/*B02*')[0])
    b2_10m = band2_10m.read(1, masked = True) * scale2reflect
    band3_10m = rio.open(glob(path_10 + '/*B03*')[0])
    b3_10m = band3_10m.read(1, masked = True) * scale2reflect
    band4_10m = rio.open(glob(path_10 + '/*B04*')[0])
    b4_10m = band4_10m.read(1, masked = True) * scale2reflect
    band8_10m = rio.open(glob(path_10 + '/*B08*')[0])
    b8_10m = band8_10m.read(1, masked = True) * scale2reflect
    # 20m bands resampled
    path_20 = s2_dir + '/R20m/'    
    b5_10m = resample_raster_10(glob(path_20 + '/*B05*')[0]) * scale2reflect
    b6_10m = resample_raster_10(glob(path_20 + '/*B06*')[0]) * scale2reflect
    b7_10m = resample_raster_10(glob(path_20 + '/*B07*')[0]) * scale2reflect
    b11_10m = resample_raster_10(glob(path_20 + '/*B11*')[0]) * scale2reflect
    b12_10m = resample_raster_10(glob(path_20 + '/*B12*')[0]) * scale2reflect
    b8a_10m = resample_raster_10(glob(path_20 + '/*B8A*')[0]) * scale2reflect


    return b2_10m, b3_10m, b4_10m, b8_10m, b5_10m, b6_10m, b7_10m, b11_10m, b12_10m, b8a_10m

@click.command()
@click.argument('s2_dir', type=click.Path())
def main(s2_dir):

    print(s2_dir)
    s2fm_dir = s2_dir + '/S2Fmask'
    scene_id = os.path.basename(s2_dir)
    year = scene_id[11:15]
    tile_id = scene_id[39:44]
    path = '/'.join(s2_dir.split('/')[:-2])
    outpath = os.path.join(s2_dir, 'indices')
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    b2_10m, b3_10m, b4_10m, b8_10m, b5_10m, b6_10m, b7_10m, b11_10m, b12_10m, b8a_10m = open_bands(s2fm_dir)

    #path to band is a path to any 10m band, allows function to aquire correct meta data. 
    band_meta_data = glob(s2fm_dir + '/R10m/*B04*')[0]
    aoi_name = rio.open(band_meta_data).name
    tile_date = scene_id.replace('.SAFE', '')
    #calculate indices
    outpath = outpath + f'/{tile_date}'
    print(outpath)