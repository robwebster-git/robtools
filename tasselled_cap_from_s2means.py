import click
import os
import rasterio as rio
import numpy as np
from glob import glob
import sys
from rasterio.warp import calculate_default_transform, reproject, Resampling


def resample_raster_10(file_name, outfile = None):
    print(f'Processing : {file_name}')
    with rio.open(file_name) as src:
        dst_crs = src.crs['init']
        data_type = src.meta['dtype']
        #print('change np.unit16 to np.', data_type)
        transform, width, height = calculate_default_transform(
                src.crs, dst_crs, src.width, src.height, *src.bounds, resolution = (10, 10))
        kwargs = src.meta.copy()
        kwargs.update({
                'crs': src.crs,
                'transform': transform,
                'width': width,
                'height': height
                })
    
        #band_10m = np.zeros((10980, 10980), np.uint16) 
        band_10m = np.full((10980, 10980), src.nodata, np.int16) 
        ## Consider reading the dimensions of target file as fed in by another agument to avoid hardcoding these
        
        for i in range(1, src.count + 1):
            reproject(
                    source=rio.band(src, i),
                    destination = band_10m, # put as array of correct dims 
                    src_transform=src.transform,
                    src_crs=src.crs,
                    src_nodata=src.nodata,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    dst_nodata=src.nodata,
                    resampling=Resampling.nearest)

        # Create masked array
        band_10m=np.ma.masked_array(band_10m, mask=(band_10m==-9999),fill_value=-9999) 

        if  outfile is not None:
            with rio.open(outfile, 'w', **kwargs) as dst:
                dst.write(band_10m.filled(-9999), 1)                      

    return band_10m

# calculate indicies 
def write_raster(data, tif, iname, outname, scaled):
    print(outname)
    meta = rio.open(tif).meta.copy()
    data, scl_off = apply_scale_and_offset(iname, scaled, data)
    meta.update({"driver": "GTiff", 
                 "height": data.shape[0],
                 "width": data.shape[1],
                 "nodata": -9999,
                 "dtype": data.dtype,
                 "compress": "deflate",
                 "tiled": True})
    with rio.open(outname, 'w', **meta) as dst:
        dst.write(data.filled(-9999), 1)
    # Add scale and offset as metadata
    arglist=['gdal_edit.py', '-scale', str(scl_off['scale']), '-offset', str(scl_off['offset']), outname] 
    res=gdal_edit.gdal_edit(arglist)

def set_output_dtype(scaled):
    """Return dtype for outputs based on scaled flag value."""
    dtype_dict = {True:np.int16, False:np.float32}
    return dtype_dict[scaled]

def set_scale_and_offset(index_name, scaled):
    """Return scale factor and offset to apply to the index."""
    scale_offset_dict = {'ndvi':         {'scale': 10000.0, 'offset': 10000.0},
                         'savi':         {'scale': 10000.0, 'offset': 10000.0},
                         'msavi':        {'scale': 10000.0, 'offset': 10000.0},
                         'lai':          {'scale': 10000.0, 'offset': 10000.0},
                         'evi':          {'scale': 10000.0, 'offset': 10000.0},
                         'dvi':          {'scale': 10000.0, 'offset': 10000.0},
                         'cirededge':    {'scale': 3000.0, 'offset': 3000.0},
                         'ind84':        {'scale': 600.0, 'offset': 0.0},
                         'cigreen':      {'scale': 1000.0, 'offset': 1000.0},
                         'tsavi':        {'scale': 10000.0, 'offset': 10000.0},
                         'wdrvi':        {'scale': 10000.0, 'offset': 10000.0},
                         'wdrvirededge': {'scale': 10000.0, 'offset': 10000.0},
                         'wdrvigreen':   {'scale': 10000.0, 'offset': 10000.0},
                         'msr670':       {'scale': 600.0, 'offset': 600.0},
                         'ndwi':         {'scale': 10000.0, 'offset': 10000.0},
                         'mndwi':        {'scale': 10000.0, 'offset': 10000.0},
                         'nbr':          {'scale': 10000.0, 'offset': 10000.0},
                         'bais2':        {'scale': 10000.0, 'offset': 10000.0},
                         'pssr':         {'scale': 600.0, 'offset': 0.0}
                         }

    scale_offset = scale_offset_dict[index_name]
    if not scaled:
        scale_offset = {'scale': 1.0, 'offset': 0.0}
    return scale_offset

def set_valid_range(index_name):
    """Return scale factor and offset to apply to the index."""
    valid_range_dict = {'ndvi':         {'min': -1.0, 'max': 1.0},
                        'savi':         {'min': -1.0, 'max': 1.0},
                        'msavi':        {'min': -1.0, 'max': 1.0},
                        'lai':          {'min': -1.0, 'max': 2.0},
                        'evi':          {'min': -1.0, 'max': 2.2},
                        'dvi':          {'min': -1.0, 'max': 1.0},
                        'cirededge':    {'min': -1.0, 'max': 9.5},
                        'ind84':        {'min': 0.0, 'max': 50.0},
                        'cigreen':      {'min': -1.0, 'max': 30.0},
                        'tsavi':        {'min': -1.0, 'max': 1.0},
                        'wdrvi':        {'min': -1.0, 'max': 2.2},
                        'wdrvirededge': {'min': -1.0, 'max': 2.2},
                        'wdrvigreen':   {'min': -1.0, 'max': 2.2},
                        'msr670':       {'min': -1.0, 'max': 50.0},
                        'ndwi':         {'min': -1.0, 'max': 1.0},
                        'mndwi':        {'min': -1.0, 'max': 1.0},
                        'nbr':          {'min': -1.0, 'max': 2.0},
                        'bais2':        {'min': -1.0, 'max': 2.0},
                        'pssr':         {'min': 0.0, 'max': 50.0}
                         }

    valid_range = valid_range_dict[index_name]
    return valid_range

def check_scale_and_limits(indices=['ndvi', 'savi', 'msavi', 'lai', 'dvi', 'cirededge', 'ind84', 
               'cigreen', 'tsavi', 'wdrvi', 'msr670', 'ndwi', 'mndwi', 'nbr', 'bais2', 'pssr']):
    """Make sure scaling and limits are OK."""
    scaled = True
    for iname in indices:
        scl_off = set_scale_and_offset(iname, scaled)
        val_rng = set_valid_range(iname)
        print(iname)
        print("SCALE_OFFSET: ",scl_off)
        print("VALID_RANGE: ",val_rng)
        print("MIN: ", (val_rng['min']*scl_off['scale']) + scl_off['offset'])
        print("MAX: ", (val_rng['max']*scl_off['scale']) + scl_off['offset'])   

def clip_to_valid_range(index_name, data):
    """Clip the data to a valid range."""
    val_rng = set_valid_range(index_name)
    mask = (data < val_rng['min']) & ~(data.mask)
    # Set values below min to min
    if mask.any():
        data[mask] = val_rng['min']
    # Set values above max to max
    mask = (data > val_rng['max']) & ~(data.mask)
    if mask.any():
        data[mask] = val_rng['max']
    return data


def apply_scale_and_offset(index_name, scaled, data):
    """Apply scaling and offsetting."""
    if scaled:
        data = clip_to_valid_range(index_name, data)
        scl_off = set_scale_and_offset(index_name, scaled)
        data = (data*scl_off['scale']) + scl_off['offset']
        data = data.round().astype(set_output_dtype(scaled))
    return data, scl_off


def tasseled_cap_tm(rast, reflectance=True):
    '''
    Applies the Tasseled Cap transformation for S2 data. Assumes that the
    data are TM reflectance data (i.e., Landsat Surface Reflectance). The
    coefficients for reflectance factor data are taken from Crist (1985) in
    Remote Sensing of Environment 17:302.
    '''
    if reflectance:
        # Reflectance factor coefficients for S2 bands 1-12 and 8a
        # Rows top down = bands 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 8a
        # Columns L-R = Brightness (TC1), Greenness (TC2), and Wetness (TC3)
        '''
        https://www.researchgate.net/profile/R_Nedkov/publication/329183655_ORTHOGONAL_TRANSFORMATION_OF_SEGMENTED_IMAGES_FROM_THE_SATELLITE_SENTINEL-2_Roumen_Nedkov/links/5c978a9645851506d728189f/ORTHOGONAL-TRANSFORMATION-OF-SEGMENTED-IMAGES-FROM-THE-SATELLITE-SENTINEL-2-Roumen-Nedkov.pdf
        
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
            (0.3882, -0.4578, -0.4064),
            (0.1366, -0.4064, -0.5602),
            (0.4750, 0.3625, -0.1389)
        ], dtype=np.float32)
        '''

        r = np.array([
            (0.0356, 0.0822, 0.1360, 0.2611, 0.2964, 0.3338, 0.3877, 0.3895, 0.0949, 0.3882, 0.1366, 0.4750),
            (-0.0635, -0.1128, -0.1680, -0.3480, -0.3303, 0.0852, 0.3302, 0.3165, 0.0467, -0.4578, -0.4064, 0.3625),
            (0.0649, 0.1363, 0.2802, 0.3072, 0.5288, 0.1379, -0.0001, -0.0807, -0.0302, -0.4064, -0.5602, -0.1389)
        ], dtype=np.float32)

    else:
        raise NotImplemented('Only support for Sentinel 2 reflectance has been implemented')

    shp = rast.shape

    # Can accept either a gdal.Dataset or numpy.array instance
    if not isinstance(rast, np.ndarray):
        x = rast.ReadAsArray().reshape(shp[0], shp[1]*shp[2])

    else:
        x = rast.reshape(shp[0], shp[1]*shp[2])

    dot_product = np.dot(r, x)
    reshaped = dot_product.reshape((3, shp[1], shp[2]))

    return reshaped.round().astype(np.int32)


def calc_tc_components(s2_stack, nodata=-9999):
    '''
    Calculates the first three components of the tasselled cap transformed
    input data - this is expected to be
    a tasseled cap-transformed raster. The NoData value is assumed to be
    negative (could never be the maximum value in a band).
    '''
    shp = s2_stack.shape

    x = s2_stack.reshape(shp[0], shp[1]*shp[2])
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

    print('Before reshaping')
    print(h.shape, v.shape, l.shape)

    h = h.reshape((1, shp[1], shp[2]))
    v = v.reshape((1, shp[1], shp[2]))
    l = l.reshape((1, shp[1], shp[2]))

    print('After reshaping')
    print(h.shape, v.shape, l.shape)

    return h, v, l


# function also resamples bands from 20 to 10m 
def open_bands(s2_dir, write_resampled = False):
    # S2 data are scaled reflectances, scale_factor=10000
    # Indices formulae are based on reflectances
    # - Rescale S2 data and calulate indices from floats
    #scale2reflect = np.float32(0.0001)
    scale2reflect = 1
    # 10m bands
    path_10 = s2_dir
    band2_10m = rio.open(glob(path_10 + '/*B02*')[0])
    b2_10m = band2_10m.read(1, masked = True) * scale2reflect
    band3_10m = rio.open(glob(path_10 + '/*B03*')[0])
    b3_10m = band3_10m.read(1, masked = True) * scale2reflect
    band4_10m = rio.open(glob(path_10 + '/*B04*')[0])
    b4_10m = band4_10m.read(1, masked = True) * scale2reflect
    band8_10m = rio.open(glob(path_10 + '/*B08*')[0])
    b8_10m = band8_10m.read(1, masked = True) * scale2reflect
    # 20m bands resampled
    path_20 = s2_dir    
    b5_10m = resample_raster_10(glob(path_20 + '/*B05*')[0]) * scale2reflect
    b6_10m = resample_raster_10(glob(path_20 + '/*B06*')[0]) * scale2reflect
    b7_10m = resample_raster_10(glob(path_20 + '/*B07*')[0]) * scale2reflect
    b11_10m = resample_raster_10(glob(path_20 + '/*B11*')[0]) * scale2reflect
    b12_10m = resample_raster_10(glob(path_20 + '/*B12*')[0]) * scale2reflect
    b8a_10m = resample_raster_10(glob(path_20 + '/*B8A*')[0]) * scale2reflect
    # 60m bands resampled
    path_60 = s2_dir    
    b1_10m = resample_raster_10(glob(path_60 + '/*B01*')[0]) * scale2reflect
    b9_10m = resample_raster_10(glob(path_60 + '/*B09*')[0]) * scale2reflect
    #b10_10m = resample_raster_10(glob(path_60 + '/*B10*')[0]) * scale2reflect


    return b1_10m, b2_10m, b3_10m, b4_10m, b5_10m, b6_10m, b7_10m, b8_10m, b9_10m, b11_10m, b12_10m, b8a_10m

@click.command()
@click.argument('s2_dir', type=click.Path())
def main(s2_dir):
    print('hello')
    print(s2_dir)
    scene_id = os.path.basename(s2_dir)
    year = scene_id[11:15]
    tile_id = scene_id[39:44]
    path = '/'.join(s2_dir.split('/')[:-2])
    outpath = os.path.join(s2_dir, 'indices')
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    b1_10m, b2_10m, b3_10m, b4_10m, b5_10m, b6_10m, b7_10m, b8_10m, b9_10m, b11_10m, b12_10m, b8a_10m = open_bands(s2_dir)

    #path to band is a path to any 10m band, allows function to aquire correct meta data. 
    band_meta_data = glob(s2fm_dir + '/R10m/*B04*')[0]
    
    with rio.open(band_meta_data) as src:
        profile = src.profile.copy()
    
    aoi_name = rio.open(band_meta_data).name
    tile_date = scene_id.replace('.SAFE', '')
    #calculate indices
    outpath = outpath + f'/{tile_date}'
    print(outpath)

    s2_stack = np.stack([b1_10m, b2_10m, b3_10m, b4_10m, b5_10m, b6_10m, b7_10m, b8_10m, b9_10m, b11_10m, b12_10m, b8a_10m])

    s2_band_sums = s2_stack.sum(axis=0)
    mask_areas = (s2_band_sums / 12).astype(np.int16) == int(profile['nodata'])
    #mask_areas = b1_10m == profile['nodata']
    mask_3d = np.repeat(mask_areas[np.newaxis, :, :], 3, axis=0)

    print(mask_areas.shape)
    print(mask_3d.shape)

    profile['dtype'] = 'float32'
    profile['count'] = 1

    tc_transformed = tasseled_cap_tm(s2_stack)

    tc_masked = np.ma.MaskedArray(tc_transformed, mask=mask_3d, fill_value=profile['nodata'])

    tc1, tc2, tc3 = calc_tc_components(tc_transformed)

    bands_dict = { 1: 'brightness', 2 : 'greeness', 3 : 'wetness'}


    for i, tc in enumerate([tc1, tc2, tc3]):
        with rio.open(f's2_tasselled_cap_TC{i+1}_{bands_dict[i+1]}.tif', 'w', **profile) as dst:
            dst.write(np.ma.MaskedArray(tc, mask=mask_areas, fill_value=profile['nodata']).filled().astype(np.float32))
    '''

    with rio.open('s2_tasselled_cap_band_test2.tif', 'w', **profile) as dst:
        dst.write(tc_masked.filled(profile['nodata']))
    '''

if __name__ == "__main__":
    main()
