import numpy as np
import rasterio as rio
import geopandas as gpd
import glob
from rasterio.features import rasterize
from pyproj import CRS



def generate_scene_counts(g_dir):
    '''
    Take in a directory containing geopackages and rasterize each cloud mask using a low resolution.
    From this, calculate the approximate scene count numbers in a given pixel to visualise the distribution
    of clear sky pixels.
    '''

    packages = glob.glob(f'{g_dir}/*.gpkg')

    transform=rio.transform.from_origin(300_000, 6200044, 100, 100) # for tile 30UUG in EPSG:32630, create 100m resolutoin transform
    #transform=rio.transform.from_origin(-690124, 7539362, 100, 100) # for tile 30UUG in EPSG:3857, create 100m resolutoin transform

    profile={}
    profile['dtype'] = np.int16
    profile['count'] = 1
    profile['width'] = 1098
    profile['height'] = 1098
    profile['nodata'] = -9999
    profile['transform'] = transform
    profile['driver'] = 'GTiff'
    profile['crs'] = CRS.from_epsg(32630)

    result=np.zeros((1098, 1098))

    for package in packages:

        print(package)

        gdf = gpd.read_file(package, layer='All Cloud Over Tile').to_crs(epsg=32630)

        shapes = [(g, 0) for g in list(gdf.geometry)]

        this_result = rasterize(shapes, fill=1, transform=transform, dtype=np.int16, out_shape=(1098, 1098))

        result += this_result
    
    print(result)
    print(result.shape)


    with rio.open('/Users/robwebster/gsi/seos/hadyard_hill_aoi/s2_discovery/test4.tif', 'w', **profile) as dst:
        dst.write(result.astype(np.int16), 1)

    return result



g_dir=f"/Users/robwebster/gsi/seos/hadyard_hill_aoi/s2_discovery/geopackages/subset/"

result = generate_scene_counts(g_dir)