'''
Script DocString should go above the import statements

Author      : Rob Webster
Date        : October 2020       

'''

import geopandas as gpd
import pandas as pd
import rasterio as rio
import fiona
import os
import sys
import glob
import click



def load_raster(filename):
    '''
    Load a raster using rasterio, return the numpy array and profile

    Parameters
    ----------
    filename : str
        The path / filename of valid raster file, usually a GeoTiff.
    
    Returns
    -------
    z : ndarray
        Description of returned array
    a : list
        Description of ...

    
    '''

    with rio.open(raster) as src:
        data = src.read(1, masked=True)
        profile = src.profile.copy()

    return data, profile


def plot_raster(filename):
    '''
    Load a raster using rasterio, return the numpy array and profile

    Parameters
    ----------
    filename : str
        The path / filename of valid raster file, usually a GeoTiff.
    
    Returns
    -------
    z : ndarray
        Description of returned array
    a : list
        Description of ...

    
    '''

    with rio.open(raster) as src:
        data = src.read(1, masked=True)
        profile = src.profile.copy()

    return data, profile


@click.command()
@click.option('--something', default=None, type=click.Path(exists=True), help="")
@click.option('--threshold', default=200, type=int, help="")
def main(input_rasters, output_directory):
    
    for raster in input_rasters:

        data, profile = load_raster(raster)

        image = plot(data, profile, raster)

        image.write(output_filename, dpi=150)

if __name__ == "__main__":
    main()