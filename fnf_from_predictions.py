'''
Script DocString should go above the import statements

Author      : Rob Webster
Date        : October 2020       

'''

import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rio
import fiona
import os
import sys
import glob
import click



def get_fnf(size, density, winter, size_threshold, density_threshold):
    '''
    Do some task

    Fuller description of task

    Parameters
    ----------
    x : type
        Description of parameter `x`.
    y
        Description of parameter `y` (with type not specified).
    
    Returns
    -------
    z : ndarray
        Description of returned array
    a : list
        Description of ...

    
    '''

    with rio.open(size) as src:
        s = src.read(1, masked=True)
        profile = src.profile.copy()

    with rio.open(density) as src:
        d = src.read(1, masked=True)

    with rio.open(winter) as src:
        w = src.read(1, masked=True)


    zeros = np.zeros(s.shape)
    non_forest = ((s < size_threshold) & (d < density_threshold)) | (w.filled(1) == 1)

    fnf = np.ma.masked_array(zeros, mask=non_forest, fill_value=1)

    return fnf, profile


@click.command()
@click.option('--size', default=None, type=click.Path(exists=True), help="")
@click.option('--density', default=None, type=click.Path(exists=True), help="")
@click.option('--winter', default=None, type=click.Path(exists=True), help="")
@click.option('--winter_mean', default=None, type=click.Path(exists=True), help="")
@click.option('--output_directory', default="", type=click.Path(exists=True), help="")
@click.option('--tile', default="", type=str, help="")
def main(size, density, winter, winter_mean, output_directory, tile):
    
    size_threshold = 25
    density_threshold = 80

    fnf, profile = get_fnf(size, density, winter, size_threshold, density_threshold)
    profile['dtype'] = np.int8

    output_filename = f'fnf_from_predictions_and_winter_fnf_size_{size_threshold}_density{density_threshold}_{tile}.tif'
    output_path = os.path.join(output_directory, output_filename)

    with rio.open(output_path, 'w', **profile) as dst:
        dst.write(fnf.filled().astype(np.int8), 1)


if __name__ == "__main__":
    main()