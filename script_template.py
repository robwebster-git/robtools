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



def some_function():
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


@click.command()
@click.option('--something', default=None, type=click.Path(exists=True), help="")
@click.option('--threshold', default=200, type=int, help="")
def main(something, threshold):
    pass


if __name__ == "__main__":
    main()