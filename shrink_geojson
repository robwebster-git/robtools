'''
Script DocString should go above the import statements

Author      : Rob Webster
Date        : October 2020       

'''

import geopandas as gpd
import pandas as pd
import rasterio as rio
import fiona
from jsmin import jsmin
import os
import sys
import glob
import click


@click.command()
@click.option('--input', default=None, type=click.Path(exists=True), help="")
@click.option('--output', default=None, type=click.Path(), help="")
def main(input, output):
    
    with open(input, 'r+') as shapefile:
        shrunk = jsmin(shapefile.read())
        shapefile.seek(0)
        shapefile.truncate()
        shapefile.write(shrunk)
        shapefile.close()
    
    '''
    with open(output, 'w') as shrunk_shapefile:
        shrunk_shapefile.write(shrunk)
    '''
if __name__ == "__main__":
    main()