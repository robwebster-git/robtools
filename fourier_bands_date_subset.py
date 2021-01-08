import glob
import os
from posix import environ
import sys
import click
import subprocess

@click.command()
@click.argument('aoi')
@click.option('--years', default=[2020], type=list, help="years to include")
@click.option('--months',default=[1,2,3,4,5,6,7,8,9,10,11,12] , type=list, help="months to include")
def main(aoi, years, months):

    years = sorted(years)
    months = sorted(months)

    years = [str(year) for year in years]
    months = [str(month) for month in months]

    months = [month.zfill(2) for month in months]

    for resolution in [10, 20, 60]:

        bands = []

        if resolution == 10:
            bands=['B02','B03','B04','B08']
        elif resolution == 20:
            bands=['B05','B06','B07','B8A','B11','B12']
        elif resolution == 60:
            bands=['B01','B09']
    

        for band in bands:

            p = subprocess.run(f"touch {os.environ['MAINDIR']}/S2AWS/{aoi}/Fourier_Files_{band}.txt", shell=True)

            for year in years:
                for month in months:
                    p = subprocess.run(f"ls {os.environ['MAINDIR']}/S2AWS/{aoi}/SAFE/*L2A_{year}{month}*/S2Fmask/R{resolution}m/*{band}*.tif >> {os.environ['MAINDIR']}/S2AWS/{aoi}/Fourier_Files_{band}.txt", shell=True)



    



if __name__ == "__main__":
    main()
