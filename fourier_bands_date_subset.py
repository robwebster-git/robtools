import csv
import os
from posix import environ
import sys
import click
import subprocess

@click.command()
@click.argument('project')
@click.option('--tiles', default=None, type=click.STRING, help='Tile IDs, comma separated - eg 13UDA,13UDC,13UDV' )
@click.option('--years', default='2020', type=click.STRING, help="years to include - comma separated eg 2018,2019,2020")
@click.option('--months',default='1,2,3,4,5,6,7,8,9,10,11,12' , type=click.STRING, help="months to include, comma separated eg 3,4,5,6,7,8,")
@click.option('--indices', type=click.Path(exists=True), default='/home/nx06/nx06/shared/RF/S2AWS/indices_to_process.csv', help='CSV of indices to process')
def main(project, tiles, years, months, indices):

    with open(indices, newline='\n') as f:
        reader = csv.reader(f)
        indices_list = list(reader)
        indices_list = [i[0] for i in indices_list]

    print(indices_list)

    tiles = tiles.split(',')

    for tile in tiles:
        aoi = f'{project}_{tile}'

        print(f'\nProcessing {aoi}\n')

        years = years.split(',')
        months = months.split(',')

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

                band_file = f"{os.environ['MAINDIR']}/S2AWS/{aoi}/Fourier_Files_{band}_subset.txt"

                subprocess.run(f"touch {band_file}", shell=True)

                for year in years:
                    for month in months:
                        subprocess.run(f"ls {os.environ['MAINDIR']}/S2AWS/{aoi}/SAFE/*L2A_{year}{month}*/S2Fmask/R{resolution}m/*{band}*.tif >> {os.environ['MAINDIR']}/S2AWS/{aoi}/Fourier_Files_{band}_subset.txt", shell=True)



            for index in indices_list:

                index_file = f"{os.environ['MAINDIR']}/S2AWS/{aoi}/Fourier_Files_{index}_subset.txt"

                subprocess.run(f"touch {index_file}", shell=True)

                for year in years:
                    for month in months:
                        subprocess.run(f"ls {os.environ['MAINDIR']}/S2AWS/{aoi}/SAFE/*L2A_{year}{month}*/indices/*_{index}.tif >> $MAINDIR/S2AWS/{aoi}/Fourier_Files_{index}_subset.txt", shell=True)

if __name__ == "__main__":
    main()
