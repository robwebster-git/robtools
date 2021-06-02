import rasterio as rio
import os
import glob
import numpy as np
from rsgislib.segmentation import segutils
from rsgislib.rastergis import ratutils

datadir = '/Users/robwebster/gsi/stratification_innovation/segmentation_tests'
outdir = '/Users/robwebster/gsi/stratification_innovation/segmentation_tests_outputs'
tmpath = '/Users/robwebster/gsi/stratification_innovation/tempdir'

indata = glob.glob(os.path.join(datadir, "*.tif"))
infile = os.path.join(outdir, 'input_data_normalised.tif')
 
num_bands = len(indata)

print(indata)

with rio.open(indata[0]) as src:
    profile = src.profile.copy()
    profile['dtype'] = 'float64'

profile['count'] = num_bands

bands = []

for band, tif in enumerate(indata, start=1):
    with rio.open(tif) as src:
        data = src.read(1, masked=True)
        data = data / (np.ma.max(data))
        data[data<0] = 0

        bands.append(data)

stack = np.stack(bands, axis=0)

print(stack.shape)

with rio.open(infile, 'w', **profile) as dst:
    dst.write(stack.astype('float64'))

numClusters=60
minPxls=101
distThres=200

outfile = f'{os.path.join(outdir, f"seg_numclusters_{numClusters}_minpxls_{minPxls}_distthresh_{distThres}.kea")}'

outfile_tif = outfile.replace('.kea', '.tif')

segutils.runShepherdSegmentation(infile,
                                 outfile,
                                 tmpath=tmpath,
                                 numClusters=numClusters,
                                 minPxls=minPxls,
                                 distThres=distThres)



with rio.open(outfile) as src:
    data = src.read(1)
    profile = src.profile.copy()
    profile['driver'] = 'GTiff'


print(profile)

with rio.open(outfile_tif, 'w', **profile) as dst:
    dst.write(data, 1)


bins = [0.1, 20, 40, 60, 80]

ratutils.populateImageStats(infile, outfile, calcMean=True)
band1 = ratutils.getColumnData(outfile, 'Band1Avg')
binned_band1 = np.digitize(band1, bins)
ratutils.setColumnData(outfile, 'CC', binned_band1)

print(band1, binned_band1)

