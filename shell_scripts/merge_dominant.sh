#!/bin/bash

#Example for FC_CL
DIR_IN=/home/nx06/nx06/shared/RF/Projects/FC_CL/Portal/species/dominant/
DIR_OUT=/home/nx06/nx06/shared/RF/Projects/FC_CL/Portal/species/dominant/

for CRS in '17' '18'
do
  DIR_OUT_CRS=${DIR_OUT}/epsg326${CRS}
  mkdir $DIR_OUT_CRS
  FPATH_OUT_CRS=${DIR_OUT_CRS}/dominant_FC_CL.tif
  echo $CRS
  echo $FPATH_OUT_CRS
  gdal_merge.py -n -9999 -a_nodata -9999 -o ${FPATH_OUT_CRS} -co TILED=YES -co COMPRESS=DEFLATE -co NUM_THREADS=4 ${DIR_IN}/dominant_${CRS}???.tif
done
