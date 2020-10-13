#!/bin/bash

#Example for FC_CL
DIR_IN=/home/nx06/nx06/shared/RF/Projects/FC_CL/FNF/FNF2/Per_Tile
DIR_OUT=/home/nx06/nx06/shared/RF/Projects/FC_CL/FNF/FNF2/Merged

for CRS in '17' '18'
do
  DIR_OUT_CRS=${DIR_OUT}/epsg326${CRS}
  mkdir $DIR_OUT_CRS
  FPATH_OUT_CRS=${DIR_OUT_CRS}/fnf2_FC_CL.tif
  echo $CRS
  echo $FPATH_OUT_CRS
  gdal_merge.py -n -9999 -a_nodata -9999 -o ${FPATH_OUT_CRS} -co TILED=YES -co COMPRESS=DEFLATE -co NUM_THREADS=4 ${DIR_IN}/fnf2_${CRS}???.tif
done
