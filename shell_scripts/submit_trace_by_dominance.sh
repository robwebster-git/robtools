#!/bin/bash

#  Set AOI as first commandline argument

AOI=$1

#  The rest of the arguments should be S2 tile IDs eg 17SMS 17SMR 18SUE etc

for TILE in ${@:2}
do
  DOMFILE=$MAINDIR/Projects/${AOI}_${TILE}/Scores/dominant_${TILE}.tif
  OUTPATH=$MAINDIR/Projects/${AOI}_${TILE}/Scores/
  SHAPEFILE=$MAINDIR/Projects/${AOI}_${TILE}/AOI/${AOI}_${TILE}_32617.geojson
  CODEDIR=/home/nx06/nx06/shared/git/ms-core/src/moonshine/validation
  python $CODEDIR/trace.py --dom $DOMFILE --out $OUTPATH --aoi ${AOI}
done
