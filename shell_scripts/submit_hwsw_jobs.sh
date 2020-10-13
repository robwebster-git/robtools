#!/bin/bash

#  Set AOI as first commandline argument

AOI=$1

#  The rest of the arguments should be S2 tile IDs eg 17SMS 17SMR 18SUE etc

for TILE in ${@:2}
do
  IN_FILES=$MAINDIR/Projects/SM_TEST_${TILE}/Scores/*pc/*pc_adj*.tif
  OUTFILE=$MAINDIR/Projects/SM_TEST_${TILE}/Scores/hwsw_${AOI}_${TILE}.tif
  SHAPEFILE=$MAINDIR/Projects/${AOI}_${TILE}/AOI/${AOI}_${TILE}.geojson
  CODEDIR=/home/nx06/nx06/rwgsi/git/ms-core/src/moonshine/validation/stratification
  python $CODEDIR/species_hwsw_stratification.py --from_scores $IN_FILES $OUTFILE
done
