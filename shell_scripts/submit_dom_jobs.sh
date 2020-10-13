#!/bin/bash

#  Set AOI as first commandline argument

AOI=$1

#  The rest of the arguments should be S2 tile IDs eg 17SMS 17SMR 18SUE etc

for TILE in ${@:2}
do
  IN_FILES=$MAINDIR/Projects/${AOI}_${TILE}/Scores/*pc/*pc.${AOI}*.tif
  OUTFILE=$MAINDIR/Projects/${AOI}_${TILE}/Scores/dominant_${TILE}.tif
  MASKFILE=$MAINDIR/Projects/${AOI}_${TILE}/AOI/${AOI}_${TILE}.geojson
  CODEDIR=/home/nx06/nx06/rwgsi/git/ms-core/src/moonshine/validation
  #echo $IN_FILES
  #echo $OUTFILE
  #echo $MASKFILE
  #echo $CODEDIR
  python $CODEDIR/dominant-species.py --threshold 75 --species_with_fia_code --mask $MASKFILE --mask-output $OUTFILE $IN_FILES &
done
