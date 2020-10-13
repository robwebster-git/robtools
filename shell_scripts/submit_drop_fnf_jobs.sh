#!/bin/bash

#  Set AOI as first commandline argument

AOI=$1

#  The rest of the arguments should be S2 tile IDs eg 17SMS 17SMR 18SUE etc

for TILE in ${@:2}
do
  IN_FILES=$MAINDIR/Projects/SM_TEST_${TILE}/Scores/*pc/*_pc.*.tif
  CODEDIR=/home/nx06/nx06/rwgsi/git/ms-core/src/moonshine/validation  
  python $CODEDIR/drop_nonforest.py $IN_FILES
done
